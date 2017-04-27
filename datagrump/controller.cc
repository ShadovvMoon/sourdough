#include <iostream>
#include <assert.h>
#include "controller.hh"
#include "timestamp.hh"
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/beta.hpp>

using namespace std;

Filter *make_filter(uint64_t size, uint64_t def)
{
  Filter *f = new Filter;
  f->size = size;
  f->start = (uint64_t*)malloc(sizeof(uint64_t) * size);
  f->clock = f->start;
  for (uint64_t i = 0; i < f->size; i++)
    *(f->start + i) = def;
  return f;
}

void push_filter(Filter *f, uint64_t elem)
{
  *f->clock = elem;
  f->clock++;
  if (f->clock >= f->start + f->size)
    {
      f->clock = f->start;
    }
}

uint64_t max_filter(Filter *f)
{
  uint64_t max = (uint64_t)0;
  for (uint64_t i = 0; i < f->size; i++)
    {
      if (*(f->start + i) > max)
        max = *(f->start + i);
    }
  return max;
}

uint64_t min_filter(Filter *f)
{
  uint64_t min = (uint64_t)(-1);
  for (uint64_t i = 0; i < f->size; i++)
    {
      if (*(f->start + i) < min)
        min = *(f->start + i);
    }
  if (BIN_VALUE(NUM_BINS) < min)
    min = BIN_VALUE(NUM_BINS);
  return min;
}


void Controller::normalise(double *prob)
{
  double sum = 0;
  for (int i = 0; i < NUM_BINS; i++)
    {
      sum += prob[i];
    }
  assert (sum > 0);
  for (int i = 0; i < NUM_BINS; i++)
    {
      prob[i] /= sum;
      assert (!std::isnan((double)prob[i]));
    }
}


uint64_t growth_tick = 0;
double poissonpdf( const double rate, int counts )
{
  if ( rate == 0 )
    {
      return ( counts == 0 );
    }
  else
    {
      return boost::math::pdf( boost::math::poisson( rate ), counts );
    }
}

/* Default constructor */
Controller::Controller( const bool debug )
  : debug_( debug )
{
  the_window_size = 4;
  rtF = make_filter(5, (uint64_t) - 1);
  bwF = make_filter(5, (uint64_t)0);

  for (int i = 0; i < NUM_BINS; i++)
    probs[i] = 1;

  evolveTime = 0.001 * TICK_LENGTH;
  //----
  stddev = BROWNIAN_MOTION_RATE * sqrt( evolveTime );
  zero_escape_probability = 1 - poissonpdf( evolveTime * OUTAGE_ESCAPE_RATE, 0 );

  assert( zero_escape_probability >= 0 );
  assert( zero_escape_probability <= 1.0 );

}

/* Get current window size, in datagrams */
unsigned int Controller::window_size( void )
{
  if ( debug_ )
    {
      cerr << "At time " << timestamp_ms()
           << " window size is " << the_window_size << endl;
    }

  return the_window_size;
}

void Controller::window_event_tick(void)
{
  int64_t min_window = 4;
  double threshold = RTT_TARGET;
  double rtt = threshold / 1000.0; // 50ms vs avg_rtt

  /* Observations */
  if (packets && prev_timestamp)
    {
      uint64_t sent = nAck - packets;
      uint64_t now = timestamp_ms();
      uint64_t time_slice = now - prev_timestamp;
      if (time_slice >= TICK_LENGTH)
        {
          double time = 0.001 * time_slice;

          // Update likelihood
          normalise(probs);
          for (int i = 0; i < NUM_BINS; i++)
            {
              probs[i] *= poissonpdf( BIN_VALUE(i) * time, sent);
            }
          require_update = true;

          // Reset
          prev_timestamp = now;
          packets = nAck;

          if (false)
            {
              double max_prob = 0.0;
              double max_rate = 0.0;

              normalise(probs);
              for (int i = 0; i < NUM_BINS; i++)
                {
                  if (max_prob < 0.5)
                    {
                      max_prob += probs[i];
                      max_rate = BIN_VALUE(i);
                    }
                }
              cerr << timestamp_ms() << "\t" << max_rate << "\t" << sent  / time << endl;
            }
        }
    }
  else
    {
      prev_timestamp = timestamp_ms();
      packets = nAck;
    }


  /* Evolution */
  if (!prev_modestamp)
    {
      cerr << "no prevmode" << endl;
      prev_modestamp = timestamp_ms();
      return;
    }

  // Update
  double forecast = min_filter(rtF);
  uint64_t time = (timestamp_ms() - prev_modestamp);
  require_update |= (prev_modestamp + TICK_LENGTH < time);
  if (!require_update)
    return;
  require_update = false;

  // Predict the future
  normalise(probs);
  for (int i = 0; i < NUM_BINS; i++)
    future_probs[i] = probs[i];
  int64_t reset = prev_modestamp;
  double qmin_rate = 1.0;

  //
  boost::math::normal diffdist( 0, stddev );
  while ( prev_modestamp + TICK_LENGTH < time + forecast)
    {
      normalise(future_probs);
      for (int i = 0; i < NUM_BINS; i++)
        new_probs[i] = 0;

      for (int i = 0; i < NUM_BINS; i++)
        {
          double old_rate = BIN_VALUE(i);
          int f = floor((old_rate - 5 * stddev) / BIN_WIDTH);
          int m = ceil((old_rate + 5 * stddev) / BIN_WIDTH);
          for (int v = f; v <= m; v++)
            {
              if (v < 0 || v >= NUM_BINS) continue;
              double a = boost::math::cdf(diffdist, (BIN_CEIL(v) - old_rate));
              double b = boost::math::cdf(diffdist, BIN_FLOOR(v) - old_rate);
              double contribution = future_probs[i] * (a - b);

              /*
              assert( !std::isnan( zfactor ) );
              assert( !std::isnan( future_probs[i] ) );
              assert( !std::isnan( BIN_VALUE(v) ) );
              assert( !std::isnan( BIN_VALUE(i) ) );
              assert( !std::isnan(a));
              assert( !std::isnan(b));
              assert( contribution >= 0.0 );
              assert( contribution <= 1.0 );
              */

              new_probs[v] += contribution;
            }
        }

      for (int i = 0; i < NUM_BINS; i++)
        {
          future_probs[i] = new_probs[i];
          if (prev_modestamp + TICK_LENGTH < time)
            probs[i] = new_probs[i];
        }

      // Estimate
      // Determine new window size
      double max_prob_future = 0.0;
      double max_rate_future = 0.0;
      normalise(future_probs);
      for (int i = 0; i < NUM_BINS; i++)
        {
          if (max_prob_future < 1.0 - RTT_CONFIDENCE)
            {
              max_prob_future += future_probs[i];
              max_rate_future = BIN_VALUE(i);
            }
        }
      if (max_rate_future >= qmin_rate)
        {
          qmin_rate = max_rate_future;
        }

      //----
      if (prev_modestamp + TICK_LENGTH < time)
          reset += TICK_LENGTH;
      prev_modestamp += TICK_LENGTH;
    }
  prev_modestamp = reset;


  // Determine new window size
  double max_prob_future = 0.0;
  double max_rate_future = 0.0;

  double max_prob = 0.0;
  double max_rate = 0.0;

  normalise(probs);
  normalise(future_probs);
  for (int i = 0; i < NUM_BINS; i++)
    {
      if (max_prob < 1.0 - RTT_CONFIDENCE)
        {
          max_prob += probs[i];
          max_rate = BIN_VALUE(i);
        }

      if (max_prob_future < 1.0 - RTT_CONFIDENCE)
        {
          max_prob_future += future_probs[i];
          max_rate_future = BIN_VALUE(i);
        }
    }

  // How many packets are currently in-flight?
  int64_t inflight = nSent - nAck;

  // How much do we think the system can handle?
  int64_t future = max(max_rate, qmin_rate); //max_rate_future);
  double mbps = future * 1542.0 * 8.0 / 1000.0 / 1000.0;

  the_window_size = future * rtt;
  the_window_size -= floor(avg_rtt / min_filter(rtF));

  // Bounds
  if (false)
    cout << the_window_size << "\t\t" << timestamp_ms() << "\t\tRATE: " << avg_rtt << "\t" << mbps << "\t" << future << "\t" << max_rate_future << "\t" << max_rate << "\t" << inflight  << "\t" << nSent  << "\t" << nAck << endl;
  if (the_window_size < min_window)
    the_window_size = min_window;
}

/* A datagram was sent */
void Controller::datagram_was_sent( const uint64_t sequence_number,
                                    /* of the sent datagram */
                                    const uint64_t send_timestamp )
/* in milliseconds */
{
  nSent = sequence_number;
  //window_event_tick();

  //if (the_window_size > 200)
  //  the_window_size = 200;
  //the_window_size = 200;

  /* Default: take no action */
  if ( debug_ )
    {
      cerr << "At time " << send_timestamp
           << " sent datagram " << sequence_number << endl;
    }
}

/* An ack was received */
void Controller::ack_received( const uint64_t sequence_number_acked,
                               /* what sequence number was acknowledged */
                               const uint64_t send_timestamp_acked,
                               /* when the acknowledged datagram was sent (sender's clock) */
                               const uint64_t recv_timestamp_acked,
                               /* when the acknowledged datagram was received (receiver's clock)*/
                               const uint64_t timestamp_ack_received )
/* when the ack was received (by sender) */
{
  nAck = sequence_number_acked;
  const uint64_t RTT = timestamp_ack_received - send_timestamp_acked;

  double N = RTT_SMOOTH;
  avg_rtt -= avg_rtt / N;
  avg_rtt += RTT / N;
  push_filter(rtF, RTT);
  window_event_tick();

  if ( debug_ )
    {
      cerr << "At time " << timestamp_ack_received
           << " received ack for datagram " << sequence_number_acked
           << " (send @ time " << send_timestamp_acked
           << ", received @ time " << recv_timestamp_acked << " by receiver's clock)"
           << endl;
    }
}

/* How long to wait (in milliseconds) if there are no acks
   before sending one more datagram */
unsigned int Controller::timeout_ms( void )
{
  return max_filter(rtF);
}
