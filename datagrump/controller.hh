#ifndef CONTROLLER_HH
#define CONTROLLER_HH

#include <cstdint>
#define MAX_ARRIVAL_RATE ((double)1000.0)
#define BROWNIAN_MOTION_RATE ((double)200)
#define OUTAGE_ESCAPE_RATE ((double)1)

#define RTT_TARGET 75.0
#define RTT_CONFIDENCE 0.50
#define RTT_SMOOTH 15

#define NUM_BINS 128
#define TICK_LENGTH 20
#define MAX_ARRIVALS_PER_TICK 30
#define NUM_TICKS 8
#define BIN_WIDTH ((double)(MAX_ARRIVAL_RATE / NUM_BINS))
#define BIN_VALUE(i) (((double)i + 0.5) * BIN_WIDTH)

#define BIN_FLOOR(i) ((i) * BIN_WIDTH)
#define BIN_CEIL(i) ((i) * BIN_WIDTH + BIN_WIDTH)

/* Congestion controller interface */
typedef enum
{
  Start,
  Drain,
  ProbeBW,
  ProbeRTT
} ControlState;

typedef struct
{
  uint64_t *start;
  uint64_t *clock;
  uint64_t size;
} Filter;

class Controller
{
private:
  bool debug_; /* Enables debugging output */

  /* Add member variables here */
  double the_window_size = 50;
  uint64_t prev_sequence_number_acked = (uint64_t) - 1;

  Filter *rtF = NULL;
  Filter *bwF = NULL;
  ControlState state = Start;

  uint64_t prev_modestamp = 0;
  uint64_t prev_timestamp = 0;
  uint64_t tick = 2.0;
  uint64_t packets = 0;
  uint64_t packet_rtt = 0;
  uint64_t MTU = 1000;
  uint64_t pRTT = 0;
  uint64_t pBW = 0;

  uint64_t nSent = 0;
  uint64_t nAck = 0;

  bool normalised = false;
  double probs[NUM_BINS];
  double new_probs[NUM_BINS];
  double future_probs[NUM_BINS];
  bool can_estimate = false;
  bool require_update = false;
  double avg_rtt = 0;

  void normalise(double *prob);
  void window_event_tick(void);


  double evolveTime = 0;
  double stddev = 0;
  double zero_escape_probability = 0;
  double cdf[NUM_BINS];



public:
  /* Public interface for the congestion controller */
  /* You can change these if you prefer, but will need to change
     the call site as well (in sender.cc) */

  /* Default constructor */
  Controller( const bool debug );

  /* Get current window size, in datagrams */
  unsigned int window_size( void );

  /* A datagram was sent */
  void datagram_was_sent( const uint64_t sequence_number,
                          const uint64_t send_timestamp );

  /* An ack was received */
  void ack_received( const uint64_t sequence_number_acked,
                     const uint64_t send_timestamp_acked,
                     const uint64_t recv_timestamp_acked,
                     const uint64_t timestamp_ack_received );

  /* How long to wait (in milliseconds) if there are no acks
     before sending one more datagram */
  unsigned int timeout_ms( void );
};

#endif
