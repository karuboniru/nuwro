#ifndef _nuwro_metropolis_h_
#define _nuwro_metropolis_h_
#include "args.h"
#include "beam.h"
#include "event1.h"
#include "geomy.h"
#include "input_data.h"
#include "metropolis.h"
#include "nucleus.h"
#include "params.h"
#include "target_mixer.h"
#include <ostream>
#include <unordered_map>
#include <string>
#include <utility>

class NuWro_metropolis {
public:
  geomy *make_detector(params &p);
  void makeevent(event *e, params &p);
  void finishevent(event *e, params &p);
  // void raport(double i, double n, const char* text, int precision=1000, int
  // k=-1, bool toFile=false);
  void init(int argc, char **argv);
	void init(const char * filename);
  // void test_events(params &p);
  // void user_events(params &p);
  // void UserAction(params& p);
  void real_events(params &p);
  void kaskada_redo(string input, string output);
  void main(int argc, char **argv);
  // inline int proces() {return _procesy.choose();}
  void set(params &p);
  void refresh_target(params &p);
  // void refresh_dyn (params &p);
  void pot_report(std::ostream &, bool format);
	event get_event();
  NuWro_metropolis();
  NuWro_metropolis(const char * filename);
  ~NuWro_metropolis();

private:
  params p;
  args a;
  // chooser _procesy;
  ofstream _progress;
  geomy *_detector;
  beam *_beam;
  nucleus *_nucleus;
  target_mixer *_mixer;
  const bool dismode = true;
  input_data input;
  Metropolis<event, size_t> sampler;
  std::vector<double> channel_sampleing_weight{};
  std::vector<double> channel_weight_sum{}, channel_weight_sum_fraction{};
  std::unordered_map<int, double> channel_count_final{};
  size_t accepted_count{};
  bool accept{false};
  std::vector<int> enabled_dyns{};
	std::mt19937 gen;
};

// extern NuWro_metropolis nuwro_metropolis;
#endif
