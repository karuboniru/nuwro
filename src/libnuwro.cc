#include "libnuwro.h"
#include "nuwro_metropolis.h"
#include <memory>

extern "C" {
void shhpythiaitokay_(void);
void youcanspeaknowpythia_(void);
}

nuwro_interface::nuwro_interface(const char *filename) {
  shhpythiaitokay_();
  instance = new NuWro_metropolis(filename);
}

nuwro_interface::~nuwro_interface() { delete instance; }

nuwro_interface::nuwro_interface(nuwro_interface &&other) {
  instance = other.instance;
  other.instance = nullptr;
}

event nuwro_interface::gen_event() { return instance->get_event(); }