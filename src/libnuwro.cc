#include "libnuwro.h"
#include "nuwro.h"
#include <memory>
#include <string>
#include <vector>

extern "C" {
void shhpythiaitokay_(void);
void youcanspeaknowpythia_(void);
}

nuwro_interface::nuwro_interface(const char *filename, std::vector<std::string> args) {
  shhpythiaitokay_();
  instance = new NuWro(filename, args);
}

nuwro_interface::nuwro_interface(const char *filename) {
  shhpythiaitokay_();
  instance = new NuWro(filename);
}

nuwro_interface::~nuwro_interface() { delete instance; }

nuwro_interface::nuwro_interface(nuwro_interface &&other) {
  instance = other.instance;
  other.instance = nullptr;
}

event nuwro_interface::gen_event() { return instance->get_event(); }