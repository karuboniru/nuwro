#pragma once
#include <string>
#include <vector>
class event;
class NuWro;
#define LIBNUWRO_SUPPORT_EXTRA_ARGS 1
// expose minimal interface
class nuwro_interface {
public:
  nuwro_interface(const char *filename, std::vector<std::string> args);
  nuwro_interface(const char *filename);
  event gen_event();
  ~nuwro_interface();
  nuwro_interface(const nuwro_interface &) = delete;
  nuwro_interface(nuwro_interface &&);

private:
  NuWro *instance;
};
