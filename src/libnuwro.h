#pragma once
#include <string>
#include <vector>
class event;
class NuWro;
// expose minimal interface
class nuwro_interface {
public:
  nuwro_interface(const char *filename, std::vector<std::string> args = {});
  event gen_event();
  ~nuwro_interface();
  nuwro_interface(const nuwro_interface &) = delete;
  nuwro_interface(nuwro_interface &&);

private:
  NuWro *instance;
};
