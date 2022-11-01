#pragma once
class event;
class NuWro_metropolis;
// expose minimal interface
class nuwro_interface {
public:
  nuwro_interface(const char *filename);
  event gen_event();
  ~nuwro_interface();
  nuwro_interface(const nuwro_interface &) = delete;
  nuwro_interface(nuwro_interface &&);

private:
  NuWro_metropolis *instance;
};
