#include "nuwro_metropolis.h"
// #include <thread>

const event & get_event(const char *filename){
  thread_local NuWro_metropolis nuwro_metropolis(filename);
  return nuwro_metropolis.get_event();
}