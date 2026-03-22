#pragma once
#include <functional>
#include <memory>
#include "generatormt.h"  // frandom() — NuWro's global MT19937 engine

template <typename T, typename ...Args> class Metropolis {
public:
  Metropolis(std::function<double(T *, Args...)> w) : get_weight(w) {};
  ~Metropolis(){};

  // Specify a new state, and sample if new state should be accepted
  // This transfer the ownership of the state to the Metropolis object
  bool update_state(T *new_state, Args... args) {
    double new_weight = get_weight(new_state, args...);
    double ratio = new_weight / current_weight;
    if (ratio > 1.0 || ratio > frandom()) {
      current_state = std::unique_ptr<T>(new_state);
      current_weight = new_weight;
      return true; // accept
    } else {
      delete new_state;
      return false; // reject
    }
  }

  const T &get_state() const { return *current_state; }

private:
  std::unique_ptr<T> current_state;
  std::function<double(T *, Args...)> get_weight;
  double current_weight{};
};
