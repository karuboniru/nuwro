#pragma once
#include <functional>
#include <memory>
#include "generatormt.h"

template <typename T, typename ...Args> class Metropolis {
public:
  Metropolis(std::function<double(T *, Args...)> w) : get_weight(w) {};
  ~Metropolis() = default;
  Metropolis(const Metropolis &) = delete;
  Metropolis &operator=(const Metropolis &) = delete;
  Metropolis(Metropolis &&) = default;
  Metropolis &operator=(Metropolis &&) = default;

  // Specify a new state, and sample if new state should be accepted
  // This transfer the ownership of the state to the Metropolis object
  bool update_state(T *new_state, Args... args) {
    std::unique_ptr<T> candidate(new_state);
    double new_weight = get_weight(new_state, args...);
    double ratio = new_weight / current_weight;
    if (ratio > 1.0 || ratio > frandom()) {
      current_state = std::move(candidate);
      current_weight = new_weight;
      return true; // accept
    }
    return false; // reject
  }

  const T &get_state() const { return *current_state; }

private:
  std::unique_ptr<T> current_state;
  std::function<double(T *, Args...)> get_weight;
  double current_weight{};
};