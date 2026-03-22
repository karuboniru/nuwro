#pragma once
#include "event1.h"
#include "params.h"
#include <memory>
#include <string>

/// Thin library interface for NuWro Metropolis-Hastings event generation.
///
/// Usage:
///   libnuwro::Generator gen;
///   gen.from_file("params.txt")
///      .set("beam_energy", "1000")
///      .initialize();
///   for (int i = 0; i < N; i++) {
///       event e = gen.gen_event();
///       // analyse e ...
///   }
///
/// Configuration methods may be called in any combination before initialize().
/// initialize() must be called exactly once.  After that, gen_event() may be
/// called any number of times from a single thread.
///
/// Data directory: NuWro needs its physics data files (from the data/ install
/// directory).  The lookup order is:
///   1. An explicit path passed to with_data_dir().
///   2. The NUWRO environment variable (set automatically by setup.sh).
/// Sourcing the installed setup.sh is the recommended approach.
///
/// Auxiliary files: unlike the nuwro executable, this interface writes no
/// auxiliary files (.par, .xsec, totals.txt, kinematic distributions).

namespace libnuwro {

class Generator {
public:
    Generator();
    ~Generator();

    // Non-copyable, movable.
    Generator(const Generator &) = delete;
    Generator &operator=(const Generator &) = delete;
    Generator(Generator &&) noexcept;
    Generator &operator=(Generator &&) noexcept;

    /// Load parameters from a NuWro params file.
    /// May be called multiple times; later calls override earlier ones.
    Generator &from_file(const std::string &params_file);

    /// Override a single parameter by name, e.g. set("beam_energy", "1000").
    /// Uses the same "name = value" syntax as params.txt.
    Generator &set(const std::string &key, const std::string &value);

    /// Replace the entire params struct at once.
    Generator &with_params(const params &p);

    /// Set the NuWro data directory explicitly.
    /// If not called, the NUWRO environment variable is used.
    Generator &with_data_dir(const std::string &path);

    /// Set the random seed passed to frandom_init_no_save().
    /// 0 (default): seed from time(NULL).
    /// >1: use value directly as the MT19937 seed for reproducible runs.
    Generator &seed(int s);

    /// Initialise all NuWro internals.  Must be called before gen_event().
    /// Calls NuWro::prepare_mh() internally; no files are written.
    void initialize();

    /// Generate one unweighted event via the Metropolis-Hastings sampler.
    /// initialize() must have been called first.
    event gen_event();

    /// Return the params as they will be (or were) used after initialize().
    const params &get_params() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

} // namespace libnuwro
