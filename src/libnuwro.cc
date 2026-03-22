#include "libnuwro.h"
#include "nuwro.h"
#include "dirs.h"
#include <sstream>
#include <stdexcept>

namespace libnuwro {

// ---------------------------------------------------------------------------
// Impl — owns the NuWro instance and accumulated configuration
// ---------------------------------------------------------------------------
struct Generator::Impl {
    NuWro    nuwro;
    params   p;
    std::string data_dir_override;
    bool     initialized = false;
};

// ---------------------------------------------------------------------------
// Generator — pimpl forwarding
// ---------------------------------------------------------------------------

Generator::Generator()  : impl(new Impl) {}
Generator::~Generator() = default;

Generator::Generator(Generator &&) noexcept            = default;
Generator &Generator::operator=(Generator &&) noexcept = default;

Generator &Generator::from_file(const std::string &params_file)
{
    impl->p.read(params_file.c_str());
    return *this;
}

Generator &Generator::set(const std::string &key, const std::string &value)
{
    std::istringstream ss(key + " = " + value + "\n");
    impl->p.read(ss, "Generator::set");
    return *this;
}

Generator &Generator::with_params(const params &p)
{
    impl->p = p;
    return *this;
}

Generator &Generator::with_data_dir(const std::string &path)
{
    impl->data_dir_override = path;
    return *this;
}

void Generator::initialize()
{
    if (impl->initialized)
        throw std::logic_error("libnuwro::Generator::initialize() called twice");

    // Establish data directory.  Explicit path wins; otherwise set_dirs()
    // falls back to the NUWRO environment variable (set by setup.sh).
    if (!impl->data_dir_override.empty()) {
        set_data_dir(impl->data_dir_override);
    } else {
        // No '/' in the name → set_dirs does not modify the string literal,
        // it just falls through to the PATH / NUWRO-env-var search.
        static char progname[] = "nuwro";
        set_dirs(progname);
    }

    impl->nuwro.prepare_mh(impl->p);
    // Sync back the final params (prepare_mh forces use_mh = 1).
    impl->p = impl->nuwro.get_params();
    impl->initialized = true;
}

event Generator::gen_event()
{
    if (!impl->initialized)
        throw std::logic_error("libnuwro::Generator::gen_event() called before initialize()");
    return impl->nuwro.get_event();
}

const params &Generator::get_params() const
{
    return impl->p;
}

} // namespace libnuwro
