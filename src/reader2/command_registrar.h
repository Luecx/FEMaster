#pragma once

#include "registry.h"

namespace fem::reader2 {

/**
 * \brief Helper that registers a command during static initialization.
 */
struct CommandRegistrar {
    /// \brief Register the provided handler immediately.
    CommandRegistrar(const Scope& scope, const std::string& name, Registry::Configure cfg);
};

} // namespace fem::reader2
