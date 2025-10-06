//
// Created by Luecx on 06.10.2025.
//
#include "command_registrar.h"
#include "registry.h"

namespace fem::reader2 {

CommandRegistrar::CommandRegistrar(const Scope& scope, const std::string& name, Registry::Configure cfg) {
    Registry::instance().register_command(scope, name, std::move(cfg));
}

}