#pragma once

#include "types_eig.h"

namespace fem {

struct Config{
    ID max_threads = 1;
};

extern Config global_config;

} // namespace fem