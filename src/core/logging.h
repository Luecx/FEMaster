#pragma once

#include <iostream>
#include <stdexcept>
#include <string>

// Forward declarations of the functions
template<typename... Args> inline void log_warning(bool condition, Args... args);
template<typename... Args> inline void log_info(bool condition, Args... args);
template<typename... Args> inline void log_error(bool condition, Args... args);

#include "logging.ipp"

