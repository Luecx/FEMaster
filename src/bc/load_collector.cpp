/**
 * @file load_collector.cpp
 * @brief Implements load collector construction.
 */

#include "load_collector.h"

namespace fem::bc {

LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Load::Ptr>(name) {}

} // namespace fem::bc
