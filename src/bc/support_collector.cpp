/**
 * @file support_collector.cpp
 * @brief Implements construction of heterogeneous support collectors.
 *
 * The collector stores mechanical and thermal supports through their common
 * interface. Type-selected equation generation remains in the header because
 * the requesting load case supplies the concrete support type as a template
 * argument.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "support_collector.h"

namespace fem::bc {

SupportCollector::SupportCollector(const std::string& name)
    : model::Collection<SupportInterface::Ptr>(name) {}

} // namespace fem::bc
