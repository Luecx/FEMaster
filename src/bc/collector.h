/**
 * @file collector.h
 * @brief Defines reusable named condition collectors.
 *
 * Loads, supports and thermal conditions share the same storage: an ordered
 * sequence of shared Condition pointers. The separate names identify the
 * intended input domain, while ModelData keeps their registries independent.
 * Collector selection and condition history belong to the reader.
 *
 * @see Condition
 * @see ConditionManager
 * @see model::Collection
 */

#pragma once

#include "condition.h"
#include "../data/collection.h"

namespace fem::bc {

// Named definition storage; physical meaning is selected by the owning registry.
using LoadCollector    = model::Collection<Condition::Ptr>;
using SupportCollector = model::Collection<Condition::Ptr>;
using ThermalCollector = model::Collection<Condition::Ptr>;

} // namespace fem::bc
