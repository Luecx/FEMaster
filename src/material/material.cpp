/**
 * @file material.cpp
 * @brief Implements the material property container.
 *
 * Provides storage for scalar properties and the associated elasticity model,
 * along with diagnostic logging helpers.
 *
 * @see src/material/material.h
 * @see src/material/elasticity.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "material.h"

#include "../core/logging.h"

#include <utility>

namespace fem {
namespace material {

/**
 * @copydoc Material::Material
 */
Material::Material(std::string name)
    : Namable(std::move(name)) {}

/**
 * @copydoc Material::has_elasticity
 */
bool Material::has_elasticity() const {
    return m_elastic != nullptr;
}

/**
 * @copydoc Material::elasticity
 */
ElasticityPtr Material::elasticity() const {
    return m_elastic;
}

/**
 * @copydoc Material::info
 */
void Material::info() const {
    logging::info(true, "Material: ", name);
    logging::info(true, "   Elasticity          : ", (has_elasticity() ? "YES" : "NO"));
    logging::info(true, "   Thermal Capacity    : ", (has_thermal_capacity() ? std::to_string(m_thermal_capacity) : "NO"));
    logging::info(true, "   Thermal Expansion   : ", (has_thermal_expansion() ? std::to_string(m_thermal_expansion) : "NO"));
    logging::info(true, "   Density             : ", (has_density() ? std::to_string(m_density) : "NO"));
}

} // namespace material
} // namespace fem
