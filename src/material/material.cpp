/**
 * @file material.cpp
 * @brief Implements the material property container.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#include "material.h"

#include "../core/logging.h"

#include <utility>

namespace fem::material {

Material::Material(std::string name)
    : Namable(std::move(name)) {}

bool Material::has_elasticity() const {
    return m_constitutive.has_elasticity();
}

bool Material::has_plasticity() const {
    return m_constitutive.has_plasticity();
}

ElasticityPtr Material::elasticity() const {
    logging::error(!has_plasticity(),
        "Direct elasticity access is invalid for plastic material '", name,
        "'; use ConstitutiveLaw for solver evaluation");
    return m_constitutive.elasticity();
}

ElasticityPtr Material::elastic_backbone() const {
    logging::error(has_elasticity(),
        "Material '", name, "' has no elastic backbone");
    return m_constitutive.elasticity();
}

PlasticityPtr Material::plasticity() const {
    return m_constitutive.plasticity();
}

void Material::set_density(Precision value) {
    logging::error(value >= Precision(0), "Material density must be non-negative");
    m_density = value;
}

void Material::info() const {
    logging::info(true, "Material: ", name);
    logging::info(true, "   Elasticity          : ", (has_elasticity() ? "YES" : "NO"));
    logging::info(true, "   Plasticity          : ", (has_plasticity() ? "YES" : "NO"));
    logging::info(true, "   Thermal Capacity    : ", (has_thermal_capacity() ? std::to_string(m_thermal_capacity) : "NO"));
    logging::info(true, "   Thermal Expansion   : ", (has_thermal_expansion() ? std::to_string(m_thermal_expansion) : "NO"));
    logging::info(true, "   Density             : ", (has_density() ? std::to_string(m_density) : "NO"));
}

} // namespace fem::material
