/******************************************************************************
* @file material.cpp
* @brief Implements the Material class for managing material properties in FEM.
*
* This file provides implementations for the Material class methods, including
* setting and querying material properties and logging material information.
*
* @see material.h
******************************************************************************/

#include "material.h"
#include "../core/logging.h"

namespace fem::material {

Material::Material(std::string name)
    : Namable(std::move(name)) {}

bool Material::has_elasticity() const {
    return m_elastic != nullptr;
}

ElasticityPtr Material::elasticity() const {
    return m_elastic;
}

void Material::info() const {
    logging::info(true, "Material: ", name);
    logging::info(true, "   Elasticity          : ", (has_elasticity() ? "YES" : "NO"));
    logging::info(true, "   Thermal Capacity    : ", (has_thermal_capacity() ? std::to_string(m_thermal_capacity) : "NO"));
    logging::info(true, "   Thermal Expansion   : ", (has_thermal_expansion() ? std::to_string(m_thermal_expansion) : "NO"));
    logging::info(true, "   Density             : ", (has_density() ? std::to_string(m_density) : "NO"));
}

} // namespace fem::material
