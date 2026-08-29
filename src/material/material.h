/**
 * @file material.h
 * @brief Declares material property containers for FEM analyses.
 *
 * A `Material` encapsulates scalar properties such as density and thermal
 * expansion and stores a polymorphic elasticity model.
 *
 * @see src/material/material.cpp
 * @see src/material/elasticity.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../data/namable.h"
#include "elasticity.h"

#include <memory>
#include <string>
#include <utility>

namespace fem {
namespace material {

/**
 * @struct Material
 * @brief Holds scalar material data and an elasticity model.
 */
struct Material : public Namable {
    using Ptr = std::shared_ptr<Material>; ///< Shared pointer alias used across the codebase.

    /**
     * @brief Constructs a material with the provided name.
     *
     * @param name Identifier of the material.
     */
    explicit Material(std::string name);

    /// Returns `true` when an elasticity model is associated with this material.
    bool has_elasticity() const;

    /// Provides access to the underlying elasticity model.
    ElasticityPtr elasticity() const;

    /// Logs material information for diagnostics.
    void info() const;

    /**
     * @brief Replaces the elasticity model with a newly constructed instance.
     *
     * @tparam T Elasticity type deriving from `Elasticity`.
     * @tparam Args Constructor argument types.
     * @param args Arguments forwarded to the elasticity constructor.
     */
    template<typename T, typename... Args>
    void set_elasticity(Args&&... args) {
        m_elastic = ElasticityPtr(new T(std::forward<Args>(args)...));
    }

    bool has_thermal_specific_heat() const { return m_thermal_specific_heat >= Precision(0); }
    bool has_thermal_conductivity() const { return m_thermal_conductivity >= Precision(0); }

    Precision get_thermal_specific_heat() const { return m_thermal_specific_heat; }
    Precision get_thermal_conductivity() const { return m_thermal_conductivity; }

    void set_thermal_specific_heat(Precision value) { m_thermal_specific_heat = value; }
    void set_thermal_conductivity(Precision value) { m_thermal_conductivity = value; }

    bool has_thermal_expansion() const { return m_thermal_expansion >= Precision(0); }
    void set_thermal_expansion(Precision value) { m_thermal_expansion = value; }
    Precision get_thermal_expansion() const { return m_thermal_expansion; }

    bool has_density() const { return m_density >= Precision(0); }
    void set_density(Precision value);
    Precision get_density() const { return m_density; }

private:
    ElasticityPtr m_elastic = nullptr; ///< Elasticity model associated with the material.
    Precision m_thermal_specific_heat   = Precision(-1); ///< Thermal capacity value.
    Precision m_thermal_conductivity    = Precision(-1); ///< Thermal conductivity value.
    Precision m_thermal_expansion       = Precision(-1); ///< Thermal expansion coefficient.
    Precision m_density                 = Precision(-1); ///< Density value.
};
} // namespace material
} // namespace fem
