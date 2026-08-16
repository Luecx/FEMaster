/**
 * @file material.h
 * @brief Declares material property containers for FEM analyses.
 *
 * A Material stores scalar thermal/inertial properties and one mechanical
 * ConstitutiveLaw composed from elasticity and optional plasticity.
 *
 * Legacy element formulations may access `elasticity()` only while the material
 * is purely elastic. Once Plasticity is attached, solver code must use
 * `constitutive_law()`; auxiliary state-neutral operations that deliberately need
 * the recoverable elastic component use `elastic_backbone()` explicitly.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#pragma once

#include "../data/namable.h"
#include "constitutive_law.h"

#include <memory>
#include <string>
#include <utility>

namespace fem::material {

struct Material : public Namable {
    using Ptr = std::shared_ptr<Material>;

    explicit Material(std::string name);

    bool has_elasticity() const;
    bool has_plasticity() const;

    /**
     * @brief Returns elasticity for legacy elasticity-only formulations.
     *
     * This accessor rejects materials carrying Plasticity so an unsupported
     * element cannot silently bypass the inelastic constitutive law.
     */
    ElasticityPtr elasticity() const;

    /**
     * @brief Returns the recoverable elastic component without solver semantics.
     *
     * Use only for auxiliary state-neutral quantities such as C3D8R hourglass
     * scaling or explicitly elastic sensitivities.
     */
    ElasticityPtr elastic_backbone() const;

    PlasticityPtr plasticity() const;

    ConstitutiveLaw& constitutive_law() { return m_constitutive; }
    const ConstitutiveLaw& constitutive_law() const { return m_constitutive; }

    void info() const;

    template<typename T, typename... Args>
    void set_elasticity(Args&&... args) {
        m_constitutive.set_elasticity(
            std::make_shared<T>(std::forward<Args>(args)...)
        );
    }

    template<typename T, typename... Args>
    void set_plasticity(Args&&... args) {
        m_constitutive.set_plasticity(
            std::make_shared<T>(std::forward<Args>(args)...)
        );
    }

    bool has_thermal_capacity() const { return m_thermal_capacity >= Precision(0); }
    void set_thermal_capacity(Precision value) { m_thermal_capacity = value; }
    Precision get_thermal_capacity() const { return m_thermal_capacity; }

    bool has_thermal_expansion() const { return m_thermal_expansion >= Precision(0); }
    void set_thermal_expansion(Precision value) { m_thermal_expansion = value; }
    Precision get_thermal_expansion() const { return m_thermal_expansion; }

    bool has_density() const { return m_density >= Precision(0); }
    void set_density(Precision value);
    Precision get_density() const { return m_density; }

private:
    ConstitutiveLaw m_constitutive;
    Precision m_thermal_capacity = Precision(-1);
    Precision m_thermal_expansion = Precision(-1);
    Precision m_density = Precision(-1);
};

using MaterialPtr = Material::Ptr;

} // namespace fem::material
