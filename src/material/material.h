/******************************************************************************
* @file material.h
* @brief Defines the Material class for managing material properties in FEM analysis.
*
* The Material class provides a mechanism to define and access material properties
* such as elasticity, density, thermal capacity, and thermal expansion. The class
* also supports setting and querying these properties dynamically.
*
* @details The class leverages polymorphism for elasticity definitions and provides
* utility methods for logging material information.
*
* @see material.cpp
* @author Finn Eggers
* @date 20.12.2024
******************************************************************************/

#pragma once  // Ensures this file is only included once during compilation

#include "../data/namable.h"
#include "elasticity.h"
#include "isotropic_elasticity.h"
#include "orthotropic_elasticity.h"
#include <memory>
#include <string>

namespace fem::material {

/******************************************************************************
* @class Material
* @brief Class for defining material properties and their elasticity in FEM.
*
* The Material class is used to define material properties such as elasticity,
* density, thermal capacity, and thermal expansion. It also supports polymorphic
* definitions of elasticity via a shared pointer.
******************************************************************************/
struct Material : public Namable {
    using Ptr = std::shared_ptr<Material>; ///< Shared pointer alias for Material.

public:
    /******************************************************************************
    * @brief Constructor for the Material class.
    *
    * Initializes a Material object with the given name.
    * @param name Name of the material.
    ******************************************************************************/
    explicit Material(std::string name);

private:
    ElasticityPtr m_elastic = nullptr; ///< Pointer to an object defining material elasticity.

#define MATERIAL_SCALAR_FIELD(name)                                          \
private:                                                                     \
    Precision m_##name = -1;                                                 \
public:                                                                      \
    /** @brief Checks if the field has been set. */                          \
    bool has_##name() const { return m_##name != -1; }                       \
    /** @brief Sets the value of the field. */                               \
    void set_##name(Precision value) { m_##name = value; }                   \
    /** @brief Retrieves the value of the field. */                          \
    Precision get_##name() const { return m_##name; }

    MATERIAL_SCALAR_FIELD(thermal_capacity) ///< Field for thermal capacity.
    MATERIAL_SCALAR_FIELD(thermal_expansion) ///< Field for thermal expansion.
    MATERIAL_SCALAR_FIELD(density) ///< Field for density.
#undef MATERIAL_SCALAR_FIELD

public:
    /******************************************************************************
    * @brief Sets the elasticity model for the material.
    *
    * This templated function constructs an elasticity model of type T with the
    * provided arguments and assigns it to the material.
    * @tparam T Elasticity type (must inherit from Elasticity).
    * @param args Arguments for constructing the elasticity model.
    ******************************************************************************/
    template<typename T, typename... Args>
    void set_elasticity(Args&&... args) {
        m_elastic = ElasticityPtr(new T(std::forward<Args>(args)...));
    }

    /******************************************************************************
    * @brief Checks if an elasticity model is defined for the material.
    * @return True if elasticity is defined, false otherwise.
    ******************************************************************************/
    bool has_elasticity() const;

    /******************************************************************************
    * @brief Retrieves the elasticity model of the material.
    * @return Pointer to the elasticity model.
    ******************************************************************************/
    ElasticityPtr elasticity() const;

    /******************************************************************************
    * @brief Logs information about the material properties.
    ******************************************************************************/
    void info() const;
};

} // namespace fem::material