#pragma once

#include "elasticity.h"
#include "isotropic_elasticity.h"
#include "orthotropic_elasticity.h"

#include <memory>
#include <string>

#define MATERIAL_SCALAR_FIELD(name) \
    private:                                \
    Precision m_##name = -1;      \
    public:                                \
    bool has_##name () const {           \
        return m_##name != -1;    \
    }                               \
                                        \
    void set_##name (Precision value ) {   \
        m_##name = value; \
    }                               \
                                    \
    Precision get_##name () const {           \
        return m_##name ;                               \
    }                                   \
\

// Defining the fem namespace (abbreviation for finite element method)
namespace fem {
// Defining the _material namespace, used to separate the classes and functions related to _material properties
namespace material {

// Material struct defines the properties and methods related to a _material in the context of finite element analysis
    struct Material {
    private:

    // Pointer to an object of type Elasticity, used to store the elasticity properties of the _material
    ElasticityPtr m_elastic = nullptr;

    MATERIAL_SCALAR_FIELD(thermal_capacity)
    MATERIAL_SCALAR_FIELD(thermal_expansion)
    MATERIAL_SCALAR_FIELD(density)

    public:

    // set_elasticity is a templated function that takes any type T (that should inherit from Elasticity)
    // and any number of arguments to construct an object of type T, then sets the ElasticityPtr to this object
    template<typename T, typename... Args>
    void set_elasticity(Args&&... args) {
        this->m_elastic = ElasticityPtr {new T(args...)};
    }

    // Function to check if the elasticity has been set
    bool has_elasticity() const {
        return this->m_elastic != nullptr;
    }

    // Function to get the elasticity of the _material
    ElasticityPtr elasticity() const {
        return m_elastic;
    }

};

#undef MATERIAL_SCALAR_FIELD

}    // namespace _material
}    // namespace fem
