#pragma once

#include "elasticity.h"
#include "isotropic_elasticity.h"
#include "orthotropic_elasticity.h"
#include "../data/namable.h"

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

namespace fem {
namespace material {

// Material struct defines the properties and methods related to a _material in the context of finite element analysis
struct Material : public Namable {
    using Ptr = std::shared_ptr<Material>;

    public:
    Material(std::string name) : Namable(name) {}

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

    void info() {
        logging::info(true, "Material: ", name);
        logging::info(true, "   Elasticity          : ", (has_elasticity() ? "YES" : "NO"));
        logging::info(true, "   Thermal Capacity    : ", (has_thermal_capacity() ? std::to_string(m_thermal_capacity) : "NO"));
        logging::info(true, "   Thermal Expansion   : ", (has_thermal_expansion() ? std::to_string(m_thermal_expansion) : "NO"));
        logging::info(true, "   Density             : ", (has_density() ? std::to_string(m_density) : "NO"));
    }

};

#undef MATERIAL_SCALAR_FIELD

}    // namespace _material
}    // namespace fem
