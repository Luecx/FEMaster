#pragma once

#include "elasticity.h"
#include "isotropic_elasticity.h"

#include <memory>
#include <string>

// Defining the fem namespace (abbreviation for finite element method)
namespace fem {
// Defining the material namespace, used to separate the classes and functions related to material properties
namespace material {

// Material struct defines the properties and methods related to a material in the context of finite element analysis
struct Material {
    private:

    // Pointer to an object of type Elasticity, used to store the elasticity properties of the material
    ElasticityPtr m_elastic = nullptr;

    // Density of the material, default value is -1 indicating it's not set
    Precision m_density = -1;

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

    // Function to get the elasticity of the material
    ElasticityPtr elasticity() const {
        return m_elastic;
    }

    // sets the density of this object. if called without any arguments, it unsets the density
    void set_density(Precision density){
        m_density = density;
    }

    // Function to check if the density has been set (value is greater than 0)
    bool has_density() const {
        return m_density > 0;
    }

    // Function to get the density of the material
    Precision density() const {
        return m_density;
    }
};
}    // namespace material
}    // namespace fem
