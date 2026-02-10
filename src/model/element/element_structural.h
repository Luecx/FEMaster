/**
 * @file element_structural.h
 * @brief Declares the abstract base class for structural FEM elements.
 *
 * Structural elements extend `ElementInterface` with routines for stiffness,
 * mass, volume, load application, and stress/strain recovery.
 *
 * @see src/model/beam/beam.h
 */

#pragma once

#include "element.h"
#include <functional>

#include "../../material/stress.h"
#include "../../material/strain.h"

namespace fem {
namespace model {

/**
 * @struct StructuralElement
 * @brief Extends `ElementInterface` with structural-specific operations.
 */
struct StructuralElement : ElementInterface {
    explicit StructuralElement(ID elem_id)
        : ElementInterface(elem_id) {}

    ~StructuralElement() override = default;

    virtual Precision volume() = 0;
    virtual MapMatrix stiffness(Precision* buffer) = 0;
    virtual MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) = 0;
    virtual MapMatrix mass(Precision* buffer) = 0;

    virtual Stresses stress(Field& displacement, std::vector<Vec3>& rst) = 0;
    virtual Strains strain(Field& displacement, std::vector<Vec3>& rst) = 0;
    virtual std::vector<Vec6> section_forces(Field& displacement) {return std::vector<Vec6>();};

    virtual void compute_stress_strain(Field& ip_stress,
                                       Field& ip_strain,
                                       Field& displacement,
                                       int ip_offset) = 0;

    virtual void compute_stress_strain_nodal(Field& displacement,
                                             Field& stress,
                                             Field& strain) = 0;

    // Field functors for volume/surface integrations (aliased to central types)
    using ScalarField = ::fem::ScalarField;
    using VecField    = ::fem::VecField;
    using TenField    = ::fem::TenField;

    // Integrate fields over the element measure (volume for solids/beams/trusses, area*thickness for shells)
    // Optionally scale by material density (per element).
    virtual Precision integrate_scalar_field(bool scale_by_density,
                                             const ScalarField& field) = 0;
    virtual Vec3      integrate_vector_field(bool scale_by_density,
                                             const VecField& field) = 0;
    virtual Mat3      integrate_tensor_field(bool scale_by_density,
                                             const TenField& field) = 0;

    // Integrate a vector field f(x) over the element volume and scatter to nodes.
    // If scale_by_density is true, multiply f(x) by the material density (per element).
    virtual void integrate_vec_field(Field& node_loads,
                                     bool scale_by_density,
                                     const VecField& field) = 0;

    virtual void apply_vload(Field& node_loads, Vec3 load) = 0;
    virtual void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) = 0;
    virtual void compute_compliance(Field& displacement, Field& result) = 0;
    virtual void compute_compliance_angle_derivative(Field& displacement, Field& result) = 0;
};

} // namespace model
} // namespace fem
