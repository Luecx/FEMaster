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

    virtual bool compute_shear_flow(Field& shear_flow,
                                    const Field& displacement,
                                    int offset) {
        (void) shear_flow;
        (void) displacement;
        (void) offset;
        return false;
    }

    virtual RowMatrix stress_strain_nodal_rst() { return RowMatrix(0, 3); }
    virtual RowMatrix stress_strain_ip_rst() { return RowMatrix(0, 3); }

    virtual void compute_stress_strain(Field* strain,
                                       Field* stress,
                                       const Field& displacement,
                                       const RowMatrix& rst,
                                       int offset,
                                       bool use_green_lagrange_nl) {
        (void) strain;
        (void) stress;
        (void) displacement;
        (void) rst;
        (void) offset;
        (void) use_green_lagrange_nl;
        logging::error(false, "compute_stress_strain is not implemented for element ", this->elem_id);
    }

    virtual void compute_stress_state(Field& stress_state,
                                      const Field& displacement,
                                      int offset,
                                      bool use_green_lagrange_nl) {
        RowMatrix rst = stress_strain_ip_rst();
        if (rst.rows() == 0) {
            return;
        }
        compute_stress_strain(nullptr, &stress_state, displacement, rst, offset, use_green_lagrange_nl);
    }

    /**
     * @brief Adds nonlinear internal nodal forces to a global NODE x 6 field.
     *
     * The model positions are interpreted as the current coordinates for this
     * call. Implementations add to node_forces instead of overwriting it.
     */
    virtual void compute_internal_force_nonlinear(Field& node_forces,
                                                  const Field& ip_stress,
                                                  int ip_offset) = 0;

    virtual bool compute_beam_section_forces(Field& section_forces,
                                             const Field& displacement,
                                             int offset) {
        (void) section_forces;
        (void) displacement;
        (void) offset;
        return false;
    }

    virtual bool compute_shell_section_forces(Field& section_forces,
                                              Field& contribution_count,
                                              const Field& displacement) {
        (void) section_forces;
        (void) contribution_count;
        (void) displacement;
        return false;
    }

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
