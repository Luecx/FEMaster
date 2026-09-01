/**
 * @file element_structural.h
 * @brief Declares the abstract base class for structural FEM elements.
 *
 * Structural elements extend `ElementInterface` with the mechanical operators
 * required by linear, prestressed and nonlinear analyses. The interface exposes
 * stiffness, geometric stiffness, nonlinear tangent/internal-force evaluation,
 * mass and structural result recovery while formulation-specific stress,
 * strain and constitutive integration remain responsibilities of the derived
 * element classes.
 *
 * The nonlinear interface is intentionally expressed in terms of nodal
 * displacement and nodal internal force. Stress and stress-resultant quantities
 * required during an element evaluation are temporary formulation data and do
 * not form part of the solver-facing interface.
 *
 * @see src/model/beam/beam.h
 *
 * @author Finn Eggers
 * @date 01.09.2026
 */

#pragma once

#include "element.h"
#include <functional>

namespace fem {
namespace model {

/**
 * @struct StructuralElement
 * @brief Common mechanical interface implemented by all structural elements.
 *
 * A structural element provides the discrete operators required by the global
 * structural analyses. Linear analyses request the material/linear stiffness
 * through `stiffness()`. Prestress-dependent analyses request the geometric
 * stiffness through `stiffness_geom()` for a supplied nodal displacement state.
 * Nonlinear equilibrium evaluations use `stiffness_tangent()`, which always
 * assembles the internal force and optionally assembles the consistent tangent.
 *
 * The nonlinear tangent routine uses the stiffness buffer itself as the request
 * for tangent assembly. A non-null buffer requests both the tangent matrix and
 * the matching internal force. A null buffer requests an internal-force-only
 * evaluation, which is useful for residual evaluations such as line-search
 * trials where no tangent matrix is required. Derived elements therefore keep a
 * single physical nonlinear evaluation path for both operations.
 *
 * Constitutive history follows a strict ownership rule. Read-only operators such
 * as `stiffness()` and `stiffness_geom()` may inspect committed material state but
 * must not modify the persistent trial state. `stiffness_tangent()` evaluates the
 * physical trial configuration from the committed state and may write the trial
 * material state. Promotion of trial state to committed state belongs to the
 * nonlinear increment controller, not to an element.
 *
 * Stress, strain, constitutive tangents and stress resultants used while building
 * element operators are local intermediate quantities. Derived formulations may
 * keep them in local fixed-size storage or other formulation-local evaluation
 * data, but they are not exchanged with the global solver as scratch fields.
 * Result recovery remains available separately through `compute_stress_strain()`
 * and the formulation-specific recovery callbacks below.
 */
struct StructuralElement : ElementInterface {
    explicit StructuralElement(ID elem_id)
        : ElementInterface(elem_id) {}

    ~StructuralElement() override = default;

    // Fundamental mechanical measures and operators. stiffness() is the linear
    // material operator. stiffness_geom() evaluates the stress-dependent
    // geometric contribution for the supplied displacement state without
    // advancing persistent material history.
    virtual Precision volume() = 0;
    virtual MapMatrix stiffness(Precision* buffer) = 0;
    virtual MapMatrix stiffness_geom(Precision* buffer, const Field& displacement) = 0;

    // Physical nonlinear equilibrium evaluation. The internal force is always
    // accumulated into nodal_forces. If buffer is non-null, the same material
    // evaluation also assembles and returns the consistent tangent matrix. If
    // buffer is null, the element skips tangent assembly and returns an empty
    // mapped matrix while still evaluating the complete internal force and trial
    // material state for the supplied displacement.
    virtual MapMatrix stiffness_tangent(Precision* buffer,
                                        NodeData&  nodal_forces,
                                        const Field& displacement) = 0;

    // Inertial operator and formulation classification
    virtual MapMatrix mass() = 0;
    virtual bool      is_shell() const { return false; }
    virtual bool      is_solid() const { return false; }

    // Natural-coordinate locations used by generic structural result recovery
    virtual RowMatrix stress_strain_nodal_rst() { return RowMatrix(0, 3); }
    virtual RowMatrix stress_strain_ip_rst() { return RowMatrix(0, 3); }

    // Some elements precompute step-local data and free it after the load case.
    virtual void step_begin() {}
    virtual void step_end()   {}

    // Field functors for volume/surface integrations (aliased to central types)
    using ScalarField = ::fem::ScalarField;
    using VecField    = ::fem::VecField;
    using TenField    = ::fem::TenField;

    // Integrate fields over the element measure (volume for solids/beams/trusses, area*thickness for shells)
    // Optionally scale by material density (per element).
    virtual Precision integrate_scalar_field(bool               scale_by_density,
                                             const ScalarField& field) = 0;
    virtual Vec3      integrate_vector_field(bool             scale_by_density,
                                             const VecField& field) = 0;
    virtual void      integrate_vector_field(Field&          node_loads,
                                             bool            scale_by_density,
                                             const VecField& field) = 0;
    virtual Mat3      integrate_tensor_field(bool             scale_by_density,
                                             const TenField& field) = 0;

    virtual void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) = 0;

    // Structural result recovery. These routines are independent of temporary
    // stresses or resultants used during operator assembly and explicitly write
    // user-visible result fields at requested natural-coordinate locations.
    virtual void compute_stress_strain(
        Field*           strain,
        Field*           stress,
        const Field&     displacement,
        const RowMatrix& rst,
        int              offset,
        bool             use_green_lagrange_nl
    ) = 0;

    virtual void compute_compliance(Field& displacement, Field& result) {
        (void) displacement;
        (void) result;
    }
    virtual void compute_compliance_angle_derivative(Field& displacement, Field& result) {
        (void) displacement;
        (void) result;
    }
    virtual bool compute_shear_flow(Field& shear_flow, const Field& displacement, int offset) {
        (void) shear_flow;
        (void) displacement;
        (void) offset;
        return false;
    }
    virtual bool compute_beam_section_forces(Field& section_forces, const Field& displacement, int offset) {
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
};
} // namespace model
} // namespace fem
