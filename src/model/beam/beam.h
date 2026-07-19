/**
 * @file beam.h
 * @brief Declares the common templated base for structural beam elements.
 *
 * `BeamElement<N>` stores fixed-size beam connectivity and provides the shared
 * interface used by concrete three-dimensional beam formulations. The class
 * centralizes access to beam sections, profiles and isotropic material data,
 * constructs the local section frames, handles rigid section offsets and
 * adapts element-specific stiffness and mass matrices to the generic
 * `StructuralElement` interface.
 *
 * The template definitions are implemented in `beam.inl`, which is included at
 * the end of this header so that every specialization remains available at the
 * point of instantiation.
 *
 * @tparam N Number of nodes in the beam element.
 *
 * @see StructuralElement
 * @see BeamSection
 * @see beam.inl
 * @see b33.h
 */

#pragma once

#include "../../core/core.h"
#include "../../material/isotropic_elasticity.h"
#include "../../material/stress/beam_stress_resultants.h"
#include "../../section/section_beam.h"
#include "../element/element_structural.h"

#include <array>
#include <cmath>

namespace fem {
namespace model {

/**
 * @brief Common base class for structural beam elements with fixed connectivity.
 *
 * `BeamElement<N>` derives from `StructuralElement` and implements the geometry,
 * section access, coordinate transformations and generic result interfaces
 * shared by concrete beam formulations. Derived element types provide their
 * formulation-specific elastic stiffness, geometric stiffness and mass
 * matrices through the three formulation-specific implementation hooks exposed here.
 *
 * Each beam node carries three translational and three rotational degrees of
 * freedom. Local element matrices are expressed either in the original section
 * frame or in the principal bending frame and are transformed to the global
 * coordinate system before assembly.
 *
 * @tparam N Number of nodes in the beam element.
 */
template<Index N>
struct BeamElement : StructuralElement {
    // Global identifiers of the beam nodes in element-local interpolation
    // order. Geometry evaluation, degree-of-freedom gathering and load
    // distribution consistently use this fixed-size connectivity.
    std::array<ID, N> node_ids{};

    // Optional global node identifier used to define the section n1 direction.
    // A negative value selects the orientation vector stored in the assigned
    // beam section instead.
    ID orientation_node_id_ = static_cast<ID>(-1);

    // Construction from the element identifier, ordered connectivity and an
    // optional orientation node. The structural base stores the element ID;
    // this class retains the topology and orientation definition.
    BeamElement(
        ID                elem_id,
        std::array<ID, N> node_ids_in,
        ID                orientation_node_id = static_cast<ID>(-1)
    );

    // Polymorphic destruction through the structural-element interface.
    ~BeamElement() override;

    // Access to the assigned beam section, its geometric profile, its material
    // and the isotropic elastic law required by the current beam formulations.
    // Each accessor validates the corresponding model assignment before
    // returning the typed object.
    BeamSection*                   get_section();
    Profile*                       get_profile();
    material::MaterialPtr          get_material();
    material::IsotropicElasticity* get_elasticity();

    // Orientation metadata and normalized section n1 direction. An explicit
    // orientation node takes precedence over the direction stored in the beam
    // section; at least one valid non-zero definition must be available.
    ID   orientation_node() const;
    bool has_orientation_node() const;
    Vec3 orientation_direction();

    // Geometric measures evaluated from the current model coordinates. Length
    // is accumulated over consecutive beam nodes and volume is obtained from
    // this centerline length and the assigned section area.
    Precision length();
    Precision volume() override;

    // Local section frames. The base frame follows the beam axis and supplied
    // n1 direction. The principal frame additionally rotates its local y-z
    // basis according to the section inertia tensor. `rotation_matrix()` keeps
    // the established internal alias for the principal frame.
    Mat3      base_rotation_matrix();
    Precision principal_angle();
    Mat3      principal_rotation_matrix();
    Mat3      rotation_matrix();

    // Rigid-offset operators for translating beam kinematics between parallel
    // reference axes in the local section frame. The single-node operator acts
    // on [u, theta], while the element operator repeats it for all N nodes.
    static StaticMatrix<6, 6>         rigid_offset_6(Precision dy, Precision dz);
    static StaticMatrix<N * 6, N * 6> rigid_offset_N(Precision dy, Precision dz);

    // Rotate a local y-z offset from the base section frame into the principal
    // bending frame using the same angle employed by the section rotation.
    static void rotate_yz_to_principal(Precision phi, Precision& y, Precision& z);

    // Formulation-specific local-to-global element operators implemented by
    // concrete beam types. All matrices contain six degrees of freedom per
    // node and are returned in the global element ordering expected by the
    // structural assembly interface.
    virtual StaticMatrix<N * 6, N * 6> stiffness_impl() = 0;
    virtual StaticMatrix<N * 6, N * 6> stiffness_geom_impl(const Field& ip_stress, int offset) = 0;
    virtual StaticMatrix<N * 6, N * 6> mass_impl() = 0;

    // Block-diagonal degree-of-freedom transformations from global coordinates
    // to the principal section frame and to the original base section frame.
    // Each node receives identical 3x3 rotations for translations and rotations.
    StaticMatrix<N * 6, N * 6> transformation();
    StaticMatrix<N * 6, N * 6> transformation_base();

    // Generic matrix assembly adapters. The returned maps reference the caller-
    // supplied storage and are populated from the formulation-specific matrix
    // implementations above.
    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override;
    MapMatrix mass(Precision* buffer) override;

    // Nonlinear internal-force and stress-state interfaces. The base class does
    // not provide a nonlinear force formulation; stress-state evaluation uses
    // the integration points and recovery implementation of the concrete beam.
    void compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) override;
    void compute_stress_state(
        Field&          stress_state,
        const Field&    displacement,
        int             offset,
        bool            use_green_lagrange_nl
    ) override;

    // Optional structural-element operations that currently have no generic
    // beam implementation. They deliberately preserve the neutral behavior of
    // the previous inline definitions.
    void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) override;
    void compute_compliance(Field& displacement, Field& result) override;
    void compute_compliance_angle_derivative(Field& displacement, Field& result) override;
    bool compute_shear_flow(Field& shear_flow, const Field& displacement, int offset) override;

    // Fixed topology and degree-of-freedom metadata required by the common
    // element interface, together with direct connectivity access. Beam
    // elements expose no surface representation through this base class.
    ElDofs     dofs() const override;
    Dim        dimensions() const override;
    Dim        n_nodes() const override;
    Dim        num_ip() const override;
    const ID*  nodes() const override;
    SurfacePtr surface(ID surface_id) override;

    // Recover nodal beam resultants in the original section frame. Shell force
    // recovery remains unsupported for beam elements and returns false.
    bool compute_beam_section_forces(
        Field&       section_forces,
        const Field& displacement,
        int          offset
    ) override;
    bool compute_shell_section_forces(
        Field&       section_forces,
        Field&       contribution_count,
        const Field& displacement
    ) override;

    // One-point midpoint integration of scalar, vector and tensor fields over
    // the beam volume. Density scaling is optional; the nodal vector overload
    // distributes the integrated resultant equally across all beam nodes.
    Precision integrate_scalar_field(
        bool scale_by_density,
        const ScalarField& field) override;
    Vec3 integrate_vector_field(
        bool scale_by_density,
        const VecField& field) override;
    void integrate_vector_field(
        Field& node_loads,
        bool scale_by_density,
        const VecField& field) override;
    Mat3 integrate_tensor_field(
        bool scale_by_density,
        const TenField& field) override;
};

} // namespace model
} // namespace fem

#include "beam.inl"
