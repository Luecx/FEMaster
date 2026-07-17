#pragma once

#include "shell.h"

#include "../geometry/surface/surface4.h"

namespace fem::model {

struct QSPT : ShellElement<4> {
    using NodeCoords      = StaticMatrix<4, 3>;
    using ShapeFunction   = StaticVector<4>;
    using ShapeDerivative = StaticMatrix<4, 2>;
    using MassMatrix      = StaticMatrix<12, 12>;
    using StiffnessMatrix = StaticMatrix<12, 12>;

    struct ShearState {
        StaticVector<4>         q = StaticVector<4>::Zero();
        StaticVector<12>        fn = StaticVector<12>::Zero();
        std::array<Precision, 4> edge_lengths {};
        Precision               flexibility = Precision(0);
    };

    Surface4               geometry;
    math::quadrature::Quadrature integration_scheme_;

    QSPT(ID p_elem_id, std::array<ID, 4> p_node_ids);

    std::string type_name() const override { return "QSPT"; }

    SurfacePtr surface(ID surface_id) override;
    Precision  volume() override;
    MapMatrix  stiffness(Precision* buffer) override;
    MapMatrix  stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override;
    MapMatrix  mass(Precision* buffer) override;
    const math::quadrature::Quadrature& integration_scheme() const override;

    Dim    num_ip() const override { return 0; }
    ElDofs dofs() const override { return ElDofs{true, true, true, false, false, false}; }

    void compute_internal_force_nonlinear(Field& node_forces,
                                          const Field& ip_stress,
                                          int ip_offset) override;

    bool compute_shear_flow(Field& shear_flow,
                            const Field& displacement,
                            int offset) override;

    Precision integrate_scalar_field(bool scale_by_density, const ScalarField& field) override;
    Vec3      integrate_vector_field(bool scale_by_density, const VecField& field) override;
    void      integrate_vector_field(Field& node_loads,
                                     bool scale_by_density,
                                     const VecField& field) override;
    Mat3      integrate_tensor_field(bool scale_by_density, const TenField& field) override;

    void compute_compliance(Field& displacement, Field& result) override;
    void compute_compliance_angle_derivative(Field& displacement, Field& result) override;

private:
    struct GeometryData {
        NodeCoords                coords = NodeCoords::Zero();
        std::array<Vec3, 4>       midpoints {};
        std::array<Vec3, 4>       edges {};
        std::array<Precision, 4>  edge_lengths {};
        Precision                 h12 = Precision(0);
        Precision                 h13 = Precision(0);
        Precision                 h23 = Precision(0);
        Precision                 h24 = Precision(0);
        Precision                 h34 = Precision(0);
        Precision                 h31 = Precision(0);
        Precision                 area = Precision(0);
    };

    GeometryData     geometry_data();
    ShearState       shear_state();
    StiffnessMatrix  stiffness_impl();
    MassMatrix       mass_impl();
    Precision        effective_density();
    Precision        effective_shear_modulus();
    StaticVector<12> displacement_vector(const Field& displacement);
};
} // namespace fem::model
