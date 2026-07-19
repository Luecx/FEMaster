/**
 * @file c3d8r.h
 * @brief Declares the reduced-integration C3D8 solid with hourglass control.
 */

#pragma once

#include "c3d8.h"

namespace fem::model {

class C3D8R final : public C3D8 {
public:
    C3D8R(ID elem_id, const std::array<ID, 8>& node_ids);

    std::string type_name() const override { return "C3D8R"; }

    const math::quadrature::Quadrature& integration_scheme() const override;

    RowMatrix stress_strain_nodal_rst() override;

    MapMatrix stiffness(Precision* buffer) override;

    void compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) override;

private:
    using HourglassModes = StaticMatrix<8, 4>;
    using Matrix24       = StaticMatrix<24, 24>;
    using Vector24       = StaticVector<24>;

    Matrix24 hourglass_stiffness();
};

} // namespace fem::model
