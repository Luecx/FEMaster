// DEBUG diagnostics only: compare the old Total-Lagrangian B matrix against
// the explicit infinitesimal-strain B matrix while keeping the old stiffness
// path active. This is intended only to diagnose shell_solid_tie.

#include "element_solid.h"

namespace fem::model {

#define FEMASTER_DEBUG_B_DIAGNOSTICS(NN)                                                        \
    template<>                                                                                   \
    MapMatrix SolidElement<NN>::stiffness(Precision* buffer) {                                   \
        StaticMatrix<NN, D> reference_coords = this->node_coords_reference();                    \
        StaticMatrix<NN, D> current_coords   = this->node_coords_current();                      \
                                                                                                 \
        Index ip = 0;                                                                            \
                                                                                                 \
        std::function<StaticMatrix<D * NN, D * NN>(Precision, Precision, Precision)> func =      \
            [this, &reference_coords, &current_coords, &ip](Precision r, Precision s, Precision t) \
                -> StaticMatrix<D * NN, D * NN> {                                                \
                Precision det0;                                                                  \
                const StaticMatrix<NN, D> dN_dX =                                               \
                    this->shape_derivatives_reference(reference_coords, r, s, t, det0);          \
                const Mat3 F = this->deformation_gradient(                                      \
                    reference_coords, current_coords, r, s, t);                                 \
                                                                                                 \
                const StaticMatrix<n_strain, D * NN> B_gl =                                     \
                    this->green_lagrange_strain_displacement(dN_dX, F);                          \
                const StaticMatrix<n_strain, D * NN> B_lin =                                    \
                    this->strain_displacement(dN_dX);                                            \
                                                                                                 \
                if (ip == 0) {                                                                   \
                    const Precision coord_diff =                                                 \
                        (current_coords - reference_coords).cwiseAbs().maxCoeff();                \
                    const Precision f_diff =                                                     \
                        (F - Mat3::Identity()).cwiseAbs().maxCoeff();                             \
                    const Precision b_diff =                                                     \
                        (B_gl - B_lin).cwiseAbs().maxCoeff();                                    \
                    const Mat3 J0 = this->jacobian(reference_coords, r, s, t);                   \
                    const Precision cond_est = J0.norm() * J0.inverse().norm();                  \
                    if (coord_diff > Precision(1e-14) ||                                         \
                        f_diff     > Precision(1e-14) ||                                         \
                        b_diff     > Precision(1e-14)) {                                         \
                        logging::info(true,                                                      \
                            "DEBUG_SOLID_B elem=", this->elem_id,                               \
                            " N=", NN,                                                           \
                            " coord_diff=", coord_diff,                                          \
                            " F_minus_I=", f_diff,                                               \
                            " Bgl_minus_Blin=", b_diff,                                          \
                            " detJ0=", det0,                                                     \
                            " cond_est=", cond_est);                                             \
                    }                                                                            \
                }                                                                                \
                                                                                                 \
                const VolumeStrainGreenLagrange strain =                                        \
                    VolumeStrainGreenLagrange::from_deformation_gradient(F);                     \
                VolumeStressPK2 stress;                                                          \
                Mat6            C;                                                               \
                                                                                                 \
                const Index      state_row = this->mp_index(ip++);                               \
                const Precision* old_state =                                                     \
                    &(*this->_model_data->material_state_old)(state_row, 0);                      \
                Precision* new_state =                                                           \
                    &(*this->_model_data->material_state_new)(state_row, 0);                      \
                evaluate_material(r, s, t, strain, old_state, new_state, stress, C);             \
                                                                                                 \
                return StaticMatrix<D * NN, D * NN>(B_gl.transpose() * (C * B_gl) * det0);      \
            };                                                                                   \
                                                                                                 \
        StaticMatrix<D * NN, D * NN> stiffness =                                                \
            integration_scheme_stiffness().integrate(func);                                     \
        stiffness = Precision(0.5) * (stiffness + stiffness.transpose());                        \
                                                                                                 \
        MapMatrix mapped{buffer, D * NN, D * NN};                                               \
        mapped = stiffness;                                                                      \
        return mapped;                                                                           \
    }

FEMASTER_DEBUG_B_DIAGNOSTICS(4)
FEMASTER_DEBUG_B_DIAGNOSTICS(5)
FEMASTER_DEBUG_B_DIAGNOSTICS(6)
FEMASTER_DEBUG_B_DIAGNOSTICS(8)
FEMASTER_DEBUG_B_DIAGNOSTICS(10)
FEMASTER_DEBUG_B_DIAGNOSTICS(15)
FEMASTER_DEBUG_B_DIAGNOSTICS(20)

#undef FEMASTER_DEBUG_B_DIAGNOSTICS

} // namespace fem::model
