// DEBUG A/B only: replace the Total-Lagrangian B matrix in the old linear
// solid stiffness path by the explicit infinitesimal-strain B matrix.
//
// This file is intentionally isolated from the production implementation and
// can be deleted after the shell_solid_tie regression has been diagnosed.

#include "element_solid.h"

namespace fem::model {

#define FEMASTER_DEBUG_LINEAR_B_STIFFNESS(NN)                                                   \
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
                /* A/B change under test: use B_lin instead of B_GL(F). */                       \
                const StaticMatrix<n_strain, D * NN> B = this->strain_displacement(dN_dX);      \
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
                return StaticMatrix<D * NN, D * NN>(B.transpose() * (C * B) * det0);            \
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

FEMASTER_DEBUG_LINEAR_B_STIFFNESS(4)
FEMASTER_DEBUG_LINEAR_B_STIFFNESS(5)
FEMASTER_DEBUG_LINEAR_B_STIFFNESS(6)
FEMASTER_DEBUG_LINEAR_B_STIFFNESS(8)
FEMASTER_DEBUG_LINEAR_B_STIFFNESS(10)
FEMASTER_DEBUG_LINEAR_B_STIFFNESS(15)
FEMASTER_DEBUG_LINEAR_B_STIFFNESS(20)

#undef FEMASTER_DEBUG_LINEAR_B_STIFFNESS

} // namespace fem::model
