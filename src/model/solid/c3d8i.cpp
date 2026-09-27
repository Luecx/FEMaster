/**
 * @file c3d8i.cpp
 * @brief Nonlinear incompatible-mode C3D8I formulation.
 */

#include "c3d8i.h"

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <utility>

namespace fem::model {

namespace {

constexpr Precision local_tolerance = Precision(1e-10);
constexpr Index local_max_iterations = 30;
constexpr Index local_max_line_search = 12;
} // namespace

C3D8I::C3D8I(ID elem_id, const std::array<ID, N>& node_ids)
    : C3D8(elem_id, node_ids) {
    for (auto& F : committed_F_) F.setIdentity();
    for (auto& F : trial_F_) F.setIdentity();
}

Vec6 C3D8I::strain_variation(const Mat3& F, const Mat3& dF) {
    const Mat3 A = F.transpose() * dF;

    Vec6 result;
    result << A(0, 0), A(1, 1), A(2, 2),
              A(1, 2) + A(2, 1),
              A(2, 0) + A(0, 2),
              A(0, 1) + A(1, 0);
    return result;
}

Vec6 C3D8I::linearized_strain(const Mat3& gradient) {
    Vec6 result;
    result << gradient(0, 0), gradient(1, 1), gradient(2, 2),
              gradient(1, 2) + gradient(2, 1),
              gradient(2, 0) + gradient(0, 2),
              gradient(0, 1) + gradient(1, 0);
    return result;
}

Precision C3D8I::stress_contraction(const Mat3& stress, const Mat3& variation) {
    return (stress.array() * variation.array()).sum();
}

void C3D8I::step_begin() {
    committed_coords_ = node_coords_current();
    trial_coords_     = committed_coords_;
    trial_internal_.setZero();
    trial_valid_ = false;

    for (auto& F : committed_F_) F.setIdentity();
    for (auto& F : trial_F_) F.setIdentity();

    Precision reference_volume = Precision(0);
    const auto& scheme = integration_scheme_stiffness();

    for (Index q = 0; q < scheme.count(); ++q) {
        const auto point = scheme.get_point(q);
        const Precision det0 = jacobian(
            committed_coords_, point.r, point.s, point.t).determinant();

        logging::error(std::isfinite(det0) && det0 > Precision(0),
            "C3D8I: invalid reference determinant in element ", elem_id);
        reference_volume += det0 * point.w;
    }

    logging::error(std::isfinite(reference_volume) && reference_volume > Precision(0),
        "C3D8I: invalid reference volume in element ", elem_id);

    characteristic_length_ = std::cbrt(reference_volume);
    nonlinear_state_initialized_ = true;
}

void C3D8I::step_end() {
    nonlinear_state_initialized_ = false;
    trial_internal_.setZero();
    trial_valid_ = false;
}

void C3D8I::nonlinear_begin_increment() {
    logging::error(nonlinear_state_initialized_,
        "C3D8I: nonlinear state is not initialized in element ", elem_id);

    trial_coords_   = committed_coords_;
    trial_F_        = committed_F_;
    trial_internal_.setZero();
    trial_valid_ = false;
}

void C3D8I::nonlinear_commit_increment() {
    logging::error(nonlinear_state_initialized_,
        "C3D8I: nonlinear state is not initialized in element ", elem_id);

    logging::error(trial_valid_,
        "C3D8I: committing increment without a valid local state in element ", elem_id);

    committed_coords_ = trial_coords_;
    committed_F_      = trial_F_;
    trial_internal_.setZero();
    trial_valid_ = false;
}

void C3D8I::nonlinear_rollback_increment() {
    if (!nonlinear_state_initialized_) return;

    trial_coords_ = committed_coords_;
    trial_F_      = committed_F_;
    trial_internal_.setZero();
    trial_valid_ = false;
}

StaticMatrix<6, C3D8I::internal_dofs>
C3D8I::linear_internal_B(
    Precision r,
    Precision s,
    Precision t,
    const StaticMatrix<N, D>& reference_coords,
    Precision& det0
) {
    const Mat3 J      = jacobian(reference_coords, r, s, t);
    const Mat3 J0     = jacobian(reference_coords, Precision(0), Precision(0), Precision(0));
    const Precision j = J.determinant();
    const Precision j0 = J0.determinant();

    logging::error(std::isfinite(j) && j > Precision(0)
                && std::isfinite(j0) && j0 > Precision(0),
        "C3D8I: invalid reference mapping in element ", elem_id);

    det0 = j;
    const Precision ratio = j0 / j;
    const Mat3 transform = J0.inverse().transpose();
    const std::array<Precision, 3> xi {r, s, t};
    const std::array<Precision, 4> h {s * t, r * t, r * s, r * s * t};

    StaticMatrix<6, internal_dofs> B =
        StaticMatrix<6, internal_dofs>::Zero();

    // Wilson modes: d[0.5*(xi_i^2-1)]/d xi_i = xi_i.
    for (Index natural = 0; natural < 3; ++natural) {
        for (Index component = 0; component < 3; ++component) {
            Mat3 parametric = Mat3::Zero();
            parametric(component, natural) =
                characteristic_length_ * xi[static_cast<std::size_t>(natural)];

            const Mat3 gradient = ratio * parametric * transform;
            B.col(3 * natural + component) = linearized_strain(gradient);
        }
    }

    // Four scalar dilatational modes:
    // [d(rst)/dr, d(rst)/ds, d(rst)/dt, rst] I.
    for (Index mode = 0; mode < volumetric_modes; ++mode) {
        const Mat3 gradient =
            ratio * h[static_cast<std::size_t>(mode)] * Mat3::Identity();
        B.col(principal_modes + mode) = linearized_strain(gradient);
    }

    return B;
}

C3D8I::Matrix37 C3D8I::linear_full_stiffness() {
    logging::error(nonlinear_state_initialized_,
        "C3D8I: step_begin() must be called before stiffness evaluation");

    const auto reference_coords = node_coords_reference();
    const auto& scheme = integration_scheme_stiffness();
    const VolumeStrainLinearized zero_strain;

    Matrix37 K = Matrix37::Zero();

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        Precision det0 = Precision(0);
        const auto dN_dX = shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const StaticMatrix<6, ndof> Bu = strain_displacement(dN_dX);

        Precision incompatible_det = Precision(0);
        const auto Bi = linear_internal_B(
            point.r, point.s, point.t, reference_coords, incompatible_det);

        logging::error(std::abs(det0 - incompatible_det)
                       <= Precision(1e-10) * std::max(Precision(1), std::abs(det0)),
            "C3D8I: inconsistent reference determinant in element ", elem_id);

        BMatrix B = BMatrix::Zero();
        B.template leftCols<ndof>() = Bu;
        B.template rightCols<internal_dofs>() = Bi;

        const Index state_row = mp_index(ip);
        const Precision* old_state =
            &(*_model_data->material_state_old)(state_row, 0);

        VolumeStressCauchy stress;
        Mat6 C;
        evaluate_material(
            point.r, point.s, point.t,
            zero_strain, old_state, nullptr, stress, C);

        K.noalias() += point.w * det0 * B.transpose() * C * B;
    }

    return K;
}

C3D8I::Vector13 C3D8I::linear_internal_parameters(const Field& displacement) {
    const Matrix37 K = linear_full_stiffness();
    const Matrix13 Kii =
        K.template block<internal_dofs, internal_dofs>(ndof, ndof);
    const StaticMatrix<internal_dofs, ndof> Kiu =
        K.template block<internal_dofs, ndof>(ndof, 0);

    Vector24 u = Vector24::Zero();
    const auto local = nodal_data<D>(displacement);
    for (Index a = 0; a < N; ++a) {
        for (Dim d = 0; d < D; ++d) {
            u(D * a + d) = local(a, d);
        }
    }

    Eigen::FullPivLU<Matrix13> lu(Kii);
    logging::error(lu.isInvertible(),
        "C3D8I: singular incompatible-mode block in element ", elem_id);

    return -lu.solve(Kiu * u);
}

MapMatrix C3D8I::stiffness(Precision* buffer) {
    const Matrix37 K = linear_full_stiffness();

    const Matrix24 Kuu = K.template block<ndof, ndof>(0, 0);
    const StaticMatrix<ndof, internal_dofs> Kui =
        K.template block<ndof, internal_dofs>(0, ndof);
    const StaticMatrix<internal_dofs, ndof> Kiu =
        K.template block<internal_dofs, ndof>(ndof, 0);
    const Matrix13 Kii =
        K.template block<internal_dofs, internal_dofs>(ndof, ndof);

    Eigen::FullPivLU<Matrix13> lu(Kii);
    logging::error(lu.isInvertible(),
        "C3D8I: singular incompatible-mode stiffness in element ", elem_id);

    const Matrix24 condensed = Kuu - Kui * lu.solve(Kiu);
    logging::error(condensed.allFinite(),
        "C3D8I: non-finite condensed stiffness in element ", elem_id);

    MapMatrix mapped{buffer, ndof, ndof};
    mapped = Precision(0.5) * (condensed + condensed.transpose());
    return mapped;
}

MapMatrix C3D8I::stiffness_geom(
    Precision*   buffer,
    const Field& displacement
) {
    const auto reference_coords = node_coords_reference();
    const Vector13 internal = linear_internal_parameters(displacement);

    Vector24 u = Vector24::Zero();
    const auto local = nodal_data<D>(displacement);
    for (Index a = 0; a < N; ++a) {
        for (Dim d = 0; d < D; ++d) {
            u(D * a + d) = local(a, d);
        }
    }

    Matrix24 geometric = Matrix24::Zero();
    const auto& scheme = integration_scheme_stiffness();

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        Precision det0 = Precision(0);
        const auto dN_dX = shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const auto Bu = strain_displacement(dN_dX);

        Precision incompatible_det = Precision(0);
        const auto Bi = linear_internal_B(
            point.r, point.s, point.t,
            reference_coords, incompatible_det);

        const VolumeStrainLinearized strain(Bu * u + Bi * internal);
        const Index state_row = mp_index(ip);
        const Precision* old_state =
            &(*_model_data->material_state_old)(state_row, 0);

        VolumeStressCauchy stress;
        Mat6 C;
        evaluate_material(
            point.r, point.s, point.t,
            strain, old_state, nullptr, stress, C);

        const Mat3 sigma = stress.tensor();
        const Precision measure = point.w * det0;

        for (Index a = 0; a < N; ++a) {
            const Vec3 dNa = dN_dX.row(a).transpose();
            for (Index b = 0; b < N; ++b) {
                const Vec3 dNb = dN_dX.row(b).transpose();
                const Precision value = dNa.dot(sigma * dNb) * measure;
                for (Dim d = 0; d < D; ++d) {
                    geometric(D * a + d, D * b + d) += value;
                }
            }
        }
    }

    MapMatrix mapped{buffer, ndof, ndof};
    mapped = Precision(0.5) * (geometric + geometric.transpose());
    return mapped;
}

C3D8I::NonlinearEvaluation C3D8I::evaluate_nonlinear(
    const Field&    displacement,
    const Vector13& internal,
    bool            with_tangent,
    bool            write_material_state
) {
    logging::error(nonlinear_state_initialized_,
        "C3D8I: nonlinear state is not initialized in element ", elem_id);

    NonlinearEvaluation result;
    const auto reference_coords = node_coords_reference();
    const auto local_displacement = nodal_data<D>(displacement);
    const StaticMatrix<N, D> current_coords =
        reference_coords + local_displacement;

    const Mat3 J_center = jacobian(
        committed_coords_, Precision(0), Precision(0), Precision(0));
    const Precision j_center = J_center.determinant();

    logging::error(std::isfinite(j_center) && j_center > Precision(0),
        "C3D8I: non-positive start-of-increment center Jacobian in element ", elem_id);

    const Mat3 transform = J_center.inverse().transpose();
    const auto& scheme = integration_scheme_stiffness();

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        const std::array<Precision, 3> xi {point.r, point.s, point.t};
        const std::array<Precision, 4> h {
            point.s * point.t,
            point.r * point.t,
            point.r * point.s,
            point.r * point.s * point.t
        };

        Precision det0 = Precision(0);
        shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);

        Precision det_start = Precision(0);
        const auto dN_dx_start = shape_derivatives_reference(
            committed_coords_, point.r, point.s, point.t, det_start);

        const Mat3 J_current = jacobian(
            current_coords, point.r, point.s, point.t);
        const Precision det_current = J_current.determinant();

        logging::error(std::isfinite(det_start) && det_start > Precision(0),
            "C3D8I: invalid start-of-increment Jacobian in element ", elem_id);
        logging::error(std::isfinite(det_current) && det_current > Precision(0),
            "C3D8I: non-positive current Jacobian in element ", elem_id,
            "\ndet: ", det_current);

        const Precision ratio = j_center / det_start;
        const Mat3 delta_F =
            J_current.transpose() *
            jacobian(committed_coords_, point.r, point.s, point.t)
                .inverse().transpose();

        std::array<Mat3, principal_modes> principal {};
        for (Index natural = 0; natural < 3; ++natural) {
            for (Index component = 0; component < 3; ++component) {
                Mat3 parametric = Mat3::Zero();
                parametric(component, natural) =
                    characteristic_length_ * xi[static_cast<std::size_t>(natural)];
                principal[static_cast<std::size_t>(3 * natural + component)] =
                    ratio * parametric * transform;
            }
        }

        std::array<Precision, volumetric_modes> theta {};
        for (Index mode = 0; mode < volumetric_modes; ++mode) {
            theta[static_cast<std::size_t>(mode)] =
                ratio * h[static_cast<std::size_t>(mode)];
        }

        Mat3 principal_F = delta_F;
        for (Index mode = 0; mode < principal_modes; ++mode) {
            principal_F.noalias() +=
                internal(mode) * principal[static_cast<std::size_t>(mode)];
        }

        Precision phi = Precision(0);
        for (Index mode = 0; mode < volumetric_modes; ++mode) {
            phi += internal(principal_modes + mode)
                 * theta[static_cast<std::size_t>(mode)];
        }

        logging::error(std::isfinite(phi) && std::abs(phi) < Precision(50),
            "C3D8I: incompatible volumetric mode diverged in element ", elem_id);

        const Precision volumetric_scale = std::exp(phi);
        const Mat3 delta_F_bar = volumetric_scale * principal_F;
        const Mat3 F_bar =
            delta_F_bar * committed_F_[static_cast<std::size_t>(ip)];

        const Precision J_bar = F_bar.determinant();
        logging::error(std::isfinite(J_bar) && J_bar > Precision(0),
            "C3D8I: non-positive enhanced deformation gradient in element ", elem_id,
            "\ndet(Fbar): ", J_bar);

        result.enhanced_F[static_cast<std::size_t>(ip)] = F_bar;

        std::array<Mat3, total_dofs> derivatives {};

        // Nodal displacement derivatives of the incremental compatible gradient.
        for (Index a = 0; a < N; ++a) {
            const Vec3 gradient = dN_dx_start.row(a).transpose();

            for (Dim p = 0; p < D; ++p) {
                Mat3 d_delta_F = Mat3::Zero();
                d_delta_F.row(p) = gradient.transpose();

                derivatives[static_cast<std::size_t>(D * a + p)] =
                    volumetric_scale * d_delta_F
                    * committed_F_[static_cast<std::size_t>(ip)];
            }
        }

        for (Index mode = 0; mode < principal_modes; ++mode) {
            derivatives[static_cast<std::size_t>(ndof + mode)] =
                volumetric_scale
                * principal[static_cast<std::size_t>(mode)]
                * committed_F_[static_cast<std::size_t>(ip)];
        }

        for (Index mode = 0; mode < volumetric_modes; ++mode) {
            derivatives[static_cast<std::size_t>(
                ndof + principal_modes + mode)] =
                theta[static_cast<std::size_t>(mode)] * F_bar;
        }

        BMatrix B = BMatrix::Zero();
        for (Index a = 0; a < total_dofs; ++a) {
            B.col(a) = strain_variation(
                F_bar, derivatives[static_cast<std::size_t>(a)]);
        }

        const VolumeStrainGreenLagrange strain =
            VolumeStrainGreenLagrange::from_deformation_gradient(F_bar);

        const Index state_row = mp_index(ip);
        const Precision* old_state =
            &(*_model_data->material_state_old)(state_row, 0);
        Precision* new_state = write_material_state
            ? &(*_model_data->material_state_new)(state_row, 0)
            : nullptr;

        VolumeStressPK2 stress;
        Mat6 C;
        evaluate_material(
            point.r, point.s, point.t,
            strain, old_state, new_state, stress,
            with_tangent ? &C : nullptr);

        const Precision measure = point.w * det0;
        result.residual.noalias() +=
            measure * B.transpose() * stress.voigt();

        if (!with_tangent) continue;

        result.tangent.noalias() +=
            measure * B.transpose() * C * B;

        const Mat3 S = stress.tensor();

        for (Index a = 0; a < total_dofs; ++a) {
            const bool a_beta = a >= ndof + principal_modes;
            const Index a_beta_index =
                a_beta ? a - ndof - principal_modes : Index(0);

            for (Index b = 0; b < total_dofs; ++b) {
                const bool b_beta = b >= ndof + principal_modes;
                const Index b_beta_index =
                    b_beta ? b - ndof - principal_modes : Index(0);

                Mat3 second_F = Mat3::Zero();

                if (a_beta && b_beta) {
                    second_F =
                        theta[static_cast<std::size_t>(a_beta_index)]
                        * theta[static_cast<std::size_t>(b_beta_index)]
                        * F_bar;
                } else if (a_beta) {
                    second_F =
                        theta[static_cast<std::size_t>(a_beta_index)]
                        * derivatives[static_cast<std::size_t>(b)];
                } else if (b_beta) {
                    second_F =
                        theta[static_cast<std::size_t>(b_beta_index)]
                        * derivatives[static_cast<std::size_t>(a)];
                }

                const Mat3 second_E =
                    derivatives[static_cast<std::size_t>(a)].transpose()
                    * derivatives[static_cast<std::size_t>(b)]
                    + F_bar.transpose() * second_F;

                result.tangent(a, b) +=
                    measure * stress_contraction(S, second_E);
            }
        }
    }

    logging::error(result.residual.allFinite(),
        "C3D8I: non-finite local residual in element ", elem_id);
    logging::error(!with_tangent || result.tangent.allFinite(),
        "C3D8I: non-finite local tangent in element ", elem_id);

    return result;
}

C3D8I::InternalSolution C3D8I::solve_internal(
    const Field& displacement,
    bool         with_tangent,
    bool         write_material_state,
    bool         update_trial_cache
) {
    Vector13 internal =
        update_trial_cache ? trial_internal_ : Vector13::Zero();

    NonlinearEvaluation evaluation;

    auto evaluate_candidate = [&](const Vector13& candidate) {
        return evaluate_nonlinear(displacement, candidate, true, false);
    };

    try {
        evaluation = evaluate_candidate(internal);
    } catch (...) {
        internal.setZero();
        evaluation = evaluate_candidate(internal);
    }

    Precision residual_norm =
        evaluation.residual.template tail<internal_dofs>().norm();
    const Precision force_energy_scale =
        characteristic_length_
        * evaluation.residual.template head<ndof>().norm();
    const Precision reference =
        std::max({residual_norm, force_energy_scale, Precision(1e-20)});

    bool converged = residual_norm <= local_tolerance * reference;

    for (Index iteration = 0;
         !converged && iteration < local_max_iterations;
         ++iteration) {
        const Vector13 r =
            evaluation.residual.template tail<internal_dofs>();
        const Matrix13 Kii =
            evaluation.tangent.template block<internal_dofs, internal_dofs>(
                ndof, ndof);

        Eigen::FullPivLU<Matrix13> lu(Kii);
        logging::error(lu.isInvertible(),
            "C3D8I: singular local incompatible-mode tangent in element ", elem_id);

        const Vector13 correction = -lu.solve(r);
        logging::error(correction.allFinite(),
            "C3D8I: non-finite incompatible-mode correction in element ", elem_id);

        bool accepted = false;
        Precision step = Precision(1);

        for (Index line_search = 0;
             line_search < local_max_line_search;
             ++line_search) {
            const Vector13 candidate = internal + step * correction;

            try {
                auto candidate_evaluation = evaluate_candidate(candidate);
                const Precision candidate_norm =
                    candidate_evaluation.residual
                        .template tail<internal_dofs>().norm();

                if (candidate_norm < residual_norm
                    || step <= Precision(1e-3)) {
                    internal = candidate;
                    evaluation = std::move(candidate_evaluation);
                    residual_norm = candidate_norm;
                    accepted = true;
                    break;
                }
            } catch (...) {
                // Reject an incompatible-mode iterate that creates an invalid
                // enhanced gradient and continue with a shorter local step.
            }

            step *= Precision(0.5);
        }

        logging::error(accepted,
            "C3D8I: local incompatible-mode Newton line search failed in element ",
            elem_id);

        converged = residual_norm <= local_tolerance * reference;
    }

    logging::error(converged,
        "C3D8I: incompatible-mode Newton did not converge in element ", elem_id,
        "\nresidual: ", residual_norm,
        "\nreference: ", reference);

    // Re-evaluate the converged local state with exactly the outputs requested by
    // the global assembly. This is the only evaluation allowed to write physical
    // trial material history.
    evaluation = evaluate_nonlinear(
        displacement, internal, with_tangent, write_material_state);

    if (update_trial_cache) {
        trial_internal_ = internal;
        trial_coords_ = node_coords_reference() + nodal_data<D>(displacement);
        trial_F_ = evaluation.enhanced_F;
        trial_valid_ = true;
    }

    return InternalSolution{internal, std::move(evaluation)};
}

void C3D8I::scatter_force(
    NodeData& nodal_forces,
    const Vector24& local_force
) {
    logging::error(nodal_forces.domain == FieldDomain::NODE
                && nodal_forces.components >= D,
        "C3D8I: invalid nonlinear force output");

    for (Index a = 0; a < N; ++a) {
        const Index node = static_cast<Index>(node_ids[a]);
        for (Dim d = 0; d < D; ++d) {
            nodal_forces(node, d) += local_force(D * a + d);
        }
    }
}

MapMatrix C3D8I::stiffness_tangent(
    Precision*   buffer,
    NodeData&    nodal_forces,
    const Field& displacement
) {
    const bool with_tangent = buffer != nullptr;
    auto solution = solve_internal(
        displacement, with_tangent, true, true);

    const Vector24 local_force =
        solution.evaluation.residual.template head<ndof>();
    scatter_force(nodal_forces, local_force);

    if (!with_tangent) {
        return MapMatrix(nullptr, 0, 0);
    }

    const Matrix37& K = solution.evaluation.tangent;
    const Matrix24 Kuu = K.template block<ndof, ndof>(0, 0);
    const StaticMatrix<ndof, internal_dofs> Kui =
        K.template block<ndof, internal_dofs>(0, ndof);
    const StaticMatrix<internal_dofs, ndof> Kiu =
        K.template block<internal_dofs, ndof>(ndof, 0);
    const Matrix13 Kii =
        K.template block<internal_dofs, internal_dofs>(ndof, ndof);

    Eigen::FullPivLU<Matrix13> lu(Kii);
    logging::error(lu.isInvertible(),
        "C3D8I: singular converged incompatible-mode tangent in element ", elem_id);

    const Matrix24 condensed = Kuu - Kui * lu.solve(Kiu);
    logging::error(condensed.allFinite(),
        "C3D8I: non-finite condensed nonlinear tangent in element ", elem_id);

    MapMatrix mapped{buffer, ndof, ndof};
    mapped = condensed;
    return mapped;
}

void C3D8I::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     displacement,
    const RowMatrix& rst,
    int              offset,
    bool             use_green_lagrange_nl
) {
    logging::error(strain != nullptr || stress != nullptr,
        "C3D8I: stress/strain recovery requires an output field");

    const auto& scheme = integration_scheme_stiffness();
    const RowMatrix ip_rst = stress_strain_ip_rst();
    const bool output_at_ip =
        rst.rows() == ip_rst.rows() && rst.leftCols(3).isApprox(ip_rst);
    const bool output_at_nodes = rst.rows() == static_cast<Eigen::Index>(N);

    logging::error(output_at_ip || output_at_nodes,
        "C3D8I: stress/strain output must use integration points or element nodes");

    RowMatrix ip_strain = RowMatrix::Zero(scheme.count(), 6);
    RowMatrix ip_stress = RowMatrix::Zero(scheme.count(), 6);

    if (use_green_lagrange_nl) {
        std::array<Mat3, N> enhanced_F;
        const auto reference_coords = node_coords_reference();
        const auto current_coords =
            reference_coords + nodal_data<D>(displacement);

        if (nonlinear_state_initialized_
            && current_coords.isApprox(
                committed_coords_, Precision(1e-10))) {
            enhanced_F = committed_F_;
        } else {
            enhanced_F = solve_internal(
                displacement, false, false, false).evaluation.enhanced_F;
        }

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto point = scheme.get_point(ip);
            const Mat3& F = enhanced_F[static_cast<std::size_t>(ip)];
            const VolumeStrainGreenLagrange E =
                VolumeStrainGreenLagrange::from_deformation_gradient(F);

            const Index state_row = mp_index(ip);
            const Precision* old_state =
                &(*_model_data->material_state_old)(state_row, 0);

            VolumeStressPK2 S;
            evaluate_material(
                point.r, point.s, point.t,
                E, old_state, nullptr, S, nullptr);
            const VolumeStressCauchy sigma = S.to_cauchy(F);

            ip_strain.row(ip) = E.voigt().transpose();
            ip_stress.row(ip) = sigma.voigt().transpose();
        }
    } else {
        const auto reference_coords = node_coords_reference();
        const Vector13 internal = linear_internal_parameters(displacement);

        Vector24 u = Vector24::Zero();
        const auto local = nodal_data<D>(displacement);
        for (Index a = 0; a < N; ++a) {
            for (Dim d = 0; d < D; ++d) {
                u(D * a + d) = local(a, d);
            }
        }

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto point = scheme.get_point(ip);

            Precision det0 = Precision(0);
            const auto dN_dX = shape_derivatives_reference(
                reference_coords, point.r, point.s, point.t, det0);
            const auto Bu = strain_displacement(dN_dX);

            Precision incompatible_det = Precision(0);
            const auto Bi = linear_internal_B(
                point.r, point.s, point.t,
                reference_coords, incompatible_det);

            const Vec6 values = Bu * u + Bi * internal;
            const VolumeStrainLinearized eps(values);

            const Index state_row = mp_index(ip);
            const Precision* old_state =
                &(*_model_data->material_state_old)(state_row, 0);

            VolumeStressCauchy sigma;
            Mat6 C;
            evaluate_material(
                point.r, point.s, point.t,
                eps, old_state, nullptr, sigma, C);

            ip_strain.row(ip) = eps.voigt().transpose();
            ip_stress.row(ip) = sigma.voigt().transpose();
        }
    }

    if (output_at_ip) {
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const Index row = static_cast<Index>(offset) + ip;
            for (Dim component = 0; component < 6; ++component) {
                if (strain) (*strain)(row, component) = ip_strain(ip, component);
                if (stress) (*stress)(row, component) = ip_stress(ip, component);
            }
        }
        return;
    }

    const RowMatrix& E = extrapolation_matrix();
    const RowMatrix nodal_strain = E * ip_strain;
    const RowMatrix nodal_stress = E * ip_stress;

    for (Index node = 0; node < N; ++node) {
        const Index row = static_cast<Index>(offset) + node;
        for (Dim component = 0; component < 6; ++component) {
            if (strain) (*strain)(row, component) = nodal_strain(node, component);
            if (stress) (*stress)(row, component) = nodal_stress(node, component);
        }
    }
}

} // namespace fem::model
