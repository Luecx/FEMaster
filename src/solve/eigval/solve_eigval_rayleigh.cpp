#include "solve_eigval_rayleigh.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace fem::solver {

Precision estimate_lambda_scale_from_geometry(
    const SparseMatrix&              A,
    const SparseMatrix&              B,
    const model::Field&              positions,
    const IndexMatrix&               active_dof_idx_mat,
    const constraint::ConstraintMap& map)
{
    if (A.rows() == 0 || A.rows() != A.cols() || B.rows() != A.rows() || B.cols() != A.cols())
        return Precision(0);

    if (positions.rows != active_dof_idx_mat.rows() || positions.components < 3)
        return Precision(0);

    if (map.n_master() != static_cast<Index>(A.rows()))
        return Precision(0);

    const Precision eps    = std::numeric_limits<Precision>::epsilon();
    const Precision A_norm = A.norm();
    const Precision B_norm = B.norm();
    if (!std::isfinite(A_norm) || !std::isfinite(B_norm) || A_norm <= Precision(0) || B_norm <= Precision(0))
        return Precision(0);

    // Determine the geometric center and spans. Physical coordinates relative
    // to this center are used for rigid rotations, while normalized coordinates
    // are used for the quadratic Ritz fields.
    std::array<Precision, 3> min_pos{
        std::numeric_limits<Precision>::infinity(),
        std::numeric_limits<Precision>::infinity(),
        std::numeric_limits<Precision>::infinity()
    };
    std::array<Precision, 3> max_pos{
        -std::numeric_limits<Precision>::infinity(),
        -std::numeric_limits<Precision>::infinity(),
        -std::numeric_limits<Precision>::infinity()
    };

    for (Index node = 0; node < positions.rows; ++node) {
        for (int axis = 0; axis < 3; ++axis) {
            const Precision value = positions(node, axis);
            if (!std::isfinite(value)) continue;
            min_pos[axis] = std::min(min_pos[axis], value);
            max_pos[axis] = std::max(max_pos[axis], value);
        }
    }

    std::array<Precision, 3> center{};
    std::array<Precision, 3> half_span{};
    for (int axis = 0; axis < 3; ++axis) {
        if (!std::isfinite(min_pos[axis]) || !std::isfinite(max_pos[axis]))
            return Precision(0);

        center[axis]    = Precision(0.5) * (min_pos[axis] + max_pos[axis]);
        half_span[axis] = Precision(0.5) * (max_pos[axis] - min_pos[axis]);
    }

    const Precision length_scale = std::max({
        Precision(2) * half_span[0],
        Precision(2) * half_span[1],
        Precision(2) * half_span[2]
    });
    if (!std::isfinite(length_scale) || length_scale <= eps)
        return Precision(0);

    // Per node: normalized quadratic monomials and their gradients with respect
    // to the physical coordinates. The gradients are needed to add compatible
    // nodal rotations theta = 1/2 curl(u) for beam/shell rotational DOFs.
    using GradientSet = std::array<std::array<Precision, 3>, 6>;
    std::vector<std::array<Precision, 3>> relative_positions(static_cast<std::size_t>(positions.rows));
    std::vector<std::array<Precision, 6>> monomials(static_cast<std::size_t>(positions.rows));
    std::vector<GradientSet>              gradients(static_cast<std::size_t>(positions.rows));

    for (Index node = 0; node < positions.rows; ++node) {
        std::array<Precision, 3> r{};
        std::array<Precision, 3> p{};
        std::array<Precision, 3> inv_span{};

        for (int axis = 0; axis < 3; ++axis) {
            r[axis] = positions(node, axis) - center[axis];
            if (half_span[axis] > eps) {
                p[axis]        = r[axis] / half_span[axis];
                inv_span[axis] = Precision(1) / half_span[axis];
            } else {
                p[axis]        = Precision(0);
                inv_span[axis] = Precision(0);
            }
        }

        relative_positions[static_cast<std::size_t>(node)] = r;

        const Precision x = p[0];
        const Precision y = p[1];
        const Precision z = p[2];

        monomials[static_cast<std::size_t>(node)] = {
            x * x,
            y * y,
            z * z,
            x * y,
            x * z,
            y * z
        };

        gradients[static_cast<std::size_t>(node)] = {{
            {Precision(2) * x * inv_span[0], Precision(0),                         Precision(0)},
            {Precision(0),                         Precision(2) * y * inv_span[1], Precision(0)},
            {Precision(0),                         Precision(0),                         Precision(2) * z * inv_span[2]},
            {y * inv_span[0],                     x * inv_span[1],                     Precision(0)},
            {z * inv_span[0],                     Precision(0),                         x * inv_span[2]},
            {Precision(0),                         z * inv_span[1],                     y * inv_span[2]}
        }};
    }

    // Reverse active DOF numbering so every reduced master coordinate can be
    // associated with its nodal position and physical DOF direction.
    const int full_size = active_dof_idx_mat.size() > 0
        ? active_dof_idx_mat.maxCoeff() + 1
        : 0;
    if (full_size <= 0)
        return Precision(0);

    std::vector<int> node_of_gid(static_cast<std::size_t>(full_size), -1);
    std::vector<int> dof_of_gid (static_cast<std::size_t>(full_size), -1);

    for (int node = 0; node < active_dof_idx_mat.rows(); ++node) {
        for (int dof = 0; dof < active_dof_idx_mat.cols(); ++dof) {
            const int gid = active_dof_idx_mat(node, dof);
            if (gid < 0 || gid >= full_size) continue;
            node_of_gid[static_cast<std::size_t>(gid)] = node;
            dof_of_gid [static_cast<std::size_t>(gid)] = dof;
        }
    }

    // Build the admissible rigid-body nullspace in reduced coordinates. A
    // candidate is kept only if it is numerically in the nullspace of A, which
    // prevents supported rigid motions from removing flexible Ritz content.
    std::vector<DynamicVector> rigid_basis;
    std::vector<DynamicVector> rigid_basis_B;
    rigid_basis.reserve(6);
    rigid_basis_B.reserve(6);

    constexpr Precision rigid_residual_tol = static_cast<Precision>(1e-8);

    for (int mode = 0; mode < 6; ++mode) {
        DynamicVector q = DynamicVector::Zero(A.rows());

        for (Index master = 0; master < map.n_master(); ++master) {
            const Index gid = map.masters[static_cast<std::size_t>(master)];
            if (gid < 0 || gid >= full_size) continue;

            const int node = node_of_gid[static_cast<std::size_t>(gid)];
            const int dof  = dof_of_gid [static_cast<std::size_t>(gid)];
            if (node < 0 || dof < 0 || dof >= 6) continue;

            const auto& r = relative_positions[static_cast<std::size_t>(node)];
            const Precision x = r[0];
            const Precision y = r[1];
            const Precision z = r[2];

            Precision value = Precision(0);
            switch (mode) {
                case 0: value = dof == 0 ? Precision(1) : Precision(0); break; // Tx
                case 1: value = dof == 1 ? Precision(1) : Precision(0); break; // Ty
                case 2: value = dof == 2 ? Precision(1) : Precision(0); break; // Tz
                case 3: // Rx: omega x r = (0, -z, y), plus nodal rx
                    if      (dof == 1) value = -z;
                    else if (dof == 2) value =  y;
                    else if (dof == 3) value = Precision(1);
                    break;
                case 4: // Ry: omega x r = (z, 0, -x), plus nodal ry
                    if      (dof == 0) value =  z;
                    else if (dof == 2) value = -x;
                    else if (dof == 4) value = Precision(1);
                    break;
                case 5: // Rz: omega x r = (-y, x, 0), plus nodal rz
                    if      (dof == 0) value = -y;
                    else if (dof == 1) value =  x;
                    else if (dof == 5) value = Precision(1);
                    break;
            }

            q[master] = value;
        }

        const Precision q_norm = q.norm();
        if (!std::isfinite(q_norm) || q_norm <= eps)
            continue;
        q /= q_norm;

        const DynamicVector Aq = A * q;
        const Precision relative_residual = Aq.norm() / A_norm;
        if (!std::isfinite(relative_residual) || relative_residual > rigid_residual_tol)
            continue;

        DynamicVector Bq = B * q;

        // Modified Gram-Schmidt in the B inner product. Two passes keep the
        // small rigid basis stable for mixed translational/rotational units.
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t j = 0; j < rigid_basis.size(); ++j) {
                const Precision projection = rigid_basis[j].dot(Bq);
                q  -= projection * rigid_basis[j];
                Bq -= projection * rigid_basis_B[j];
            }
        }

        const Precision b_norm_sq = q.dot(Bq);
        if (!std::isfinite(b_norm_sq) || b_norm_sq <= Precision(100) * eps * B_norm)
            continue;

        const Precision inv_b_norm = Precision(1) / std::sqrt(b_norm_sq);
        q  *= inv_b_norm;
        Bq *= inv_b_norm;

        rigid_basis.push_back(std::move(q));
        rigid_basis_B.push_back(std::move(Bq));
    }

    // Build a compact quadratic Ritz basis. The six scalar monomials x², y²,
    // z², xy, xz and yz are applied in each translational direction. Beam/shell
    // rotations follow from theta = 1/2 curl(u). Every vector is first deflated
    // against admissible rigid modes and then B-orthonormalized against the
    // already accepted Ritz vectors.
    std::vector<DynamicVector> ritz_basis;
    std::vector<DynamicVector> ritz_basis_B;
    ritz_basis.reserve(18);
    ritz_basis_B.reserve(18);

    for (int polynomial = 0; polynomial < 6; ++polynomial) {
        for (int direction = 0; direction < 3; ++direction) {
            DynamicVector q = DynamicVector::Zero(A.rows());

            for (Index master = 0; master < map.n_master(); ++master) {
                const Index gid = map.masters[static_cast<std::size_t>(master)];
                if (gid < 0 || gid >= full_size) continue;

                const int node = node_of_gid[static_cast<std::size_t>(gid)];
                const int dof  = dof_of_gid [static_cast<std::size_t>(gid)];
                if (node < 0 || dof < 0 || dof >= 6) continue;

                const Precision f = monomials[static_cast<std::size_t>(node)][polynomial];
                const auto& grad  = gradients[static_cast<std::size_t>(node)][polynomial];

                Precision value = Precision(0);

                // Translation u = L f e_direction.
                if (dof == direction) {
                    value = length_scale * f;
                }

                // Compatible nodal rotation theta = 1/2 curl(u).
                if (direction == 0) {
                    if      (dof == 4) value =  Precision(0.5) * length_scale * grad[2];
                    else if (dof == 5) value = -Precision(0.5) * length_scale * grad[1];
                } else if (direction == 1) {
                    if      (dof == 3) value = -Precision(0.5) * length_scale * grad[2];
                    else if (dof == 5) value =  Precision(0.5) * length_scale * grad[0];
                } else {
                    if      (dof == 3) value =  Precision(0.5) * length_scale * grad[1];
                    else if (dof == 4) value = -Precision(0.5) * length_scale * grad[0];
                }

                q[master] = value;
            }

            if (q.squaredNorm() <= eps)
                continue;

            DynamicVector Bq = B * q;

            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t j = 0; j < rigid_basis.size(); ++j) {
                    const Precision projection = rigid_basis[j].dot(Bq);
                    q  -= projection * rigid_basis[j];
                    Bq -= projection * rigid_basis_B[j];
                }

                for (std::size_t j = 0; j < ritz_basis.size(); ++j) {
                    const Precision projection = ritz_basis[j].dot(Bq);
                    q  -= projection * ritz_basis[j];
                    Bq -= projection * ritz_basis_B[j];
                }
            }

            const Precision b_norm_sq = q.dot(Bq);
            if (!std::isfinite(b_norm_sq) || b_norm_sq <= Precision(100) * eps * B_norm)
                continue;

            const Precision inv_b_norm = Precision(1) / std::sqrt(b_norm_sq);
            q  *= inv_b_norm;
            Bq *= inv_b_norm;

            ritz_basis.push_back(std::move(q));
            ritz_basis_B.push_back(std::move(Bq));
        }
    }

    if (ritz_basis.empty())
        return Precision(0);

    const Eigen::Index m = static_cast<Eigen::Index>(ritz_basis.size());
    DynamicMatrix Q(A.rows(), m);
    DynamicMatrix BQ(B.rows(), m);

    for (Eigen::Index j = 0; j < m; ++j) {
        Q.col(j)  = ritz_basis[static_cast<std::size_t>(j)];
        BQ.col(j) = ritz_basis_B[static_cast<std::size_t>(j)];
    }

    // Rayleigh-Ritz in the quadratic subspace. Because Q is already
    // B-orthonormal, Mr is close to identity, but forming it explicitly keeps
    // the projected problem correct despite finite-precision orthogonalization.
    const DynamicMatrix AQ = A * Q;
    DynamicMatrix Kr = Q.transpose() * AQ;
    DynamicMatrix Mr = Q.transpose() * BQ;

    Kr = Precision(0.5) * (Kr + Kr.transpose()).eval();
    Mr = Precision(0.5) * (Mr + Mr.transpose()).eval();

    Eigen::GeneralizedSelfAdjointEigenSolver<DynamicMatrix> eig(Kr, Mr, Eigen::EigenvaluesOnly);
    if (eig.info() != Eigen::Success)
        return Precision(0);

    for (Eigen::Index i = 0; i < eig.eigenvalues().size(); ++i) {
        const Precision lambda = eig.eigenvalues()[i];
        if (std::isfinite(lambda) && lambda > Precision(0))
            return lambda;
    }

    return Precision(0);
}

} // namespace fem::solver
