#include "solve_eigval_rayleigh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace fem::solver {

Precision estimate_lambda_scale_from_geometry(
    const SparseMatrix&               A,
    const SparseMatrix&               B,
    const model::Field&               positions,
    const IndexMatrix&                active_dof_idx_mat,
    const constraint::ConstraintMap&  map)
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
    // are used for the quadratic trial fields.
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

    std::vector<std::array<Precision, 3>> relative_positions(static_cast<std::size_t>(positions.rows));
    std::vector<std::array<Precision, 6>> monomials(static_cast<std::size_t>(positions.rows));

    for (Index node = 0; node < positions.rows; ++node) {
        std::array<Precision, 3> r{};
        std::array<Precision, 3> p{};

        for (int axis = 0; axis < 3; ++axis) {
            r[axis] = positions(node, axis) - center[axis];
            p[axis] = half_span[axis] > eps ? r[axis] / half_span[axis] : Precision(0);
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
    // prevents supported or otherwise constrained rigid motions from removing
    // flexible content from the Rayleigh trial fields.
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
        // small six-vector basis sufficiently orthogonal even for mixed units.
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

    Precision lambda_min = std::numeric_limits<Precision>::infinity();

    // Six quadratic scalar fields, each applied independently in ux, uy and uz.
    // Before evaluating the Rayleigh quotient, remove every admissible rigid-body
    // component in the B inner product so translations and rotations cannot pull
    // the estimate down toward the numerical zero modes.
    for (int polynomial = 0; polynomial < 6; ++polynomial) {
        for (int direction = 0; direction < 3; ++direction) {
            DynamicVector q = DynamicVector::Zero(A.rows());

            for (Index master = 0; master < map.n_master(); ++master) {
                const Index gid = map.masters[static_cast<std::size_t>(master)];
                if (gid < 0 || gid >= full_size) continue;

                const int node = node_of_gid[static_cast<std::size_t>(gid)];
                const int dof  = dof_of_gid [static_cast<std::size_t>(gid)];
                if (node < 0 || dof != direction) continue;

                q[master] = monomials[static_cast<std::size_t>(node)][polynomial];
            }

            if (q.squaredNorm() <= eps)
                continue;

            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t j = 0; j < rigid_basis.size(); ++j) {
                    const Precision projection = rigid_basis_B[j].dot(q);
                    q -= projection * rigid_basis[j];
                }
            }

            if (q.squaredNorm() <= eps)
                continue;

            const DynamicVector Aq = A * q;
            const DynamicVector Bq = B * q;
            const Precision num = q.dot(Aq);
            const Precision den = q.dot(Bq);

            if (!std::isfinite(num) || !std::isfinite(den) || num <= Precision(0) || den <= Precision(0))
                continue;

            const Precision lambda = num / den;
            if (std::isfinite(lambda) && lambda > Precision(0))
                lambda_min = std::min(lambda_min, lambda);
        }
    }

    return std::isfinite(lambda_min) ? lambda_min : Precision(0);
}

} // namespace fem::solver
