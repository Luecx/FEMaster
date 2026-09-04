#include "solve_eigval_rayleigh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
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

    // Normalize coordinates to approximately [-1, 1] per active geometric axis.
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

    std::vector<std::array<Precision, 6>> monomials(static_cast<std::size_t>(positions.rows));
    for (Index node = 0; node < positions.rows; ++node) {
        std::array<Precision, 3> p{};
        for (int axis = 0; axis < 3; ++axis) {
            const Precision scale = half_span[axis];
            p[axis] = scale > std::numeric_limits<Precision>::epsilon()
                ? (positions(node, axis) - center[axis]) / scale
                : Precision(0);
        }

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

    // Reverse active DOF numbering so each reduced master coordinate can be
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

    Precision lambda_min = std::numeric_limits<Precision>::infinity();

    // Six quadratic scalar fields, each applied independently in ux, uy and uz.
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

            if (q.squaredNorm() <= std::numeric_limits<Precision>::epsilon())
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
