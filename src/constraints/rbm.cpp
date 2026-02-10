/**
 * @file rbm.cpp
 * @brief Implements the RBM (rigid-body mode) constraint builder.
 */

#include "rbm.h"

#include "../core/logging.h"
#include "../model/model_data.h"

#include <algorithm>
#include <limits>

namespace fem {
namespace constraint {

namespace {

// Collect unique node ids from the configured region(s)
static std::vector<ID> collect_nodes(const Rbm& rbm, model::ModelData& /*model_data*/) {
    std::vector<ID> ids;
    ids.reserve(256);

    if (rbm.node_region) {
        for (ID id : *rbm.node_region) ids.push_back(id);
    }

    // unique and sorted
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}

// Farthest-point sampling on node coordinates
static std::vector<ID> farthest_point_sample(const std::vector<ID>& ids,
                                             const model::Field& X,
                                             std::size_t k) {
    std::vector<ID> out;
    out.reserve(std::min(k, ids.size()));
    if (ids.empty() || k == 0) return out;

    // Compute centroid
    Vec3 c = Vec3::Zero();
    for (ID id : ids) c += X.row_vec3(static_cast<Index>(id));
    c /= Precision(std::max<std::size_t>(1, ids.size()));

    // Pick first: farthest from centroid
    auto argmax_from_c = [&]() -> std::size_t {
        Precision best = -std::numeric_limits<Precision>::infinity();
        std::size_t idx = 0;
        for (std::size_t i = 0; i < ids.size(); ++i) {
            const Vec3 p = X.row_vec3(static_cast<Index>(ids[i]));
            const Precision d2 = (p - c).squaredNorm();
            if (d2 > best) { best = d2; idx = i; }
        }
        return idx;
    }();

    out.push_back(ids[argmax_from_c]);

    // Running min distance to selected set
    std::vector<Precision> dmin(ids.size(), std::numeric_limits<Precision>::infinity());
    auto update_dmin = [&](ID sel) {
        const Vec3 ps = X.row_vec3(static_cast<Index>(sel));
        for (std::size_t i = 0; i < ids.size(); ++i) {
            const Vec3 pi = X.row_vec3(static_cast<Index>(ids[i]));
            const Precision d2 = (pi - ps).squaredNorm();
            dmin[i] = std::min(dmin[i], d2);
        }
    };
    update_dmin(out.back());

    while (out.size() < std::min(k, ids.size())) {
        // pick idx with largest dmin
        Precision best = -std::numeric_limits<Precision>::infinity();
        std::size_t j = 0;
        for (std::size_t i = 0; i < ids.size(); ++i) {
            if (dmin[i] > best) { best = dmin[i]; j = i; }
        }
        const ID pick = ids[j];
        // Avoid duplicates
        if (std::find(out.begin(), out.end(), pick) == out.end()) {
            out.push_back(pick);
            update_dmin(pick);
        } else {
            // advance to next distinct candidate
            dmin[j] = -std::numeric_limits<Precision>::infinity();
        }
        if (out.size() == ids.size()) break;
    }
    return out;
}

} // namespace

Equations Rbm::get_equations(SystemDofIds& /*system_nodal_dofs*/, model::ModelData& model_data) const {
    Equations eqs;

    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    const auto& X = *model_data.positions;

    // Collect candidate nodes
    std::vector<ID> ids = collect_nodes(*this, model_data);
    if (ids.empty()) return eqs;

    // Decide sampling budget (at least 3 points to define rotations)
    const std::size_t K = static_cast<std::size_t>(
        std::max<Index>(3, (max_points > 0 ? max_points : 16)));

    std::vector<ID> use = ids;
    if (ids.size() > K) {
        use = farthest_point_sample(ids, X, K);
    }

    // Centroid x0 of selected nodes
    Vec3 x0 = Vec3::Zero();
    for (ID i : use) x0 += X.row_vec3(static_cast<Index>(i));
    x0 /= Precision(std::max<std::size_t>(1, use.size()));

    // Helper to add entry
    auto add_entry = [](Equation& e, ID nid, Dim dof, Precision c) {
        e.entries.push_back(EquationEntry{nid, dof, c});
    };

    // 3 translation equations: sum u = 0
    Equation ex, ey, ez;
    for (ID i : use) {
        add_entry(ex, i, 0, Precision(1));
        add_entry(ey, i, 1, Precision(1));
        add_entry(ez, i, 2, Precision(1));
    }
    eqs.emplace_back(std::move(ex));
    eqs.emplace_back(std::move(ey));
    eqs.emplace_back(std::move(ez));

    // 3 rotation equations: sum r x u = 0
    Equation erx, ery, erz; // components of sum (r x u)
    for (ID i : use) {
        const Vec3 xi = X.row_vec3(static_cast<Index>(i));
        const Vec3 r  = xi - x0;
        const Precision rx = r(0), ry = r(1), rz = r(2);

        // (r x u)_x =  ry * u_z - rz * u_y
        add_entry(erx, i, 2,  ry);
        add_entry(erx, i, 1, -rz);
        // (r x u)_y =  rz * u_x - rx * u_z
        add_entry(ery, i, 0,  rz);
        add_entry(ery, i, 2, -rx);
        // (r x u)_z =  rx * u_y - ry * u_x
        add_entry(erz, i, 1,  rx);
        add_entry(erz, i, 0, -ry);
    }
    eqs.emplace_back(std::move(erx));
    eqs.emplace_back(std::move(ery));
    eqs.emplace_back(std::move(erz));

    return eqs;
}

} // namespace constraint
} // namespace fem
