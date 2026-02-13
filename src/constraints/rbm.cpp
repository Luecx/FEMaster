/**
 * @file rbm.cpp
 * @brief Implements the RBM (rigid-body mode) constraint builder.
 */

#include "rbm.h"

#include "../core/logging.h"
#include "../core/types_eig.h"
#include "../model/model_data.h"
#include "../model/element/element_structural.h"

#include <algorithm>

namespace fem {
namespace constraint {

namespace {

// Collect unique node ids from all elements in the configured element set.
static std::vector<ID> collect_nodes(const Rbm& rbm, model::ModelData& model_data) {
    std::vector<ID> ids;
    ids.reserve(256);

    if (rbm.element_region) {
        for (ID elem_id : *rbm.element_region) {
            if (elem_id < 0 || static_cast<std::size_t>(elem_id) >= model_data.elements.size()) {
                continue;
            }

            auto& eptr = model_data.elements[static_cast<std::size_t>(elem_id)];
            if (!eptr) {
                continue;
            }

            for (ID node_id : *eptr) {
                ids.push_back(node_id);
            }
        }
    }

    // unique and sorted
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}

// Compute center of gravity of all structural elements in the configured set.
static Vec3 compute_center_of_gravity(const Rbm& rbm, model::ModelData& model_data) {
    Precision m = Precision(0);
    Vec3 c_num  = Vec3::Zero();

    if (rbm.element_region) {
        for (ID elem_id : *rbm.element_region) {
            if (elem_id < 0 || static_cast<std::size_t>(elem_id) >= model_data.elements.size()) {
                continue;
            }

            auto& eptr = model_data.elements[static_cast<std::size_t>(elem_id)];
            if (!eptr) {
                continue;
            }

            auto* se = eptr->as<model::StructuralElement>();
            if (!se) {
                continue;
            }

            m += se->integrate_scalar_field(true, [](const Vec3&) { return Precision(1); });
            c_num += se->integrate_vector_field(true, [](const Vec3& x) { return x; });
        }
    }

    logging::error(m > Precision(0),
                   "RBM: selected element set has zero integrated mass; center of gravity cannot be computed");
    return c_num / m;
}

static bool has_any_translation_dof(const SystemDofIds& dofs, ID node_id) {
    return dofs(node_id, 0) >= 0 || dofs(node_id, 1) >= 0 || dofs(node_id, 2) >= 0;
}

static bool add_if_available(Equation& equation,
                             const SystemDofIds& dofs,
                             ID node_id,
                             Dim dof,
                             Precision coeff) {
    if (dofs(node_id, dof) < 0) {
        return false;
    }
    equation.entries.push_back(EquationEntry{node_id, dof, coeff});
    return true;
}

} // namespace

Equations Rbm::get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) const {
    Equations eqs;

    logging::error(model_data.positions != nullptr, "RBM: positions field not set in model data");
    const auto& X = *model_data.positions;

    // Collect candidate nodes from all elements in the region.
    std::vector<ID> ids = collect_nodes(*this, model_data);
    if (ids.empty()) return eqs;

    // Keep only nodes that contribute at least one translational DOF.
    std::vector<ID> use;
    use.reserve(ids.size());
    for (ID id : ids) {
        if (has_any_translation_dof(system_nodal_dofs, id)) {
            use.push_back(id);
        }
    }
    if (use.empty()) return eqs;

    // Use center of gravity of the selected element set as rotation center.
    const Vec3 x0 = compute_center_of_gravity(*this, model_data);

    // Optional normalization keeps coefficients comparable for large sets.
    const Precision w = Precision(1) / Precision(std::max<std::size_t>(1, use.size()));

    // 3 translation equations: sum u = 0 (scaled)
    Equation ex, ey, ez;
    ex.entries.reserve(use.size());
    ey.entries.reserve(use.size());
    ez.entries.reserve(use.size());

    for (ID i : use) {
        add_if_available(ex, system_nodal_dofs, i, 0, w);
        add_if_available(ey, system_nodal_dofs, i, 1, w);
        add_if_available(ez, system_nodal_dofs, i, 2, w);
    }

    if (!ex.entries.empty()) eqs.emplace_back(std::move(ex));
    if (!ey.entries.empty()) eqs.emplace_back(std::move(ey));
    if (!ez.entries.empty()) eqs.emplace_back(std::move(ez));

    // 3 rotation equations: sum (r x u) = 0, with r measured about x0.
    Equation erx, ery, erz;
    erx.entries.reserve(2 * use.size());
    ery.entries.reserve(2 * use.size());
    erz.entries.reserve(2 * use.size());

    for (ID i : use) {
        const Vec3 xi = X.row_vec3(static_cast<Index>(i));
        const Vec3 r  = xi - x0;
        const Precision rx = r(0), ry = r(1), rz = r(2);

        // (r x u)_x =  ry * u_z - rz * u_y
        add_if_available(erx, system_nodal_dofs, i, 2,  w * ry);
        add_if_available(erx, system_nodal_dofs, i, 1, -w * rz);

        // (r x u)_y =  rz * u_x - rx * u_z
        add_if_available(ery, system_nodal_dofs, i, 0,  w * rz);
        add_if_available(ery, system_nodal_dofs, i, 2, -w * rx);

        // (r x u)_z =  rx * u_y - ry * u_x
        add_if_available(erz, system_nodal_dofs, i, 1,  w * rx);
        add_if_available(erz, system_nodal_dofs, i, 0, -w * ry);
    }

    if (!erx.entries.empty()) eqs.emplace_back(std::move(erx));
    if (!ery.entries.empty()) eqs.emplace_back(std::move(ery));
    if (!erz.entries.empty()) eqs.emplace_back(std::move(erz));

    return eqs;
}

} // namespace constraint
} // namespace fem
