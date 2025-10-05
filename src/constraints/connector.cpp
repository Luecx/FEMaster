/******************************************************************************
 * @file connector.cpp
 * @brief Implements connector constraints that tie node DOFs.
 *
 * Connector objects convert bit-mask specifications into linear constraint
 * equations that act on global degrees of freedom.
 *
 * @see src/constraints/connector.h
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "connector.h"

#include "../core/logging.h"
#include "../model/model_data.h"

#include <cmath>
#include <utility>

namespace fem {
namespace constraint {
namespace {
constexpr Precision kEps = static_cast<Precision>(1e-8); ///< Numerical tolerance for projections.
}

/******************************************************************************
 * @copydoc Connector::Connector
 ******************************************************************************/
Connector::Connector(ID node_1_id,
                     ID node_2_id,
                     fem::cos::CoordinateSystem::Ptr coordinate_system,
                     ConnectorType connector_type)
    : coordinate_system_(std::move(coordinate_system))
    , connector_type_(connector_type)
    , constrained_dofs_({
          static_cast<bool>(connector_type_ & 0b100000),
          static_cast<bool>(connector_type_ & 0b010000),
          static_cast<bool>(connector_type_ & 0b001000),
          static_cast<bool>(connector_type_ & 0b000100),
          static_cast<bool>(connector_type_ & 0b000010),
          static_cast<bool>(connector_type_ & 0b000001)})
    , node_1_id_(node_1_id)
    , node_2_id_(node_2_id) {}

/******************************************************************************
 * @copydoc Connector::get_equations
 ******************************************************************************/
Equations Connector::get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) {
    (void)model_data;
    Equations equations{};

    const char axis_names[3] = {'X', 'Y', 'Z'};

    for (int i = 0; i < 6; ++i) {
        if (!constrained_dofs_(i)) {
            continue;
        }

        const Vec3 local_direction = (i < 3) ? Vec3::Unit(i) : Vec3::Unit(i - 3);
        const Vec3 global_direction = coordinate_system_->to_global(local_direction);

        std::vector<EquationEntry> entries;
        const int block = (i < 3) ? 0 : 1; // 0: translations, 1: rotations

        for (int j = 0; j < 3; ++j) {
            const Precision coeff = global_direction(j);
            if (std::abs(coeff) <= kEps) {
                continue;
            }

            const Dim dof = static_cast<Dim>(3 * block + j);
            const auto dof_id_n1 = system_nodal_dofs(node_1_id_, dof);
            const auto dof_id_n2 = system_nodal_dofs(node_2_id_, dof);
            if (dof_id_n1 < 0 || dof_id_n2 < 0) {
                continue;
            }

            entries.emplace_back(EquationEntry{node_1_id_, dof, coeff});
            entries.emplace_back(EquationEntry{node_2_id_, dof, -coeff});
        }

        if (!entries.empty()) {
            equations.emplace_back(Equation{std::move(entries)});
        }

        logging::warning(entries.empty(),
                         "Connector (", node_1_id_, " <-> ", node_2_id_,
                         ") generated no equation for local ",
                         (i < 3 ? "T" : "R"), axis_names[(i < 3) ? i : i - 3],
                         " (", global_direction(0), ", ", global_direction(1), ", ", global_direction(2), "). ",
                         "Coordinate system '", coordinate_system_->name, "'.");
    }

    return equations;
}

/******************************************************************************
 * @copydoc Connector::dofs
 ******************************************************************************/
Dofs Connector::dofs() const {
    Dofs impacted;
    impacted.setZero();

    for (int i = 0; i < 6; ++i) {
        if (!constrained_dofs_(i)) {
            continue;
        }

        const Vec3 local_dir = (i < 3) ? Vec3::Unit(i) : Vec3::Unit(i - 3);
        const Vec3 global_dir = coordinate_system_->to_global(local_dir);
        const int block = (i < 3) ? 0 : 1;

        for (int j = 0; j < 3; ++j) {
            if (std::abs(global_dir(j)) > kEps) {
                impacted(static_cast<int>(3 * block + j)) = true;
            }
        }
    }

    return impacted;
}

} // namespace constraint
} // namespace fem
