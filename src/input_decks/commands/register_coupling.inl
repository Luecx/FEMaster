// register_coupling.inl — registers *COUPLING
#pragma once
/**
 * @file register_coupling.inl
 * @brief Register *COUPLING command.
 *
 * Couplings tie DOFs of a SLAVE (node set or surface) to a MASTER (usually one node).
 * Types:
 *  • KINEMATIC  — slaves follow the master kinematically (rigid “spider” behavior).
 *  • STRUCTURAL — loads/moments are distributed/collected work-equivalently.
 *
 * Notes:
 *  • Exactly one of SLAVE (node set) or SFSET (surface) must be given.
 *  • MASTER should usually contain exactly one node; larger masters can over-constrain.
 *  • Only couple DOFs that exist on involved nodes (e.g., truss nodes have no rotations).
 *  • Near-collinear axis/geometry definitions can yield unrealistic behavior and must be avoided.
 */

#include <array>
#include <memory>
#include <string>

#include "../../constraints/coupling.h"
#include "../../core/logging.h"
#include "../../core/types_eig.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_coupling(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("COUPLING", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Create kinematic or structural couplings between a MASTER and a SLAVE (node set) "
            "or SFSET (surface). Provide exactly one of SLAVE or SFSET. DOF entries > 0 enable "
            "coupling for [Ux, Uy, Uz, Rx, Ry, Rz]. KINEMATIC makes slaves follow the master; "
            "STRUCTURAL distributes/collects loads in a work-equivalent way. Avoid near-collinear "
            "axis/geometry definitions as they can cause unrealistic results."
        );

        // Persistent per-command state
        auto master  = std::make_shared<std::string>();
        auto slave   = std::make_shared<std::string>();
        auto surface = std::make_shared<std::string>();
        auto type    = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MASTER").required().doc("Master node set (commonly exactly one node).")
                .key("TYPE").required().allowed({"KINEMATIC", "STRUCTURAL"})
                .key("SLAVE").optional().doc("Slave node set.")
                .key("SFSET").optional().doc("Slave surface set.")
        );

        command.on_enter([&model, master, type, slave, surface](const fem::dsl::Keys& keys) {
            *master  = keys.raw("MASTER");
            *type    = keys.raw("TYPE");
            *slave   = keys.has("SLAVE") ? keys.raw("SLAVE") : std::string{};
            *surface = keys.has("SFSET") ? keys.raw("SFSET") : std::string{};

            // Exactly one target (SLAVE xor SFSET)
            if ( (*slave).empty() == (*surface).empty() )
                logging::error(false, "COUPLING: requires exactly one of SLAVE or SFSET.");

            // MASTER existence & non-empty (use ->size())
            if (!model._data->node_sets.has(*master))
                logging::error(false, "COUPLING: MASTER node set '" + *master + "' does not exist.");
            if (model._data->node_sets.has(*master) && model._data->node_sets.get(*master)->size() == 0)
                logging::error(false, "COUPLING: MASTER node set '" + *master + "' is empty.");

            // If SLAVE is used: existence & non-empty (use ->size())
            if (!(*slave).empty() && !model._data->node_sets.has(*slave))
                logging::error(false, "COUPLING: SLAVE node set '" + *slave + "' does not exist.");
            if (!(*slave).empty() && model._data->node_sets.has(*slave) &&
                model._data->node_sets.get(*slave)->size() == 0)
                logging::error(false, "COUPLING: SLAVE node set '" + *slave + "' is empty.");

            // TYPE validation
            if (*type != "KINEMATIC" && *type != "STRUCTURAL")
                logging::error(false, "COUPLING: unsupported TYPE='" + *type + "'.");
        });


        command.variant(
            fem::dsl::Variant::make()
                .doc("One data line: six DOF flags (1/0) for [Ux, Uy, Uz, Rx, Ry, Rz]. Values > 0 enable coupling.")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .fixed<fem::Precision, 6>()
                                .name("DOF")
                                .desc("Coupled degrees of freedom (1=on, 0=off).")
                                .on_missing(fem::Precision{0})
                                .on_empty  (fem::Precision{0})
                        )
                        .bind([&model, master, type, slave, surface](const std::array<fem::Precision, 6>& dofs_raw) {
                            fem::Dofs mask;
                            for (int i = 0; i < 6; ++i) mask(i) = dofs_raw[i] > fem::Precision{0};

                            if (!mask.any())
                                logging::error(true, "COUPLING: all DOF flags are zero; coupling has no effect.");

                            const bool is_surface = !surface->empty();
                            const std::string& slave_ref = is_surface ? *surface : *slave;

                            const constraint::CouplingType ctype =
                                (*type == "KINEMATIC") ? constraint::CouplingType::KINEMATIC
                                                        : constraint::CouplingType::STRUCTURAL;

                            model.add_coupling(*master, slave_ref, mask, ctype, is_surface);
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
