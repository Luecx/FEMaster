/**
 * @file register_coupling.inl
 * @brief Registers kinematic and structural couplings.
 *
 * Coupling syntax is resolved completely inside the parser. The parser chooses
 * the node- or surface-based concrete constructor, assigns the semantic master
 * region used for diagnostics/output and stores the finished Coupling directly
 * in ModelData. Model itself does not dispatch constraint types.
 *
 * The data line selects the generalized translational and rotational DOFs that
 * participate in the coupling. Equation generation and distribution over the
 * resolved slave region are deferred to the concrete constraint object.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../../../constraints/types/coupling.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_coupling(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("COUPLING", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));

        auto master  = std::make_shared<std::string>();
        auto slave   = std::make_shared<std::string>();
        auto surface = std::make_shared<std::string>();
        auto type    = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MASTER").required()
                .key("TYPE").required().allowed({"KINEMATIC", "STRUCTURAL"})
                .key("SLAVE").optional()
                .key("SFSET").optional()
        );
        command.on_enter([&model, master, slave, surface, type](const fem::io::dsl::Keys& keys) {
            *master  = keys.raw("MASTER");
            *type    = keys.raw("TYPE");
            *slave   = keys.has("SLAVE") ? keys.raw("SLAVE") : std::string{};
            *surface = keys.has("SFSET") ? keys.raw("SFSET") : std::string{};

            logging::error(model._data->compiled,
                "COUPLING: constraints require a compiled model");
            logging::error(slave->empty() != surface->empty(),
                "COUPLING: exactly one of SLAVE or SFSET is required");
            logging::error(model._data->node_sets.has(*master),
                "COUPLING: master node set ", *master, " does not exist");
            logging::error(model._data->node_sets.get(*master) && model._data->node_sets.get(*master)->size() == 1,
                "COUPLING: master node set ", *master, " must contain exactly one node");

            if (!slave->empty()) {
                logging::error(model._data->node_sets.has(*slave),
                    "COUPLING: slave node set ", *slave, " does not exist");
                logging::error(model._data->node_sets.get(*slave) && model._data->node_sets.get(*slave)->size() > 0,
                    "COUPLING: slave node set ", *slave, " is empty");
            } else {
                logging::error(model._data->surface_sets.has(*surface),
                    "COUPLING: slave surface set ", *surface, " does not exist");
                logging::error(model._data->surface_sets.get(*surface) && model._data->surface_sets.get(*surface)->size() > 0,
                    "COUPLING: slave surface set ", *surface, " is empty");
            }
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 6>().name("DOF")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, master, slave, surface, type](const std::array<fem::Precision, 6>& raw) {
                    Dofs dofs;
                    for (Index i = 0; i < 6; ++i) {
                        dofs(i) = raw[static_cast<std::size_t>(i)] > Precision(0);
                    }

                    const auto coupling_type = *type == "KINEMATIC"
                        ? constraint::CouplingType::KINEMATIC
                        : constraint::CouplingType::STRUCTURAL;
                    const auto master_region = model._data->node_sets.get(*master);
                    const ID master_node = master_region->first();

                    if (!slave->empty()) {
                        constraint::Coupling coupling{
                            master_node,
                            model._data->node_sets.get(*slave),
                            dofs,
                            coupling_type
                        };
                        coupling.master_region = master_region;
                        model._data->couplings.emplace_back(std::move(coupling));
                    } else {
                        constraint::Coupling coupling{
                            master_node,
                            model._data->surface_sets.get(*surface),
                            dofs,
                            coupling_type
                        };
                        coupling.master_region = master_region;
                        model._data->couplings.emplace_back(std::move(coupling));
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
