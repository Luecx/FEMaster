/**
 * @file register_coupling.inl
 * @brief Registers kinematic and structural couplings against named node or surface regions.
 *
 * Coupling syntax is resolved completely inside the parser. `SURFACE` is the
 * canonical slave-target keyword, while the legacy FEMaster spellings `SLAVE`
 * and `SFSET` are accepted as aliases. All three names therefore share one
 * target-resolution path: a matching compiled `SurfaceRegion` is preferred and
 * a matching `NodeRegion` is used otherwise. This allows element-based and
 * node-based `*SURFACE` definitions to be consumed through the same coupling
 * syntax while preserving existing FEMaster decks.
 *
 * `MASTER` remains the canonical reference-node keyword and accepts the Abaqus
 * spelling `REF NODE` through the normalized `REFNODE` alias. The reference may
 * identify either a compiled node set containing exactly one node or one direct
 * node reference. The parser assigns the semantic master region when available
 * for diagnostics/output and stores the finished Coupling directly in ModelData;
 * Model itself does not dispatch constraint types.
 *
 * The data line selects the generalized translational and rotational DOFs that
 * participate in the coupling. Equation generation and structural load
 * distribution over the resolved slave region are deferred to the concrete
 * constraint object.
 *
 * @see constraint::Coupling
 * @see model::NodeRegion
 * @see model::SurfaceRegion
 *
 * @author Finn Eggers
 * @date 24.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../reference.h"
#include "../../../constraints/types/coupling.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers the FEMaster coupling grammar and resolves its compiled regions.
 *
 * Keyword aliases are normalized by the DSL before the callbacks execute, so
 * `SURFACE`, `SLAVE` and `SFSET` arrive as one canonical target string and
 * `MASTER`/`REF NODE` arrive as one canonical master string. Target resolution
 * deliberately checks compiled surface sets before node sets. If both namespaces
 * contain the same name, the geometric surface therefore wins deterministically.
 *
 * The master may be a named node set containing exactly one node or a direct
 * compiled node reference such as `42` or `INSTANCE.42`. Slave targets must be
 * non-empty. The final data line converts six positive values into the active
 * coupling-DOF mask before constructing the existing kinematic or structural
 * Coupling implementation.
 *
 * @param registry DSL registry receiving the `COUPLING` command.
 * @param model Compiled model whose regions and constraints are modified.
 */
inline void register_coupling(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("COUPLING", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));

        struct Context {
            std::string master;
            std::string target;
            std::string type;
            model::NodeRegion::Ptr master_region = nullptr;
            ID master_node = -1;
            bool target_is_surface = false;
        };
        auto ctx = std::make_shared<Context>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MASTER").alternative("REFNODE").required()
                .key("TYPE").required().allowed({"KINEMATIC", "STRUCTURAL"})
                .key("SURFACE").alternative("SLAVE").alternative("SFSET").required()
        );
        command.on_enter([&model, ctx](const fem::io::dsl::Keys& keys) {
            ctx->master            = keys.raw("MASTER");
            ctx->target            = keys.raw("SURFACE");
            ctx->type              = keys.raw("TYPE");
            ctx->master_region     = nullptr;
            ctx->master_node       = -1;
            ctx->target_is_surface = false;

            // Coupling equations and target regions operate in compiled assembly space
            logging::error(model._data->compiled,
                "COUPLING: constraints require a compiled model");

            // Resolve the reference node from either a one-node set or a direct
            // assembly node reference such as ID or INSTANCE.ID
            if (model._data->node_sets.has(ctx->master)) {
                ctx->master_region = model._data->node_sets.get(ctx->master);
                logging::error(ctx->master_region != nullptr && ctx->master_region->size() == 1,
                    "COUPLING: master node set ", ctx->master, " must contain exactly one node");
                ctx->master_node = ctx->master_region->first();
            } else {
                ctx->master_node = io::reader::compiled_node_id(model, ctx->master);
            }

            // Resolve one common target namespace. Geometric surfaces take
            // precedence when the same name also exists as a node set.
            if (model._data->surface_sets.has(ctx->target)) {
                const auto target = model._data->surface_sets.get(ctx->target);
                logging::error(target != nullptr && target->size() > 0,
                    "COUPLING: slave surface ", ctx->target, " is empty");
                ctx->target_is_surface = true;
                return;
            }

            logging::error(model._data->node_sets.has(ctx->target),
                "COUPLING: target ", ctx->target, " is neither a surface nor a node set");
            const auto target = model._data->node_sets.get(ctx->target);
            logging::error(target != nullptr && target->size() > 0,
                "COUPLING: slave node set ", ctx->target, " is empty");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 6>().name("DOF")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, ctx](const std::array<fem::Precision, 6>& raw) {
                    // Convert the six deck values into the generalized coupling mask
                    Dofs dofs;
                    for (Index i = 0; i < 6; ++i) {
                        dofs(i) = raw[static_cast<std::size_t>(i)] > Precision(0);
                    }

                    // Keep the existing coupling mechanics unchanged; this parser
                    // change only normalizes names and target-region resolution.
                    const auto coupling_type = ctx->type == "KINEMATIC"
                        ? constraint::CouplingType::KINEMATIC
                        : constraint::CouplingType::STRUCTURAL;

                    if (ctx->target_is_surface) {
                        constraint::Coupling coupling{
                            ctx->master_node,
                            model._data->surface_sets.get(ctx->target),
                            dofs,
                            coupling_type
                        };
                        coupling.master_region = ctx->master_region;
                        model._data->couplings.emplace_back(std::move(coupling));
                        return;
                    }

                    constraint::Coupling coupling{
                        ctx->master_node,
                        model._data->node_sets.get(ctx->target),
                        dofs,
                        coupling_type
                    };
                    coupling.master_region = ctx->master_region;
                    model._data->couplings.emplace_back(std::move(coupling));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
