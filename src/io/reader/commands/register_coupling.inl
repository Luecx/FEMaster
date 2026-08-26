/**
 * @file register_coupling.inl
 * @brief Registers direct and scoped kinematic/structural coupling definitions.
 *
 * `COUPLING` owns the common reference-node and slave-region definition. The
 * slave target is addressed through the canonical `SURFACE` keyword while the
 * legacy FEMaster spellings `SLAVE` and `SFSET` remain aliases. Resolution first
 * checks compiled surface regions and falls back to compiled node regions, so
 * element-based and node-based `*SURFACE` definitions share one coupling path.
 *
 * `MASTER` is the canonical reference-node keyword and accepts the normalized
 * Abaqus spelling `REFNODE` (`REF NODE`) as an alias. The reference may name a
 * one-node `NodeRegion` or one direct compiled node reference such as `42` or
 * `INSTANCE.42`.
 *
 * Two equivalent syntax styles are supported. Existing FEMaster decks may keep
 * `TYPE=KINEMATIC|STRUCTURAL` on `*COUPLING` followed by the six-value DOF mask.
 * If `TYPE` is omitted, `*COUPLING` becomes a scope and exactly one child command
 * must select the behavior:
 *
 * @code
 * *COUPLING, MASTER=RP, SURFACE=CLOUD
 * *KINEMATIC
 * 1, 3
 * 6,
 * @endcode
 *
 * or
 *
 * @code
 * *COUPLING, MASTER=RP, SURFACE=CLOUD
 * *DISTRIBUTING, COUPLING=STRUCTURAL
 * @endcode
 *
 * `KINEMATIC` and `DISTRIBUTING` share the same coupling context. DOF data uses
 * Abaqus-style inclusive `first,last` ranges; omitting all range lines selects
 * all six DOFs and omitting `last` selects only `first`. `DISTRIBUTING` maps
 * deliberately to FEMaster `CouplingType::STRUCTURAL`. Its optional Abaqus
 * `COUPLING=CONTINUUM|STRUCTURAL` selector is accepted, but both values use the
 * FEMaster structural load-distribution implementation by design.
 *
 * Constraint equations and structural load distribution remain implemented by
 * `constraint::Coupling`; this file only resolves deck syntax and model regions.
 *
 * @see constraint::Coupling
 * @see constraint::CouplingType
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
#include <utility>

#include "../../../constraints/types/coupling.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers the coupling parent scope together with its kinematic and distributing
 * child commands.
 *
 * One shared context is created for the three DSL commands. Entering a new
 * `COUPLING` resets and resolves that context. A direct `TYPE=...` definition
 * constructs the constraint from the legacy six-value DOF mask, whereas a parent
 * without `TYPE` consumes no data itself and leaves construction to exactly one
 * `KINEMATIC` or `DISTRIBUTING` child. The parent `on_exit` hook validates that a
 * complete coupling was emitted before the scope closes.
 *
 * Surface targets take precedence over node sets when both namespaces contain the
 * same name. The resolved region pointers are retained in the shared context so
 * child commands do not repeat target lookup or depend on parser-global state.
 *
 * @param registry DSL registry receiving `COUPLING`, `KINEMATIC` and `DISTRIBUTING`.
 * @param model Compiled model whose regions are resolved and whose coupling list is modified.
 */
inline void register_coupling(fem::io::dsl::Registry& registry, model::Model& model) {
    struct Context {
        std::string master;
        std::string surface;
        std::string type;

        ID master_node = -1;

        model::NodeRegion::Ptr    master_region  = nullptr;
        model::NodeRegion::Ptr    slave_nodes    = nullptr;
        model::SurfaceRegion::Ptr slave_surfaces = nullptr;

        bool defined = false;
    };

    auto ctx = std::make_shared<Context>();

    // Construct exactly one Coupling from the already resolved parent context.
    // Both legacy TYPE syntax and child scopes use this path.
    const auto emit = [&model, ctx](const Dofs& dofs, constraint::CouplingType type) {
        logging::error(!ctx->defined,
            "COUPLING: coupling behavior is already defined");
        logging::error(ctx->master_node >= 0,
            "COUPLING: reference node is not resolved");
        logging::error(ctx->slave_surfaces != nullptr || ctx->slave_nodes != nullptr,
            "COUPLING: slave region is not resolved");

        if (ctx->slave_surfaces) {
            constraint::Coupling coupling{ctx->master_node, ctx->slave_surfaces, dofs, type};
            coupling.master_region = ctx->master_region;
            model._data->couplings.emplace_back(std::move(coupling));
        } else {
            constraint::Coupling coupling{ctx->master_node, ctx->slave_nodes, dofs, type};
            coupling.master_region = ctx->master_region;
            model._data->couplings.emplace_back(std::move(coupling));
        }

        ctx->defined = true;
    };

    // Apply one inclusive Abaqus-style DOF range to a child coupling mask. The
    // first explicit range replaces the default "all six" selection.
    const auto add_dof_range = [](Dofs& dofs,
                                  bool& has_range,
                                  Index first,
                                  Index last,
                                  const char* command) {
        logging::error(first >= 1 && first <= 6,
            command, ": first DOF must be between 1 and 6");
        if (last >= static_cast<Index>(-1)) last = first;
        logging::error(last >= first && last <= 6,
            command, ": last DOF must be between first DOF and 6");

        if (!has_range) {
            dofs = Dofs{false, false, false, false, false, false};
            has_range = true;
        }
        for (Index dof = first; dof <= last; ++dof) {
            dofs(dof - 1) = true;
        }
    };

    registry.command("COUPLING", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT", "ASSEMBLY"));
        command.doc("Define a coupling reference node and slave surface or node set.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MASTER").alternative("REFNODE").required()
                    .doc("Reference node or one-node node set")
                .key("SURFACE").alternative("SLAVE").alternative("SFSET").required()
                    .doc("Slave surface or node set")
                .key("TYPE").optional().allowed({"KINEMATIC", "STRUCTURAL"})
                    .doc("Legacy direct coupling type; omit to use a child scope")
        );

        command.on_enter([&model, ctx](const fem::io::dsl::Keys& keys) {
            // Every parent definition owns a fresh context; child commands retain
            // this exact object until the parent scope closes.
            *ctx = Context{};
            ctx->master  = keys.raw("MASTER");
            ctx->surface = keys.raw("SURFACE");
            ctx->type    = keys.has("TYPE") ? keys.raw("TYPE") : std::string{};

            logging::error(model._data->compiled,
                "COUPLING: constraints require a compiled model");

            // Resolve the reference node from either a one-node set or a direct
            // compiled node reference such as ID or INSTANCE.ID.
            if (model._data->node_sets.has(ctx->master)) {
                ctx->master_region = model._data->node_sets.get(ctx->master);
                logging::error(ctx->master_region != nullptr && ctx->master_region->size() == 1,
                    "COUPLING: master node set ", ctx->master, " must contain exactly one node");
                ctx->master_node = ctx->master_region->first();
            } else {
                ctx->master_node = model.compiled_node_id(ctx->master);
            }

            // Resolve the common slave target once. A geometric surface wins if
            // the same name is also present in the node-set namespace.
            if (model._data->surface_sets.has(ctx->surface)) {
                ctx->slave_surfaces = model._data->surface_sets.get(ctx->surface);
                logging::error(ctx->slave_surfaces != nullptr && ctx->slave_surfaces->size() > 0,
                    "COUPLING: slave surface ", ctx->surface, " is empty");
                return;
            }

            logging::error(model._data->node_sets.has(ctx->surface),
                "COUPLING: target ", ctx->surface, " is neither a surface nor a node set");
            ctx->slave_nodes = model._data->node_sets.get(ctx->surface);
            logging::error(ctx->slave_nodes != nullptr && ctx->slave_nodes->size() > 0,
                "COUPLING: slave node set ", ctx->surface, " is empty");
        });

        command.on_exit([ctx](const fem::io::dsl::Keys&) {
            logging::error(ctx->defined,
                "COUPLING: requires TYPE or one KINEMATIC/DISTRIBUTING child");
        });

        // Legacy FEMaster syntax: TYPE is specified on the parent and one line
        // carries the six boolean-like generalized DOF selectors.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_present("TYPE"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 6>().name("DOF")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([ctx, emit](const std::array<fem::Precision, 6>& raw) {
                    Dofs dofs;
                    for (Index i = 0; i < 6; ++i) {
                        dofs(i) = raw[static_cast<std::size_t>(i)] > Precision(0);
                    }

                    const auto type = ctx->type == "KINEMATIC"
                        ? constraint::CouplingType::KINEMATIC
                        : constraint::CouplingType::STRUCTURAL;
                    emit(dofs, type);
                })
            )
        );

        // Scoped syntax: the parent intentionally consumes no data. The next
        // admitted child command selects the coupling behavior and DOFs.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::negate(
                fem::io::dsl::Condition::key_present("TYPE")
            ))
        );
    });

    // `KINEMATIC` is valid only as a child of a COUPLING without direct TYPE.
    // Multiple range lines can select non-contiguous DOFs; no lines means all six.
    auto kinematic_dofs      = std::make_shared<Dofs>();
    auto kinematic_has_range = std::make_shared<bool>(false);

    registry.command("KINEMATIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::all_of({
            fem::io::dsl::Condition::parent_is("COUPLING"),
            fem::io::dsl::Condition::negate(
                fem::io::dsl::Condition::parent_has_key("TYPE")
            )
        }));
        command.doc("Select a kinematic coupling and optional inclusive DOF ranges.");

        command.on_enter([ctx, kinematic_dofs, kinematic_has_range](const fem::io::dsl::Keys&) {
            logging::error(!ctx->defined,
                "KINEMATIC: coupling behavior is already defined");
            *kinematic_dofs      = Dofs{true, true, true, true, true, true};
            *kinematic_has_range = false;
        });

        command.on_exit([kinematic_dofs, emit](const fem::io::dsl::Keys&) {
            emit(*kinematic_dofs, constraint::CouplingType::KINEMATIC);
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Index>().name("FIRST")
                    .one<Index>().name("LAST").on_missing(static_cast<Index>(-1)).on_empty(static_cast<Index>(-1))
                )
                .bind([kinematic_dofs, kinematic_has_range, add_dof_range](Index first, Index last) {
                    add_dof_range(*kinematic_dofs, *kinematic_has_range, first, last, "KINEMATIC");
                })
            )
        );
    });

    // `DISTRIBUTING` deliberately maps to FEMaster STRUCTURAL behavior. The
    // Abaqus continuum/structural selector is syntax-compatible but does not
    // change the internal implementation.
    auto distributing_dofs      = std::make_shared<Dofs>();
    auto distributing_has_range = std::make_shared<bool>(false);

    registry.command("DISTRIBUTING", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::all_of({
            fem::io::dsl::Condition::parent_is("COUPLING"),
            fem::io::dsl::Condition::negate(
                fem::io::dsl::Condition::parent_has_key("TYPE")
            )
        }));
        command.doc("Select FEMaster structural load distribution for a coupling.");
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("COUPLING").optional("CONTINUUM").allowed({"CONTINUUM", "STRUCTURAL"})
                    .doc("Accepted Abaqus distributing subtype; both map to STRUCTURAL")
        );

        command.on_enter([ctx, distributing_dofs, distributing_has_range](const fem::io::dsl::Keys&) {
            logging::error(!ctx->defined,
                "DISTRIBUTING: coupling behavior is already defined");
            *distributing_dofs      = Dofs{true, true, true, true, true, true};
            *distributing_has_range = false;
        });

        command.on_exit([distributing_dofs, emit](const fem::io::dsl::Keys&) {
            emit(*distributing_dofs, constraint::CouplingType::STRUCTURAL);
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Index>().name("FIRST")
                    .one<Index>().name("LAST").on_missing(static_cast<Index>(-1)).on_empty(static_cast<Index>(-1))
                )
                .bind([distributing_dofs, distributing_has_range, add_dof_range](Index first, Index last) {
                    add_dof_range(*distributing_dofs, *distributing_has_range, first, last, "DISTRIBUTING");
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
