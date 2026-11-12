/**
 * @file register_coupling.cpp
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

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <memory>
#include <string>
#include <utility>

#include "../../../constraints/types/coupling.h"
#include "../../../model/element/element_structural.h"
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
void register_coupling(fem::io::dsl::Registry& registry, model::Model& model) {
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

    // Common compiled reference and slave-region lookup for COUPLING and RIGID BODY.
    const auto resolve = [&model, ctx](const std::string& master, const std::string& surface,
                                        bool only_nodesets = false) {
        *ctx = Context{};
        ctx->master  = master;
        ctx->surface = surface;
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
            if (!only_nodesets && model._data->surface_sets.has(ctx->surface)) {
                ctx->slave_surfaces = model._data->surface_sets.get(ctx->surface);
                logging::error(ctx->slave_surfaces != nullptr && ctx->slave_surfaces->size() > 0,
                    "COUPLING: slave surface ", ctx->surface, " is empty");
                return;
            }

            logging::error(model._data->node_sets.has(ctx->surface),
                "COUPLING: target ", ctx->surface,
                only_nodesets ? " is not a node set" : " is neither a surface nor a node set");
            ctx->slave_nodes = model._data->node_sets.get(ctx->surface);
            logging::error(ctx->slave_nodes != nullptr && ctx->slave_nodes->size() > 0,
                "COUPLING: slave node set ", ctx->surface, " is empty");
    };

    registry.command("COUPLING", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
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

        command.on_enter([ctx, resolve](const fem::io::dsl::Keys& keys) {
            resolve(keys.raw("MASTER"), keys.raw("SURFACE"));
            ctx->type = keys.has("TYPE") ? keys.raw("TYPE") : std::string{};
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


    // Use the Abaqus rigid-body region selectors. Legacy NSET silently aliases
    // TIE NSET, retaining compatibility with existing FEMaster decks.
    registry.command("RIGIDBODY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Define a rigid body from element and/or PIN/TIE node sets.");
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("REFNODE").required().doc("Reference node or one-node node set")
                .key("ELSET").optional().doc("Rigid elements")
                .key("PINNSET").optional().doc("Translation-only nodes")
                .key("TIENSET").alternative("NSET").optional().doc("Translation and rotation nodes")
                .key("POSITION").optional("INPUT").allowed({"INPUT", "CENTER OF MASS"})
                .key("ROTNODE").optional().doc("CalculiX rotation node (unsupported)")
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            logging::error(!keys.has("ROTNODE"), "RIGIDBODY: ROT NODE is not supported");
            logging::error(model._data->compiled, "RIGIDBODY: model must be compiled");
            logging::error(keys.has("ELSET") || keys.has("PINNSET") || keys.has("TIENSET"),
                "RIGIDBODY: ELSET, PIN NSET or TIE NSET must be supplied");

            const std::string ref = keys.raw("REFNODE");
            model::NodeRegion::Ptr master_region;
            ID master = -1;
            if (model._data->node_sets.has(ref)) {
                master_region = model._data->node_sets.get(ref);
                logging::error(master_region && master_region->size() == 1,
                    "RIGIDBODY: REF NODE set must contain one node: ", ref);
                master = master_region->first();
            } else {
                master = model.compiled_node_id(ref);
            }

            model::RigidBody body;
            body.reference_node = master;
            body.center_of_mass = keys.raw("POSITION") == "CENTER OF MASS";

            // ELSET members contribute their mass and element loads, but not
            // their elastic stiffness. Translational elements imply PIN and
            // elements with rotational DOFs (shell/beam) imply TIE.
            std::map<ID, bool> node_is_tie;
            if (keys.has("ELSET")) {
                const std::string name = keys.raw("ELSET");
                logging::error(model._data->elem_sets.has(name),
                    "RIGIDBODY: unknown ELSET ", name);
                auto region = model._data->elem_sets.get(name);
                logging::error(region && region->size() > 0,
                    "RIGIDBODY: empty ELSET ", name);
                for (ID element_id : *region) {
                    logging::error(element_id >= 0
                        && static_cast<std::size_t>(element_id) < model._data->elements.size()
                        && model._data->elements[static_cast<std::size_t>(element_id)] != nullptr,
                        "RIGIDBODY: invalid element ", element_id);
                    logging::error(!model._data->is_rigid_element(element_id),
                        "RIGIDBODY: element in more than one rigid body: ", element_id);
                    const auto& element = model._data->elements[static_cast<std::size_t>(element_id)];
                    logging::error(element->as<model::StructuralElement>() != nullptr,
                        "RIGIDBODY: non-structural element ", element_id);
                    const auto dofs = element->dofs();
                    const bool tie = dofs(0, 3) || dofs(0, 4) || dofs(0, 5);
                    body.rigid_elements.push_back(element_id);
                    for (ID node : *element) {
                        auto [it, inserted] = node_is_tie.emplace(node, tie);
                        if (!inserted) it->second = it->second || tie;
                    }
                }
            }

            std::set<ID> explicit_pin;
            std::set<ID> explicit_tie;
            auto resolve_set = [&](const char* key, std::set<ID>& ids) {
                if (!keys.has(key)) return;
                const std::string name = keys.raw(key);
                logging::error(model._data->node_sets.has(name),
                    "RIGIDBODY: unknown node set ", name);
                auto region = model._data->node_sets.get(name);
                logging::error(region && region->size() > 0,
                    "RIGIDBODY: empty node set ", name);
                for (ID node : *region) ids.insert(node);
            };
            resolve_set("PINNSET", explicit_pin);
            resolve_set("TIENSET", explicit_tie);

            // Explicit PIN/TIE assignments supersede ELSET-derived types.
            for (ID node : explicit_pin) {
                logging::error(explicit_tie.count(node) == 0,
                    "RIGIDBODY: node occurs in both PIN NSET and TIE NSET: ", node);
                node_is_tie[node] = false;
            }
            for (ID node : explicit_tie) node_is_tie[node] = true;

            for (const auto& [node, tie] : node_is_tie) {
                if (node == master) continue; // No self-coupling equation.
                if (tie) body.tie_nodes.push_back(node);
                else     body.pin_nodes.push_back(node);
            }

            for (const auto& previous : model._data->rigid_bodies) {
                logging::error(previous.reference_node != master,
                    "RIGIDBODY: reference node reused by another rigid body");
                std::set<ID> used(previous.pin_nodes.begin(), previous.pin_nodes.end());
                used.insert(previous.tie_nodes.begin(), previous.tie_nodes.end());
                for (const auto& [node, tie] : node_is_tie) {
                    (void) tie;
                    logging::error(used.count(node) == 0,
                        "RIGIDBODY: node in more than one rigid body: ", node);
                }
            }

            if (body.center_of_mass) {
                Precision total_mass = Precision(0);
                Vec3 first_moment = Vec3::Zero();
                auto integrate_mass = [&](const model::ElementPtr& element) {
                    auto* structural = element->as<model::StructuralElement>();
                    if (!structural) return;
                    const Precision mass = structural->integrate_scalar_field(
                        true, [](const Vec3&) { return Precision(1); });
                    logging::error(std::isfinite(mass) && mass >= Precision(0),
                        "RIGIDBODY: non-finite or negative rigid mass");
                    total_mass += mass;
                    first_moment += structural->integrate_vector_field(
                        true, [](const Vec3& x) { return x; });
                };

                for (ID id : body.rigid_elements) {
                    integrate_mass(model._data->elements[static_cast<std::size_t>(id)]);
                }
                // Include point masses attached to any rigid node or the
                // reference point. Compiled MASS elements are distinct from
                // post-compile POINTMASS elements, and both are supported.
                auto is_rigid_node = [&](ID id) {
                    return id == master || node_is_tie.count(id) != 0;
                };
                for (const auto& element : model._data->elements) {
                    if (!element || element->n_nodes() != 1
                        || !is_rigid_node(element->nodes()[0])) continue;
                    if (std::find(body.rigid_elements.begin(), body.rigid_elements.end(),
                                  element->elem_id) == body.rigid_elements.end()) {
                        integrate_mass(element);
                    }
                }
                for (const auto& point : model._data->point_elements) {
                    if (point && point->n_nodes() == 1
                        && is_rigid_node(point->nodes()[0])) {
                        integrate_mass(point);
                    }
                }
                logging::error(total_mass > Precision(0)
                    && std::isfinite(total_mass) && first_moment.allFinite(),
                    "RIGIDBODY: CENTER OF MASS requires positive physical mass");
                logging::error(model._data->positions != nullptr,
                    "RIGIDBODY: nodal positions not initialized");
                const Vec3 center = first_moment / total_mass;
                for (Index d = 0; d < 3; ++d) {
                    (*model._data->positions)(master, d) = center[d];
                    if (model._data->positions_reference
                        && model._data->positions_reference != model._data->positions) {
                        (*model._data->positions_reference)(master, d) = center[d];
                    }
                }
            }

            // Retain a separate rigid body record; reuse the tested linear
            // kinematic Coupling implementation for each resolved node type.
            auto couple = [&](const std::vector<ID>& ids, const Dofs& dofs,
                              const std::string& suffix) {
                if (ids.empty()) return;
                auto region = std::make_shared<model::NodeRegion>(
                    "RIGIDBODY_" + suffix + "_" + std::to_string(master));
                for (ID id : ids) region->add(id);
                constraint::Coupling coupling{master, region, dofs,
                    constraint::CouplingType::KINEMATIC};
                coupling.master_region = master_region;
                model._data->couplings.emplace_back(std::move(coupling));
            };
            couple(body.pin_nodes, Dofs{true, true, true, false, false, false}, "PIN");
            couple(body.tie_nodes, Dofs{true, true, true, true, true, true}, "TIE");

            if (model._data->rigid_element_mask.size() < model._data->elements.size()) {
                model._data->rigid_element_mask.resize(model._data->elements.size(), 0);
            }
            for (ID id : body.rigid_elements) {
                model._data->rigid_element_mask[static_cast<std::size_t>(id)] = 1;
            }
            model._data->rigid_bodies.push_back(std::move(body));
        });

        command.variant(fem::io::dsl::Variant::make());
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
