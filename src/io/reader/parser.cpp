/**
 * @file parser.cpp
 * @brief Implements one-shot parsing and input-condition history resolution.
 *
 * The input file is parsed exactly once into `dsl::Deck`. Semantic model
 * construction then follows a visible dependency order: global definitions,
 * Part-local topology, Instances, `Model::compile()`, assembly materialization,
 * compiled model features and finally load cases. Scope-sensitive commands are
 * selected from their stored parent nodes instead of being replayed in parser
 * stages.
 *
 * The registered command callbacks remain the semantic implementation.
 * `ParsedCommand::enter()`, `execute()` and `leave()` control when those existing
 * callbacks run; syntax validation, variant selection and data aggregation have
 * already completed in `dsl::DeckParser`.
 *
 * In addition to syntax processing, this file owns the input-side condition
 * history. Two generations of original input identifiers may each address
 * several physical Condition pointers after a nodal TRANSFORM split. The reader
 * updates compatible values in place and retires incompatible fragments;
 * ConditionManager stores only active pointer membership per family.
 *
 * @see Parser
 * @see io::dsl::Deck
 * @see io::dsl::DeckParser
 * @see model::Model::compile
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#include "parser.h"

#include "../../bc/structural/load_c.h"
#include "../../bc/structural/load_d.h"
#include "../../bc/structural/load_p.h"
#include "../../bc/structural/load_v.h"
#include "../../bc/structural/load_inertial.h"
#include "../../bc/structural/support.h"
#include "../../bc/thermal/temperature.h"
#include "../../bc/thermal/heat_flux.h"
#include "../../bc/thermal/convection.h"
#include "../../core/logging.h"
#include "../../loadcase/loadcase.h"
#include "../../loadcase/linear_transient.h"
#include "../../loadcase/nonlinear_static.h"
#include "../../model/model.h"
#include "../dsl/deck_parser.h"
#include "../dsl/file.h"
#include "../writer/writers.h"
#include "commands/register_functions.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

namespace fem::io::reader {
namespace {

/**
 * @brief Prepares the start/target history of a physical load in the DSL.
 *
 * Existing loads are carried into the next step with start=target. New
 * definitions enter with start=0 and the prescribed target. MOD or OP=NEW
 * retires an inherited definition by setting only its target to zero; for an
 * explicit amplitude the start is its effective value at the step start.
 * Its original start value must be retained because it is still needed
 * during the current step. This is especially important when a definition
 * was introduced and replaced again within the same input step.
 *
 * Unspecified NaN vector components retain their omission mask. Supports and
 * prescribed temperatures are not additive loads and intentionally do not
 * follow this procedure; releasing a constrained DOF requires separate
 * treatment of the corresponding constraint equation.
 *
 * @param condition Physical definition whose numeric history is prepared.
 * @param newly_defined Whether the current step introduces the definition.
 * @param retired Whether the definition is being ramped out of the step.
 */
void prepare_load_history(bc::Condition& condition, bool newly_defined, bool retired,
                          Precision start_time = Precision(0)) {
    const Precision amplitude_at_start = retired && condition.amplitude_
        ? condition.amplitude_->evaluate(start_time) : Precision(1);
    const bool nonlinear_inertia = dynamic_cast<bc::InertialLoad*>(&condition) != nullptr;

    const auto prepare_value = [&](Precision& target, Precision& start) {
        if (std::isnan(target)) return;

        if (retired) {
            if (condition.amplitude_ && !nonlinear_inertia) start = amplitude_at_start * target;
            else if (std::isnan(start)) start = target;
            target = Precision(0);
        } else {
            start = newly_defined ? Precision(0) : target;
        }
    };

    const auto prepare_vector = [&](auto& target, auto& start) {
        for (Eigen::Index i = 0; i < target.size(); ++i) {
            prepare_value(target[i], start[i]);
        }
    };

    if (auto* load = dynamic_cast<bc::CLoad*>(&condition)) {
        prepare_vector(load->values_, load->values_start_);
    } else if (auto* load = dynamic_cast<bc::DLoad*>(&condition)) {
        prepare_vector(load->values_, load->values_start_);
    } else if (auto* load = dynamic_cast<bc::VLoad*>(&condition)) {
        prepare_vector(load->values_, load->values_start_);
    } else if (auto* load = dynamic_cast<bc::PLoad*>(&condition)) {
        prepare_value(load->pressure_, load->pressure_start_);
    } else if (auto* load = dynamic_cast<bc::InertialLoad*>(&condition)) {
        if (retired && condition.amplitude_) load->start_scale_ = amplitude_at_start;
        if (!retired) load->start_scale_ = Precision(1);
        prepare_vector(load->center_acc_, load->center_acc_start_);
        prepare_vector(load->omega_,      load->omega_start_);
        prepare_vector(load->alpha_,      load->alpha_start_);
    } else if (auto* load = dynamic_cast<bc::HeatFlux*>(&condition)) {
        prepare_value(load->heat_flux_, load->heat_flux_start_);
    } else {
        return;
    }

    // A retired definition must decay from its previous effective value
    // even if the original condition used a named amplitude. New definitions
    // can independently install a different amplitude on their own object.
    if (retired) condition.amplitude_.reset();
}

// Evaluate explicitly prescribed amplitude history at the start of this step.
Precision step_start_time(const loadcase::LoadCase* active) {
    const auto* transient = dynamic_cast<const loadcase::Transient*>(active);
    return transient ? static_cast<Precision>(transient->t_start) : Precision(0);
}

bool supports_load_transition(bc::ConditionFamily family) {
    return family == bc::CLOAD || family == bc::DLOAD || family == bc::DSLOAD
        || family == bc::PLOAD || family == bc::VLOAD
        || family == bc::INERTIAL_LOAD || family == bc::HEAT_FLUX;
}

/**
 * A same-step load introduced at zero may be discarded if MOD supersedes it.
 * A carried-over load retains a nonzero step-start contribution and must ramp
 * down even when its target or orientation changes.
 */
bool has_nonzero_load_start(const bc::Condition& condition) {
    const auto nonzero = [](const auto& values) {
        for (Eigen::Index i = 0; i < values.size(); ++i) {
            if (std::isfinite(values[i]) && values[i] != Precision(0)) return true;
        }
        return false;
    };
    if (const auto* a = dynamic_cast<const bc::CLoad*>(&condition)) return nonzero(a->values_start_);
    if (const auto* a = dynamic_cast<const bc::DLoad*>(&condition)) return nonzero(a->values_start_);
    if (const auto* a = dynamic_cast<const bc::VLoad*>(&condition)) return nonzero(a->values_start_);
    if (const auto* a = dynamic_cast<const bc::PLoad*>(&condition)) {
        return std::isfinite(a->pressure_start_) && a->pressure_start_ != Precision(0);
    }
    if (const auto* a = dynamic_cast<const bc::InertialLoad*>(&condition)) {
        return nonzero(a->center_acc_start_) || nonzero(a->omega_start_) || nonzero(a->alpha_start_);
    }
    if (const auto* a = dynamic_cast<const bc::HeatFlux*>(&condition)) {
        return std::isfinite(a->heat_flux_start_) && a->heat_flux_start_ != Precision(0);
    }
    return false;
}

bool has_zero_load_target(const bc::Condition& condition) {
    const auto zero_vector = [](const auto& values) {
        for (Eigen::Index i = 0; i < values.size(); ++i) {
            if (std::isfinite(values[i]) && values[i] != Precision(0)) return false;
        }
        return true;
    };

    if (const auto* load = dynamic_cast<const bc::CLoad*>(&condition)) return zero_vector(load->values_);
    if (const auto* load = dynamic_cast<const bc::DLoad*>(&condition)) return zero_vector(load->values_);
    if (const auto* load = dynamic_cast<const bc::VLoad*>(&condition)) return zero_vector(load->values_);
    if (const auto* load = dynamic_cast<const bc::PLoad*>(&condition)) return load->pressure_ == Precision(0);
    if (const auto* load = dynamic_cast<const bc::InertialLoad*>(&condition)) {
        return load->center_acc_.isZero() && load->omega_.isZero() && load->alpha_.isZero();
    }
    if (const auto* load = dynamic_cast<const bc::HeatFlux*>(&condition)) return load->heat_flux_ == Precision(0);
    return false;
}

} // namespace


/**
 * Constructs an idle parser and prepares the unified command documentation grammar.
 */
Parser::Parser()
    : model_(std::make_shared<model::Model>()),
      writer_("") {
    configure_documentation_registry();
}

Parser::~Parser() = default;

/**
 * Parses one complete input deck and executes its semantic commands afterwards.
 *
 * Every run owns one registry for the entire parse/process lifetime. This is
 * required because parsed command and segment pointers refer directly to that
 * registry. The file itself is consumed only once; all later dependency ordering
 * operates on the in-memory deck representation.
 *
 * @param input_path Input deck parsed exactly once.
 * @param output_path Optional base path for result files.
 * @param writer_formats Result formats enabled for analysis output.
 */
void Parser::run(const std::string& input_path,
                 const std::string& output_path,
                 const io::writer::WriterFileFormats& writer_formats) {
    // Reset all mutable state so each run represents an independent deck.
    model_ = std::make_shared<model::Model>();
    active_loadcase_.reset();
    next_loadcase_id_ = 1;
    step_state_ = StepState{};
    node_transforms.clear();
    for (auto& index : inherited_conditions_) index.clear();
    for (auto& index : current_conditions_) index.clear();
    selected_conditions_.clear();

    // Register the complete grammar once and parse the complete source once.
    io::dsl::Registry registry;
    register_commands(registry);

    io::dsl::File       file(input_path);
    io::dsl::DeckParser deck_parser(registry);
    const io::dsl::Deck deck = deck_parser.parse(file);

    // Semantic dependencies are now explicit and independent of source order.
    process_deck(deck, input_path, output_path, writer_formats);
}

/**
 * Executes the unified FEMaster/Abaqus deck in explicit model-dependency order.
 *
 * The function deliberately spells out the semantic order from top to bottom.
 * Loops only iterate repeated input scopes such as MATERIAL, PART, ASSEMBLY,
 * COUPLING and LOADCASE; command ordering inside each scope remains visible at
 * the call site. `Model::compile()` forms the one-way boundary between sparse
 * Part/Instance topology and dense assembly materialization.
 */
void Parser::process_deck(const io::dsl::Deck&                  deck,
                          const std::string&                    input_path,
                          const std::string&                    output_path,
                          const io::writer::WriterFileFormats& writer_formats) {
    const auto& root = deck.root();

    // ---------------------------------------------------------------------
    // Global definitions required by sections and later load definitions
    // ---------------------------------------------------------------------
    root.execute_children("HEADING");

    for (const auto* material : root.children("MATERIAL")) {
        material->enter();

        material->execute_children("ELASTIC");
        material->execute_children("HYPERELASTIC");
        material->execute_children("PLASTIC");
        material->execute_children("DENSITY");
        material->execute_children("THERMALEXPANSION");
        material->execute_children("EXPANSION");

        material->leave();
    }

    root.execute_children("PROFILE");
    root.execute_children("ORIENTATION");
    root.execute_children("AMPLITUDE");

    // ---------------------------------------------------------------------
    // Default-Part topology before Model::compile()
    // ---------------------------------------------------------------------
    root.execute_children("NODE");
    root.execute_children("ELEMENT");

    root.execute_children("NSET");
    root.execute_children("ELSET");
    root.execute_children("SURFACE");
    root.execute_children("SFSET");

    root.execute_children("SOLIDSECTION");
    root.execute_children("BEAMSECTION");
    root.execute_children("TRUSSSECTION");
    root.execute_children("SHELLSECTION");

    root.execute_children("MASS");
    root.execute_children("ROTARYINERTIA");
    root.execute_children("SPRING");

    // ---------------------------------------------------------------------
    // Explicit Part topology before Model::compile()
    // ---------------------------------------------------------------------
    for (const auto* part : root.children("PART")) {
        part->enter();

        part->execute_children("NODE");
        part->execute_children("ELEMENT");

        part->execute_children("NSET");
        part->execute_children("ELSET");
        part->execute_children("SURFACE");
        part->execute_children("SFSET");

        part->execute_children("SOLIDSECTION");
        part->execute_children("BEAMSECTION");
        part->execute_children("TRUSSSECTION");
        part->execute_children("SHELLSECTION");

        part->execute_children("MASS");
        part->execute_children("ROTARYINERTIA");
        part->execute_children("SPRING");

        part->leave();
    }

    // ---------------------------------------------------------------------
    // Assembly orphan topology before Model::compile()
    // ---------------------------------------------------------------------
    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("NODE");
        assembly->execute_children("ELEMENT");
    }

    // ---------------------------------------------------------------------
    // Instances depend on completed Parts and orphan assembly topology
    // ---------------------------------------------------------------------
    root.execute_children("INSTANCE");

    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("INSTANCE");
    }

    // ---------------------------------------------------------------------
    // Transition from sparse topology to dense assembly data
    // ---------------------------------------------------------------------
    model_->compile();

    // ---------------------------------------------------------------------
    // Assembly regions, properties and dense fields
    // ---------------------------------------------------------------------
    root.execute_children("FIELD");
    root.execute_children("NORMAL");

    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("NSET");
        assembly->execute_children("ELSET");
        assembly->execute_children("SURFACE");
        assembly->execute_children("SFSET");

        assembly->execute_children("MASS");
        assembly->execute_children("ROTARYINERTIA");
        assembly->execute_children("SPRING");

        assembly->execute_children("FIELD");
        assembly->execute_children("NORMAL");
    }

    // Initial conditions bind already materialized named fields to persistent
    // model-state handles. Temperature and velocity then remain authoritative
    // until a later analysis procedure explicitly evolves them.
    root.execute_children("INITIALCONDITION");
    root.execute_children("INITIALCONDITIONS");

    // Apply initial tie adjustments before geometry-derived reference fields are completed.
    root.execute_children("TIE");
    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("TIE");
    }

    // Complete reference normals from the final initial geometry.
    model_->build_shell_element_normals();

    // ---------------------------------------------------------------------
    // Compiled model features
    // ---------------------------------------------------------------------
    root.execute_children("POINTMASS");

    // Resolve all compiled nodal bases before any CLOAD or SUPPORT is materialized.
    // Assembly sets already exist, so these assignments are independent of
    // whether their loads are declared at root, assembly or load-case scope.
    root.execute_children("TRANSFORM");
    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("TRANSFORM");
    }

    // ---------------------------------------------------------------------
    // Root-level topology-dependent constraints
    // ---------------------------------------------------------------------
    root.execute_children("RBM");
    root.execute_children("CONNECTOR");
    root.execute_children("MPC");
    root.execute_children("CONTACT");
    root.execute_children("EQUATION");

    // ---------------------------------------------------------------------
    // Assembly-level topology-dependent constraints
    // ---------------------------------------------------------------------
    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("RBM");
        assembly->execute_children("CONNECTOR");
        assembly->execute_children("MPC");
        assembly->execute_children("CONTACT");
        assembly->execute_children("EQUATION");
    }

    // ---------------------------------------------------------------------
    // Couplings
    // ---------------------------------------------------------------------
    for (const auto* coupling : root.children("COUPLING")) {
        coupling->enter();

        coupling->execute_children("KINEMATIC");
        coupling->execute_children("DISTRIBUTING");

        coupling->leave();
    }

    for (const auto* assembly : root.children("ASSEMBLY")) {
        for (const auto* coupling : assembly->children("COUPLING")) {
            coupling->enter();

            coupling->execute_children("KINEMATIC");
            coupling->execute_children("DISTRIBUTING");

            coupling->leave();
        }
    }

    // ---------------------------------------------------------------------
    // Result writers, condition definitions and sequential analysis execution
    // ---------------------------------------------------------------------
    initialize_writers(input_path, output_path, writer_formats);

    // Topology and geometric constraints are complete. Execute only condition
    // definitions and analyses in source order so later root-level history or
    // named definitions cannot change the state of an earlier solve.
    auto is_condition_command = [](const std::string& name) {
        return name == "BOUNDARY" || name == "SUPPORT" || name == "CLOAD" || name == "DLOAD"
            || name == "DSLOAD" || name == "PLOAD" || name == "VLOAD"
            || name == "INERTIALOAD" || name == "TEMPERATURE"
            || name == "HEATFLUX" || name == "CONVECTION";
    };

    for (const auto* command : root.children()) {
        const auto& name = command->command().name_;
        if (is_condition_command(name)) {
            command->execute();
        } else if (name == "ASSEMBLY") {
            // Assembly-level condition targets already resolve compiled regions
            for (const auto* child : command->children()) {
                if (is_condition_command(child->command().name_)) {
                    child->execute();
                }
            }
        } else if (name == "OVERVIEW") {
            // Report the model after preceding named definitions have been read
            command->execute();
        } else if (name == "LOADCASE") {
            // The native loadcase builds its analysis directly on entering.
            command->enter();
            command->execute_children();
            command->leave();
        } else if (name == "STEP") {
            // STEP is purely a scope. Exactly one procedure constructs the analysis,
            // and its final ENDSTEP child executes it after all conditions/output.
            const auto children = command->children();
            logging::error(children.size() == 1,
                "STEP requires exactly one analysis procedure");

            const auto* procedure = children.front();
            const auto& procedure_name = procedure->command().name_;
            logging::error(procedure_name == "STATIC" || procedure_name == "FREQUENCY"
                        || procedure_name == "BUCKLE" || procedure_name == "DYNAMIC"
                        || procedure_name == "STEADYSTATEDYNAMICS",
                "STEP: unsupported analysis procedure ", procedure_name);

            const auto entries = procedure->children();
            logging::error(!entries.empty() && entries.back()->command().name_ == "ENDSTEP",
                "STEP requires END STEP after the procedure commands");

            command->enter();
            procedure->enter();
            for (const auto* entry : entries) entry->execute();
            procedure->leave();
            command->leave();
        }
    }

    close_writers();
}

/**
 * Initializes enabled result writers from the requested output path and publishes
 * the compiled model topology before any load case starts writing frames.
 */
void Parser::initialize_writers(const std::string&                    input_path,
                                const std::string&                    output_path,
                                const io::writer::WriterFileFormats& writer_formats) {
    std::string writer_base = output_path.empty() ? input_path : output_path;
    for (const std::string& ext : { std::string(".res"),  std::string(".frd"),
                                    std::string(".femr"), std::string(".fil"), std::string(".inp")}) {
        if (writer_base.size() >= ext.size()
         && writer_base.compare(writer_base.size() - ext.size(), ext.size(), ext) == 0) {
            writer_base.resize (writer_base.size() - ext.size());
            break;
        }
    }

    writer_ = io::writer::ResultWriters(writer_base, writer_formats);
    writer_.write_model_data(*model_->_data);
}

/**
 * Flushes and closes every enabled result writer after semantic analysis processing.
 */
void Parser::close_writers() {
    writer_.close();
}

/**
 * Prints one requested view of the registered FEMaster command grammar.
 */
void Parser::document(const DocOptions& opts) const {
    using A = DocOptions::Action;
    using F = DocOptions::Format;
    using V = DocOptions::Verbosity;

    if (opts.format != F::Text) {
        std::cout << "(Note) Only text output is implemented currently. Falling back to text.\n\n";
    }

    switch (opts.action) {
        case A::List:       documentation_registry_.print_index(); break;
        case A::Show:       documentation_registry_.print_help(opts.cmd, opts.verbosity == V::Compact); break;
        case A::Tokens:     documentation_registry_.print_tokens(opts.cmd); break;
        case A::Variants:   documentation_registry_.print_variants(opts.cmd); break;
        case A::Search:     documentation_registry_.print_search(opts.query, opts.regex); break;
        case A::WhereToken: documentation_registry_.print_where_token(opts.query); break;
        case A::All:        documentation_registry_.print_help({}, false); break;
    }
}

/**
 * Returns the mutable model currently constructed or analyzed by the parser.
 */
model::Model& Parser::model() {
    logging::error(model_ != nullptr,
        "Parser: model is not initialized");
    return *model_;
}

/**
 * Returns the model currently constructed or analyzed by the parser.
 */
const model::Model& Parser::model() const {
    logging::error(model_ != nullptr,
        "Parser: model is not initialized");
    return *model_;
}

const io::dsl::Registry& Parser::registry() const { return documentation_registry_; }

/**
 * Activates a load case and supplies its parser-owned analysis context.
 */
void Parser::begin_loadcase(loadcase::LoadCase::Ptr loadcase) {
    logging::error(active_loadcase_ == nullptr,
        "Parser: nested load cases are not supported");
    logging::error(loadcase != nullptr,
        "Parser: cannot activate a null load case");
    logging::error(model_ != nullptr,
        "Parser: cannot activate a load case without a model");

    loadcase->id     = next_loadcase_id_++;
    loadcase->writer = &writer_;
    loadcase->model  = model_.get();

    // Output requests are configured while the step is parsed, but derived
    // fields are recovered later during run(). Bind the model before exposing
    // the active load case to any child output command.
    loadcase->output.bind(model_.get());

    // At entry all existing input definitions are inherited. Current-step
    // modifications move their references into current_conditions_.
    for (auto& index : current_conditions_) index.clear();

    // For the first analysis, load definitions enter from zero. Later
    // analyses inherit the previous target as their new starting load.
    // All updates and removals within the step preserve this starting level.
    const bool first_step = loadcase->id == 1;
    for (const auto family : {bc::CLOAD, bc::DLOAD, bc::DSLOAD, bc::PLOAD,
                              bc::VLOAD, bc::INERTIAL_LOAD, bc::HEAT_FLUX}) {
        for (const auto& condition : conditions().get(family)) {
            prepare_load_history(*condition, first_step, false);
        }
    }

    // Named collector selection remains local to this analysis.
    selected_conditions_.clear();
    active_loadcase_ = std::move(loadcase);
}

/**
 * @brief Runs one analysis with direct and temporarily selected conditions.
 *
 * Direct input definitions already reside in ModelData's ConditionManager.
 * Reusable named collector definitions are inserted into the corresponding
 * family only for the duration of this run. The insertion Boolean identifies
 * pointers that were not already active, so cleanup never removes an existing
 * direct condition that happens to share a collector object.
 *
 * On normal return, the reader removes just those newly inserted pointers.
 * The reusable collector keeps its shared definitions. Active loads whose
 * end magnitude has reached zero are then removed from the manager, and their
 * input identifier entries are pruned so subsequent steps cannot match them.
 *
 * The solver obtains start/target interpolation from the concrete condition
 * implementations; this method only updates history ownership after success.
 */
void Parser::end_loadcase() {
    logging::error(active_loadcase_ != nullptr,
        "Parser: cannot end a load case when none is active");

    auto loadcase = std::move(active_loadcase_);

    // Insert only missing collector pointers and remember precisely those
    // references. Previously active direct conditions must survive cleanup.
    std::vector<std::pair<bc::ConditionFamily, bc::Condition::Ptr>> temporary_conditions;
    for (const auto& [family, condition] : selected_conditions_) {
        if (conditions().add(family, condition)) {
            temporary_conditions.emplace_back(family, condition);
        }
    }

    // Model assembly reads the active family sets during this numerical solve.
    loadcase->run();

    // Restore direct history by removing only this analysis's insertions.
    for (const auto& [family, condition] : temporary_conditions) {
        conditions().remove(family, condition);
    }
    selected_conditions_.clear();

    // Physically retired loads have a zero target. Their references remained
    // active through the entire step, so remove them only after the solve
    // returns successfully. Nonzero targets form the carried-over history.
    for (const auto family : {bc::CLOAD, bc::DLOAD, bc::DSLOAD, bc::PLOAD,
                              bc::VLOAD, bc::INERTIAL_LOAD, bc::HEAT_FLUX}) {
        std::vector<bc::Condition::Ptr> retired;
        for (const auto& condition : conditions().get(family)) {
            if (has_zero_load_target(*condition)) retired.push_back(condition);
        }
        for (const auto& condition : retired) conditions().remove(family, condition);

        // Remove retired pointers from both logical generations.
        for (auto* index : {&inherited_conditions_[family], &current_conditions_[family]}) {
            for (auto it = index->begin(); it != index->end();) {
                auto& entries = it->second;
                entries.erase(std::remove_if(entries.begin(), entries.end(), [&](const auto& condition) {
                    return !conditions().contains(family, condition);
                }), entries.end());
                if (entries.empty()) it = index->erase(it);
                else ++it;
            }
        }
    }

    // Commit the modified definitions only after the solve completed.
    for (std::size_t i = 0; i < bc::N_CONDITION_FAMILIES; ++i) {
        auto& inherited = inherited_conditions_[i];
        auto& current   = current_conditions_[i];
        for (auto& [identifier, entries] : current) {
            auto& target = inherited[identifier];
            target.insert(target.end(), entries.begin(), entries.end());
        }
        current.clear();
    }
}

/**
 * Returns the load case currently configured by consecutive parser commands.
 */
loadcase::LoadCase* Parser::active_loadcase() {
    return active_loadcase_.get();
}

StepState& Parser::step_state() {
    return step_state_;
}

const StepState& Parser::step_state() const {
    return step_state_;
}

/**
 * Forwards access to model-owned active condition storage. Identifier lookup
 * and replacement semantics are maintained independently by the input reader.
 */
bc::ConditionManager& Parser::conditions() {
    return model_->_data->conditions;
}

/**
 * Forwards read-only access to the current ModelData condition history.
 */
const bc::ConditionManager& Parser::conditions() const {
    return model_->_data->conditions;
}

/**
 * @brief Inserts one physical condition and indexes its original input target.
 *
 * The ConditionManager contains the pointer in one independent family. A
 * newly inserted pointer with a non-empty input identifier is additionally
 * stored in the reader index. Several physical pointers may have the same
 * identifier because one input definition can be split across TRANSFORM bases.
 *
 * Repeated insertion of an identical pointer has no effect on either the
 * active set or the identifier index.
 *
 * @param family Input-history scope receiving the physical definition.
 * @param condition Non-null shared physical definition.
 * @param identifier Original input target, retained only by the parser.
 */
void Parser::add_condition(bc::ConditionFamily family, bc::Condition::Ptr condition,
                           const std::string& identifier) {
    // Validate before inserting and indexing a physical definition.
    logging::error(condition != nullptr,
        "Cannot register a null condition");

    // A load first defined inside an analysis ramps from zero unless a named
    // amplitude explicitly prescribes its temporal history. Definitions made
    // outside a step keep their legacy nominal target semantics.
    if (active_loadcase_) prepare_load_history(*condition, true, false);

    // Even unnamed definitions need a logical entry for OP=NEW.
    if (conditions().add(family, condition)) {
        auto& index = active_loadcase_ ? current_conditions_ : inherited_conditions_;
        index[family][identifier].push_back(std::move(condition));
    }
}

/**
 * The input reader compares semantic source and DOF identity independently
 * of physical region splitting. Compatible conditions retain their pointers
 * and step-start values; incompatible conditions preserve the old contribution
 * for a step-long ramp-down.
 */
// Logical target comparison remains in the input layer.
static bool same_condition_target(const bc::Condition& current, const bc::Condition& replacement,
                                  const std::string& identifier) {
    // Input-history identity belongs to the reader, not to the physical
    // Condition interface. Regions are equal if they reference the same
    // compiled region or contain the same ordered entity IDs.
    const auto same_region = [](const auto& left, const auto& right) {
        return left && right && (left == right || left->data() == right->data());
    };

    // NaN denotes an omitted generalized component rather than a prescribed
    // zero. Different masks therefore occupy different history targets.
    const auto same_mask = [](const auto& left, const auto& right) {
        for (Eigen::Index i = 0; i < left.size(); ++i) {
            if (std::isnan(left[i]) != std::isnan(right[i])) return false;
        }
        return true;
    };

    // Compare only target identity and active component masks. Numerical
    // magnitudes, orientations and amplitudes are replacement values.
    // The concrete type is deliberately resolved here in the input layer;
    // assembly never needs to know these keyword replacement rules.
    const auto same_target = [&](const bc::Condition& current, const bc::Condition& replacement) {
        if (const auto* a = dynamic_cast<const bc::CLoad*>(&current)) {
            const auto* b = dynamic_cast<const bc::CLoad*>(&replacement);
            if (!b) return false;

            // For transformed Abaqus CLOADs the original identifier already
            // selected the candidate fragments. Their individual node regions
            // intentionally differ; only the generalized DOF mask must match.
            return (identifier.empty() ? same_region(a->region_, b->region_) : true)
                && same_mask(a->values_, b->values_);
        }

        if (const auto* a = dynamic_cast<const bc::DLoad*>(&current)) {
            const auto* b = dynamic_cast<const bc::DLoad*>(&replacement);
            return b && same_region(a->region_, b->region_)
                     && same_mask(a->values_, b->values_);
        }
        if (const auto* a = dynamic_cast<const bc::PLoad*>(&current)) {
            const auto* b = dynamic_cast<const bc::PLoad*>(&replacement);
            return b && same_region(a->region_, b->region_);
        }
        if (const auto* a = dynamic_cast<const bc::VLoad*>(&current)) {
            const auto* b = dynamic_cast<const bc::VLoad*>(&replacement);
            return b && same_region(a->region_, b->region_)
                     && same_mask(a->values_, b->values_);
        }
        if (const auto* a = dynamic_cast<const bc::InertialLoad*>(&current)) {
            const auto* b = dynamic_cast<const bc::InertialLoad*>(&replacement);
            return b && same_region(a->region_, b->region_)
                     && a->consider_point_masses_ == b->consider_point_masses_;
        }
        if (const auto* a = dynamic_cast<const bc::Support*>(&current)) {
            const auto* b = dynamic_cast<const bc::Support*>(&replacement);
            if (!b) return false;

            // Identified native SUPPORT rows may split into multiple nodal
            // TRANSFORM fragments. The original source selects the fragments;
            // only the prescribed generalized DOF mask must still match.
            const bool same_target_region = !identifier.empty()
                || same_region(a->node_region(),    b->node_region())
                || same_region(a->element_region(), b->element_region())
                || same_region(a->surface_region(), b->surface_region());
            return same_target_region && same_mask(a->values(), b->values());
        }
        if (const auto* a = dynamic_cast<const bc::Temperature*>(&current)) {
            const auto* b = dynamic_cast<const bc::Temperature*>(&replacement);
            return b
                && (same_region(a->node_region_,    b->node_region_)
                 || same_region(a->element_region_, b->element_region_)
                 || same_region(a->surface_region_, b->surface_region_));
        }
        if (const auto* a = dynamic_cast<const bc::HeatFlux*>(&current)) {
            const auto* b = dynamic_cast<const bc::HeatFlux*>(&replacement);
            return b && same_region(a->region_, b->region_);
        }
        if (const auto* a = dynamic_cast<const bc::Convection*>(&current)) {
            const auto* b = dynamic_cast<const bc::Convection*>(&replacement);
            return b && same_region(a->region_, b->region_);
        }
        return false;
    };

    return same_target(current, replacement);
}

// Update a compatible existing definition without changing its start values.
// Supports and incompatible orientations deliberately fall back to replacement.
static bool update_load_target(bc::Condition& current, const bc::Condition& incoming) {
    if (current.amplitude_ != incoming.amplitude_) return false;
    const auto same_region = [](const auto& a, const auto& b) {
        return a && b && (a == b || a->data() == b->data());
    };
    if (auto* a = dynamic_cast<bc::CLoad*>(&current)) {
        const auto* b = dynamic_cast<const bc::CLoad*>(&incoming);
        if (!b || !same_region(a->region_, b->region_) || a->orientation_ != b->orientation_) return false;
        a->values_ = b->values_;
        return true;
    }
    if (auto* a = dynamic_cast<bc::DLoad*>(&current)) {
        const auto* b = dynamic_cast<const bc::DLoad*>(&incoming);
        if (!b || !same_region(a->region_, b->region_) || a->orientation_ != b->orientation_) return false;
        a->values_ = b->values_;
        return true;
    }
    if (auto* a = dynamic_cast<bc::VLoad*>(&current)) {
        const auto* b = dynamic_cast<const bc::VLoad*>(&incoming);
        if (!b || !same_region(a->region_, b->region_) || a->orientation_ != b->orientation_) return false;
        a->values_ = b->values_;
        return true;
    }
    if (auto* a = dynamic_cast<bc::PLoad*>(&current)) {
        const auto* b = dynamic_cast<const bc::PLoad*>(&incoming);
        if (!b || !same_region(a->region_, b->region_)) return false;
        a->pressure_ = b->pressure_;
        return true;
    }
    if (auto* a = dynamic_cast<bc::InertialLoad*>(&current)) {
        const auto* b = dynamic_cast<const bc::InertialLoad*>(&incoming);
        if (!b || !same_region(a->region_, b->region_)
               || a->consider_point_masses_ != b->consider_point_masses_
               || !a->center_.isApprox(b->center_)) return false;
        a->center_acc_ = b->center_acc_;
        a->omega_      = b->omega_;
        a->alpha_      = b->alpha_;
        return true;
    }
    if (auto* a = dynamic_cast<bc::HeatFlux*>(&current)) {
        const auto* b = dynamic_cast<const bc::HeatFlux*>(&incoming);
        if (!b || !same_region(a->region_, b->region_)) return false;
        a->heat_flux_ = b->heat_flux_;
        return true;
    }
    return false;
}

/**
 * Handles one input row's complete set of physical fragments atomically.
 * An earlier fragment cannot accidentally replace a later fragment of the
 * same NSET after nodal TRANSFORM splitting.
 */
void Parser::modify_conditions(bc::ConditionFamily family, const std::string& identifier,
                               std::vector<bc::Condition::Ptr> replacements) {
    for (const auto& condition : replacements) {
        logging::error(condition != nullptr, "Cannot modify a null condition");
    }
    if (replacements.empty()) return;

    std::vector<bc::Condition::Ptr> previous;
    for (auto* index : {&inherited_conditions_[family], &current_conditions_[family]}) {
        const auto found = index->find(identifier);
        if (found == index->end()) continue;

        auto& entries = found->second;
        entries.erase(std::remove_if(entries.begin(), entries.end(), [&](const auto& existing) {
            const bool matched = std::any_of(replacements.begin(), replacements.end(), [&](const auto& incoming) {
                return same_condition_target(*existing, *incoming, identifier);
            });
            if (matched) previous.push_back(existing);
            return matched;
        }), entries.end());
        if (entries.empty()) index->erase(found);
    }

    for (auto& incoming : replacements) {
        const auto found = std::find_if(previous.begin(), previous.end(), [&](const auto& existing) {
            return same_condition_target(*existing, *incoming, identifier)
                && update_load_target(*existing, *incoming);
        });
        if (found == previous.end()) {
            add_condition(family, std::move(incoming), identifier);
            continue;
        }
        auto& index = active_loadcase_ ? current_conditions_ : inherited_conditions_;
        index[family][identifier].push_back(*found);
        previous.erase(found);
    }

    for (const auto& existing : previous) {
        if (!active_loadcase_ || !supports_load_transition(family)) {
            conditions().remove(family, existing);
            continue;
        }

        // The previous condition may use an explicit amplitude whose
        // effective start value differs from its nominal value. Evaluate
        // retirement before deciding whether the old fragment is disposable.
        prepare_load_history(*existing, false, true, step_start_time(active_loadcase_.get()));
        if (!has_nonzero_load_start(*existing)) {
            conditions().remove(family, existing);
        }
    }
}

/**
 * Retires only the source definitions inherited at step entry.
 * Current-step MOD and new definitions remain valid. Physical additive loads
 * ramp to zero; constraints are removed immediately. Root-level NEW clears
 * all previously registered direct definitions of this family.
 */
void Parser::clear_conditions(bc::ConditionFamily family) {
    auto& inherited = inherited_conditions_[family];
    for (auto& [identifier, entries] : inherited) {
        for (const auto& condition : entries) {
            if (active_loadcase_ && supports_load_transition(family)) {
                prepare_load_history(*condition, false, true, step_start_time(active_loadcase_.get()));
            } else {
                conditions().remove(family, condition);
            }
        }
    }
    inherited.clear();
}

/**
 * @brief Schedules a named collector definition for the current analysis.
 *
 * A collector may be selected repeatedly or share a physical pointer with
 * another collector. The reader therefore deduplicates the (family, pointer)
 * pair without merging different objects at the same physical target.
 *
 * Selection does not yet change the manager's active sets. end_loadcase()
 * inserts these references immediately before analysis and removes only its
 * new insertions after normal completion.
 *
 * @param family Family in which the selected condition will participate.
 * @param condition Non-null reusable physical definition.
 */
void Parser::select_collector_condition(bc::ConditionFamily family, bc::Condition::Ptr condition) {
    // Reject invalid collector definitions before recording active selection.
    logging::error(condition != nullptr,
        "Cannot select a null condition");

    // Named load definitions may have been parsed outside any analysis. Reject
    // implicit follower pressure when it is actually selected for a nonlinear
    // solve, regardless of the input keyword or collector's stored family.
    if (active_loadcase_ && active_loadcase_->type_name() == "NONLINEARSTATIC") {
        logging::error(dynamic_cast<bc::PLoad*>(condition.get()) == nullptr,
            "LOADS: follower pressure is not supported in nonlinear steps");
    }

    const auto selected = std::make_pair(family, condition);
    if (std::find(selected_conditions_.begin(), selected_conditions_.end(), selected)
        == selected_conditions_.end()) {
        selected_conditions_.push_back(std::move(selected));
    }
}

/**
 * Selects the explicit named target for a reusable structural load definition.
 *
 * Only a non-empty user-supplied name activates collector storage. Direct
 * analysis conditions bypass this operation and update ModelData::conditions.
 * Activation for a solve belongs to LOADS and is independent of this insertion
 * target; defining a named collector does not automatically apply its entries.
 *
 * @param name Non-empty name of the collector being defined.
 * @return Name of the activated definition target.
 */
std::string Parser::activate_load_collector(const std::string& name) {
    // Reject unnamed definition targets before creating persistent model storage
    logging::error(!name.empty(),
        "Named load collector definitions require NAME");

    // Activate only the reusable definition target; solve activation belongs to the DSL
    model()._data->load_cols.activate(name);
    return name;
}

/**
 * Rebuilds the persistent registry used for command-language documentation.
 */
void Parser::configure_documentation_registry() {
    documentation_registry_ = io::dsl::Registry{};
    register_commands(documentation_registry_);
}

/**
 * Registers the unified FEMaster/Abaqus command grammar exactly once per run.
 *
 * Registration describes syntax and semantic callbacks only. No parser stage or
 * command activation mode is attached to the grammar; execution timing belongs
 * exclusively to `process_deck()`.
 */
void Parser::register_commands(io::dsl::Registry& registry) {
    logging::error(model_ != nullptr,
        "Parser: model must exist before registering commands");

    auto& mdl = *model_;

    // Semantic Part/Instance topology and scope-aware assembly commands
    commands::register_part        (registry, mdl);
    commands::register_end_part    (registry, mdl);
    commands::register_assembly    (registry);
    commands::register_end_assembly(registry);
    commands::register_instance    (registry, mdl);
    commands::register_end_instance(registry);
    commands::register_node        (registry, mdl);
    commands::register_element     (registry, mdl);
    commands::register_nset        (registry, mdl);
    commands::register_elset       (registry, mdl);
    commands::register_surface     (registry, mdl);
    commands::register_sfset       (registry, mdl);

    // Root model and native load-case terminator scopes
    registry.command("MODEL", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is("ROOT"));
        command.keyword(io::dsl::KeywordSpec::make().key("NAME").optional());
        command.variant(io::dsl::Variant::make());
    });
    registry.command("END", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is("LOADCASE"));
        command.closes_parent();
        command.variant(io::dsl::Variant::make());
    });

    // Field, material, profile and section definitions
    commands::register_heading(registry);
    commands::register_field(registry, mdl);
    commands::register_normal(registry, mdl);
    commands::register_material(registry, mdl);
    commands::register_elastic(registry, mdl);
    commands::register_hyperelastic(registry, mdl);
    commands::register_plastic(registry, mdl);
    commands::register_density(registry, mdl);
    commands::register_expansion(registry, mdl);
    commands::register_orientation(registry, mdl);
    commands::register_profile(registry, mdl);
    commands::register_solid_section(registry, mdl);
    commands::register_beam_section(registry, mdl);
    commands::register_truss_section(registry, mdl);
    commands::register_shell_section(registry, mdl);
    commands::register_mass(registry, mdl);
    commands::register_rotary_inertia(registry, mdl);
    commands::register_spring(registry, mdl);

    // Loads, constraints, features and model diagnostics
    commands::register_transform(registry, *this);
    commands::register_cload(registry, *this);
    commands::register_dload(registry, *this);
    commands::register_dsload(registry, *this);
    commands::register_pload(registry, *this);
    commands::register_vload(registry, *this);
    commands::register_inertialload(registry, *this);
    commands::register_initial_condition(registry, mdl);
    commands::register_thermal_conditions(registry, *this);
    commands::register_rbm(registry, mdl);
    commands::register_support(registry, *this);
    commands::register_boundary(registry, *this);
    commands::register_amplitude(registry, mdl);
    commands::register_connector(registry, mdl);
    commands::register_mpc(registry, mdl);
    commands::register_coupling(registry, mdl);
    commands::register_tie(registry, mdl);
    commands::register_contact(registry, mdl);
    commands::register_point_mass(registry, mdl);
    commands::register_overview(registry, mdl);
    commands::register_equation(registry, mdl);

    // Load-case creation, solver settings and result requests
    commands::register_loadcase_begin(registry, *this);
    commands::register_step(registry, *this);
    commands::register_loadcase_supports(registry, *this);
    commands::register_loadcase_loads(registry, *this);
    commands::register_loadcase_solver(registry, *this);
    commands::register_loadcase_constraintmethod(registry, *this);
    commands::register_loadcase_frequency(registry, *this);
    commands::register_loadcase_request_stiffness(registry, *this);
    commands::register_loadcase_request_stgeom(registry, *this);
    commands::register_loadcase_numeigenvalues(registry, *this);
    commands::register_loadcase_sigma(registry, *this);
    commands::register_loadcase_topodensity(registry, *this);
    commands::register_loadcase_topoorient(registry, *this);
    commands::register_loadcase_topoexponent(registry, *this);
    commands::register_loadcase_constraintsummary(registry, *this);
    commands::register_loadcase_nonlinear(registry, *this);
    commands::register_loadcase_time(registry, *this);
    commands::register_loadcase_write_every(registry, *this);
    commands::register_loadcase_damping(registry, *this);
    commands::register_loadcase_newmark(registry, *this);
    commands::register_loadcase_initialvelocity(registry, *this);
    commands::register_loadcase_inertiarelief(registry, *this);
    commands::register_loadcase_rebalance(registry, *this);
    commands::register_output(registry, *this);
}

} // namespace fem::io::reader
