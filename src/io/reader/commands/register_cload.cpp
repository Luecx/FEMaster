/**
 * @file register_cload.cpp
 * @brief Registers native and Abaqus-compatible concentrated nodal loads.
 *
 * CLOAD accepts FEMaster rows TARGET, Fx, Fy, Fz [, Mx, My, Mz] and
 * Abaqus-compatible rows TARGET, DOF, MAGNITUDE. Both formats materialize the
 * generalized load vector [Fx, Fy, Fz, Mx, My, Mz] in bc::CLoad.
 *
 * Parser::activate_load_collector() selects only an explicit named definition
 * target; unnamed analysis conditions update ModelData::conditions. Nodal
 * TRANSFORM assignments override the command-wide ORIENTATION. Nodes sharing
 * the same effective coordinate system are collected into one load region.
 *
 * Instance-qualified scalar references are mapped into compiled assembly node
 * IDs. Coordinate-system transformation and amplitude evaluation remain in
 * bc::CLoad::apply(), where the load basis is evaluated at each node position.
 *
 * @see bc::CLoad
 * @see Parser::activate_load_collector
 * @see model::Model::resolve_node_region
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>

#include "../parser.h"
#include "../util.h"
#include "../../../bc/structural/load_c.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers force and moment assignments to compiled nodes.
 *
 * Command entry resolves the optional default orientation and amplitude and
 * selects named definition storage only when NAME is supplied. Unnamed analysis
 * Abaqus-format rows update persistent source/DOF history, while unnamed
 * native vector rows apply only in the current analysis. Both formats normalize
 * to six generalized components before regions and nodal bases are resolved.
 *
 * Nodes sharing an effective coordinate system receive one common load. Nodal
 * TRANSFORM assignments override the command-wide orientation, or the global
 * basis if it is omitted. Forces and moments are stored in the selected local
 * basis without changing the model's nodal degrees of freedom.
 *
 * The DSL selects one format for the complete command occurrence. Rows with
 * TARGET, DOF, MAGNITUDE use the Abaqus-compatible form, while native rows
 * require all three force components and may optionally provide three moments.
 *
 * @param registry Parser registry receiving the command definition.
 * @param parser Parser providing the compiled model and nodal transform state.
 */
void register_cload(dsl::Registry& registry, Parser& parser) {
    registry.command("CLOAD", [&](dsl::Command& command) {
        // Loads require compiled nodes and may be declared inside an analysis.
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc(
            "Create concentrated nodal loads using FEMaster vector components or an "
            "Abaqus-compatible DOF and magnitude. NAME is required outside LOADCASE."
        );

        // Keep resolved modifiers for all rows of one occurrence. Command entry
        // resets these handles so omitted modifiers do not leak into later loads.
        auto orientation = std::make_shared<cos::CoordinateSystem::Ptr>(nullptr);
        auto amplitude   = std::make_shared<bc::Amplitude::Ptr        >(nullptr);

        // ---------------------------------------------------------------------
        // Select the collector and resolve command-wide modifiers
        // ---------------------------------------------------------------------
        // NAME selects reusable definitions. Without NAME, native vector
        // rows apply only within their analysis; Abaqus rows update history.
        auto collector = std::make_shared<bc::LoadCollector::Ptr>(nullptr);

        command.keyword(
            dsl::KeywordSpec::make()
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                .key("NAME")
                    .alternative("LOAD_COLLECTOR")
                    .alternative("LOADCOLLECTOR")
                    .optional()
                    .doc("Collector name; optional inside LOADCASE")
                .key("ORIENTATION").optional().doc("Default coordinate system for the load components")
                .key("AMPLITUDE"  ).optional().doc("Amplitude scaling the complete generalized load")
                .flag("FOLLOWER")
                .flag("REAL")
                .flag("IMAGINARY")
        );

        command.on_enter([&parser, orientation, amplitude, collector](const dsl::Keys& keys) {
            auto& model = parser.model();

            orientation->reset();
            amplitude  ->reset();
            collector  ->reset();


            const std::string orientation_name = keys.raw("ORIENTATION");
            const std::string amplitude_name   = keys.raw("AMPLITUDE");
            const std::string collector_name   = keys.raw("NAME");

            const auto* loadcase = parser.active_loadcase();
            logging::error(loadcase == nullptr || loadcase->type_name() != "EIGENFREQ",
                "CLOAD: not supported in a FREQUENCY step");
            logging::error(!keys.has("FOLLOWER"),
                "CLOAD: FOLLOWER is not supported");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "CLOAD: REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "CLOAD: IMAGINARY is not supported");
            logging::error(orientation_name.empty() || model._data->coordinate_systems.has(orientation_name),
                "CLOAD: coordinate system ", orientation_name, " does not exist");
            logging::error(amplitude_name.empty() || model._data->amplitudes.has(amplitude_name),
                "CLOAD: amplitude ", amplitude_name, " does not exist");
            logging::error(!collector_name.empty() || loadcase != nullptr,
                "CLOAD: definitions outside an analysis require NAME");

            if (!orientation_name.empty())
                *orientation = model._data->coordinate_systems.get(orientation_name);
            if (!amplitude_name.empty())
                *amplitude = model._data->amplitudes.get(amplitude_name);
            if (!collector_name.empty()) {
                parser.activate_load_collector(collector_name);
                *collector = model._data->load_cols.get();
            }
            if (collector_name.empty() && keys.raw("OP") == "NEW")
                parser.clear_conditions(bc::CLOAD);
        });

        // ---------------------------------------------------------------------
        // Materialize one logical source/DOF definition across all orientations
        // ---------------------------------------------------------------------

        const auto add_load = [&parser, orientation, amplitude, collector](
            const std::string& target,
            const Vec6&        values,
            bool               persistent
        ) {
            auto& model = parser.model();
            auto region = model.resolve_node_region(target);
            const std::string identifier = (model._data->node_sets.has(target) ? "NSET:" : "NODE:") + target;

            std::vector<bc::Condition::Ptr> replacements;
            for (auto& [load_region, load_orientation] :
                 group_by_orientation(parser, std::move(region), *orientation)) {

                if (*collector || !persistent) {
                    auto load = std::make_shared<bc::CLoad>();
                    load->region_      = std::move(load_region);
                    load->orientation_ = std::move(load_orientation);
                    load->amplitude_   = *amplitude;
                    load->values_      = values;
                    if (*collector) {
                        (*collector)->add(std::move(load));
                    } else {
                        parser.select_collector_condition(bc::CLOAD, std::move(load));
                    }
                    continue;
                }

                for (Dim dof = 0; dof < 6; ++dof) {
                    if (std::isnan(values[dof])) continue;

                    auto load = std::make_shared<bc::CLoad>();
                    load->region_      = load_region;
                    load->orientation_ = load_orientation;
                    load->amplitude_   = *amplitude;
                    load->values_[dof] = values[dof];
                    replacements.push_back(std::move(load));
                }
            }

            if (persistent && !*collector) {
                parser.modify_conditions(bc::CLOAD, identifier, std::move(replacements));
            }
        };

        // ---------------------------------------------------------------------
        // FEMaster format: TARGET, Fx, Fy, Fz [, Mx, My, Mz]
        // ---------------------------------------------------------------------
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string   >().name("TARGET").desc("Compiled node set or scalar node reference")
                    .fixed<Precision, 3>().name("FORCE" ).desc("Fx, Fy, Fz")
                    .fixed<Precision, 3>().name("MOMENT").desc("Mx, My, Mz")
                        .on_missing(Precision{0})
                        .on_empty  (Precision{0})
                )
                .bind([add_load](const std::string&              target,
                                 const std::array<Precision, 3>& force,
                                 const std::array<Precision, 3>& moment) {
                    // Normalize translational and rotational components in
                    // their common selected basis into the generalized vector.
                    Vec6 values;
                    values << force[0], force[1], force[2], moment[0], moment[1], moment[2];
                    add_load(target, values, false);
                })
            )
        );

        // ---------------------------------------------------------------------
        // Abaqus-compatible format: TARGET, DOF, MAGNITUDE
        // ---------------------------------------------------------------------
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET"   ).desc("Compiled node set or scalar node reference")
                    .one<int        >().name("DOF"      ).desc("Loaded degree of freedom, 1 through 6")
                    .one<Precision  >().name("MAGNITUDE").desc("Load magnitude")
                )
                .bind([add_load](const std::string& target, int dof, Precision magnitude) {
                    // Abaqus numbers translations 1..3 and rotations 4..6.
                    // Validate before indexing the zero-based generalized vector.
                    logging::error(dof >= 1 && dof <= 6,
                        "CLOAD: DOF must be in [1,6]");

                    // Preserve the selected DOF as condition identity. NaN marks
                    // components that are not part of this CLOAD definition and
                    // therefore must not be replaced by a later independent DOF.
                    Vec6 values = Vec6::Constant(NAN);
                    values[dof - 1] = magnitude;
                    add_load(target, values, true);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
