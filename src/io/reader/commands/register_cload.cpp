/**
 * @file register_cload.cpp
 * @brief Shared native/Abaqus concentrated load input, collector and STEP handling.
 *
 * The two physical row layouts are distinguished solely by their token count:
 * TARGET,DOF,MAGNITUDE or TARGET,Fx,Fy,Fz[,Mx,My,Mz].
 * One DSL segment normalizes each row independently, so rows may be mixed.
 */
#include "register_functions.h"
#include "../parser.h"
#include "../../dsl/registry.h"
#include "../../dsl/invoke.h"
#include "../../dsl/pattern_element.h"
#include "../../../bc/neumann/load_c.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_harmonic.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/nonlinear_static.h"
#include "../../../model/model.h"

#include <algorithm>
#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::io::reader::commands {
namespace dsl = fem::io::dsl;
namespace {
constexpr const char* omitted = "__CLOAD_OMITTED__";
constexpr const char* empty = "__CLOAD_EMPTY__";

Precision numerical_value(const std::string& token) {
    if (token == empty || token == omitted) {
        throw std::runtime_error("CLOAD: a required numeric value is missing");
    }
    // Parse the entire lexeme: accepting just a numeric prefix would be unsafe
    // for concentrated forces and moments.
    std::size_t end = 0;
    const Precision result = static_cast<Precision>(std::stod(token, &end));
    if (end != token.size()) throw std::runtime_error("CLOAD: invalid number '" + token + "'");
    return result;
}

Vec6 read_components(const std::array<std::string, 6>& components) {
    std::size_t count = 0;
    while (count < components.size() && components[count] != omitted) ++count;
    for (std::size_t i = count; i < components.size(); ++i) {
        if (components[i] != omitted)
            throw std::runtime_error("CLOAD: a component follows an omitted value");
    }

    Vec6 result = Vec6::Zero();
    if (count == 2) {
        const std::string& dof_token = components[0];
        if (!dsl::detail::is_int_token(dof_token))
            throw std::runtime_error("CLOAD: DOF must be an integer in [1,6]");
        const int dof = std::stoi(dof_token);
        if (dof < 1 || dof > 6) throw std::runtime_error("CLOAD: DOF must be in [1,6]");
        result[dof - 1] = numerical_value(components[1]);
    } else if (count >= 3 && count <= 6) {
        for (std::size_t i = 0; i < count; ++i)
            result[static_cast<Index>(i)] = components[i] == empty
                ? Precision(0) : numerical_value(components[i]);
    } else {
        throw std::runtime_error("CLOAD: expected DOF,magnitude or Fx,Fy,Fz[,Mx,My,Mz]");
    }
    return result;
}

std::vector<std::string>* active_loads(loadcase::LoadCase* base) {
    if (auto* lc = dynamic_cast<loadcase::LinearBuckling*>(base)) return &lc->loads;
    if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) return &lc->loads;
    if (auto* lc = dynamic_cast<loadcase::NonlinearStatic*>(base)) return &lc->loads;
    if (auto* lc = dynamic_cast<loadcase::Transient*>(base)) return &lc->loads;
    if (auto* lc = dynamic_cast<loadcase::LinearHarmonic*>(base)) return &lc->loads;
    return nullptr;
}

void add_concentrated_load(model::Model& model, model::NodeRegion::Ptr region, const Vec6& values,
                           cos::CoordinateSystem::Ptr orientation, bc::Amplitude::Ptr amplitude) {
    auto load = std::make_shared<bc::CLoad>();
    load->region_ = std::move(region);
    load->values_ = values;
    load->orientation_ = std::move(orientation);
    load->amplitude_ = std::move(amplitude);
    model.add_load(std::move(load));
}
} // namespace

void register_cload(dsl::Registry& registry, Parser& parser) {
    registry.command("CLOAD", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STEP"}));
        command.doc("Apply nodal loads as TARGET,DOF,MAGNITUDE or TARGET,Fx,Fy,Fz[,Mx,My,Mz]. "
                    "LOAD_COLLECTOR is required outside a LOADCASE or STEP.");

        auto orientation = std::make_shared<cos::CoordinateSystem::Ptr>(nullptr);
        auto amplitude = std::make_shared<bc::Amplitude::Ptr>(nullptr);
        auto scale = std::make_shared<Precision>(Precision(1));
        auto in_step = std::make_shared<bool>(false);

        command.keyword(
            dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").alternative("LOADCOLLECTOR").alternative("NAME").optional()
                .key("ORIENTATION").optional()
                .key("AMPLITUDE").optional()
                .flag("REAL").flag("IMAGINARY").flag("FOLLOWER")
        );
        command.on_enter([&parser, orientation, amplitude, scale, in_step](
                             const dsl::ParentInfo& parent, const dsl::Keys& keys) {
            auto& model = parser.model();
            const bool in_loadcase = parent.command == "LOADCASE";
            *in_step = parent.command == "STEP";
            const bool in_analysis = in_loadcase || *in_step;
            const std::string requested = keys.raw("LOAD_COLLECTOR");
            logging::error(in_analysis || !requested.empty(),
                "CLOAD: LOAD_COLLECTOR is required outside LOADCASE/STEP");
            logging::error(!keys.has("FOLLOWER"), "CLOAD: FOLLOWER is unsupported");
            logging::error(!keys.has("IMAGINARY"), "CLOAD: IMAGINARY is unsupported");

            if (*in_step) {
                logging::error(parser.step_state().step_active && parser.active_loadcase() != nullptr,
                    "CLOAD: an analysis procedure must precede loads inside STEP");
                logging::error(parser.active_loadcase()->type_name() != "EIGENFREQ",
                    "CLOAD: not supported inside a FREQUENCY step");
            }

            orientation->reset();
            amplitude->reset();
            *scale = Precision(1);
            const std::string basis = keys.raw("ORIENTATION");
            if (!basis.empty()) {
                logging::error(model._data->coordinate_systems.has(basis),
                    "CLOAD: coordinate system ", basis, " does not exist");
                *orientation = model._data->coordinate_systems.get(basis);
            }

            std::string amplitude_name = keys.raw("AMPLITUDE");
            if (*in_step) {
                const auto resolved = parser.resolve_load_amplitude(amplitude_name);
                *scale = resolved.first;
                amplitude_name = resolved.second;
            }
            if (!amplitude_name.empty()) {
                logging::error(model._data->amplitudes.has(amplitude_name),
                    "CLOAD: amplitude ", amplitude_name, " does not exist");
                *amplitude = model._data->amplitudes.get(amplitude_name);
            }

            std::string collector = requested;
            if (in_analysis) {
                auto* lc = parser.active_loadcase();
                logging::error(lc != nullptr, "CLOAD: no active loadcase");
                auto* names = active_loads(lc);
                logging::error(names != nullptr, "CLOAD: unsupported analysis type");
                if (collector.empty()) {
                    collector = *in_step ? "__ABQ_STEP_LOADS"
                        : "__FEMASTER_INLINE_CLOAD_" + std::to_string(lc->id);
                }
                if (std::find(names->begin(), names->end(), collector) == names->end())
                    names->push_back(collector);
            }
            model._data->load_cols.activate(collector);
        });

        command.on_exit([&parser, in_step](const dsl::Keys&) {
            // Other STEP load commands rely on the default collector being active.
            if (*in_step) parser.model()._data->load_cols.activate("__ABQ_STEP_LOADS");
        });

        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .fixed<std::string, 6>().name("VALUES")
                        .on_missing(std::string{omitted}).on_empty(std::string{empty})
                )
                .bind([&parser, orientation, amplitude, scale](
                          const std::string& target, const std::array<std::string, 6>& components) {
                    auto& model = parser.model();
                    const Vec6 values = read_components(components) * *scale;
                    const bool has_nodal_transform = !parser.step_state().node_transforms.empty();

                    const auto add_node = [&](ID node_id) {
                        auto basis = *orientation;
                        if (!basis && has_nodal_transform) {
                            const auto& transforms = parser.step_state().node_transforms;
                            const auto it = transforms.find(node_id);
                            if (it != transforms.end()) {
                                logging::error(model._data->coordinate_systems.has(it->second),
                                    "CLOAD: coordinate system ", it->second, " does not exist");
                                basis = model._data->coordinate_systems.get(it->second);
                            }
                        }
                        auto region = std::make_shared<model::NodeRegion>("INTERNAL");
                        region->add(node_id);
                        add_concentrated_load(model, std::move(region), values, std::move(basis), *amplitude);
                    };

                    if (model._data->node_sets.has(target)) {
                        // Retain the native shared-region representation unless nodes
                        // have Abaqus-specific, potentially different local transforms.
                        if (!has_nodal_transform || *orientation) {
                            add_concentrated_load(model, model._data->node_sets.get(target),
                                                  values, *orientation, *amplitude);
                        } else {
                            for (const ID node : *model._data->node_sets.get(target)) add_node(node);
                        }
                    } else {
                        add_node(model.compiled_node_id(target));
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
