/**
 * @file model.cpp
 * @brief Implements model-level registration, lifecycle and compact reporting.
 *
 * This file contains the small behavioral operations that act directly on
 * `ModelData`: nonlinear step-cache coordination, registration of loads,
 * supports and amplitudes, construction of point-mass features and the compact
 * stream representation of a model.
 *
 * Part/instance flattening remains in `model_compile.cpp`, global matrix and
 * load construction remain in `model_build.cpp`, and the detailed hierarchical
 * diagnostic report is implemented in `model_overview.cpp`.
 *
 * @see Model
 * @see ModelData
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "model.h"

#include "../bc/load_collector.h"
#include "../bc/support_collector.h"
#include "../feature/point_mass.h"

#include <iterator>

namespace fem::model {

void Model::step_begin() {
    try {
        for (auto& elem : _data->elements) {
            if (elem) {
                if (auto* structural = elem->as<StructuralElement>()) {
                    structural->step_begin();
                }
            }
        }
    } catch (...) {
        step_end();
        throw;
    }
}

void Model::step_end() {
    for (auto& elem : _data->elements) {
        if (elem) {
            if (auto* structural = elem->as<StructuralElement>()) {
                structural->step_end();
            }
        }
    }
}

void Model::add_load(bc::Load::Ptr load) {
    logging::error(_data != nullptr && _data->compiled,
        "Model: loads require a compiled model");
    logging::error(load != nullptr,
        "Model: cannot add a null load");
    logging::error(_data->load_cols.has_any() && _data->load_cols.get() != nullptr,
        "Model: no load collector is active");

    _data->load_cols.get()->add(std::move(load));
}

void Model::add_amplitude(bc::Amplitude::Ptr amplitude) {
    logging::error(amplitude != nullptr,
        "Model: cannot add a null amplitude");
    logging::error(!_data->amplitudes.has(amplitude->name),
        "Model: amplitude ", amplitude->name, " is already defined");

    _data->amplitudes.add(std::move(amplitude));
}

void Model::add_support(bc::Support support) {
    logging::error(_data != nullptr && _data->compiled,
        "Model: supports require a compiled model");
    logging::error(_data->supp_cols.has_any() && _data->supp_cols.get() != nullptr,
        "Model: no support collector is active");

    _data->supp_cols.get()->add(std::move(support));
}

void Model::add_point_mass_feature(const std::string& nset,
                                   Precision mass,
                                   Vec3 rotary_inertia,
                                   Vec3 spring,
                                   Vec3 rotary_spring) {
    logging::error(_data != nullptr && _data->compiled,
        "Model: point-mass features require a compiled model");
    logging::error(_data->node_sets.has(nset),
        "Node set ", nset, " is not a defined node set");

    auto feature = std::make_shared<feature::PointMass>();
    feature->region_                  = _data->node_sets.get(nset);
    feature->mass_                    = mass;
    feature->rotary_inertia_          = rotary_inertia;
    feature->spring_constants_        = spring;
    feature->rotary_spring_constants_ = rotary_spring;
    _data->features.push_back(std::move(feature));
}

/**
 * Writes a compact model summary to the supplied stream.
 *
 * The representation contains only stable high-level counts and the topology
 * compilation state. Detailed hierarchical diagnostics remain the
 * responsibility of `Model::print_overview()` and are not emitted through the
 * logger as a side effect of streaming.
 *
 * @param ostream Destination stream receiving the model summary.
 * @param model Model whose semantic and compiled entity counts are reported.
 * @return Reference to `ostream` for chained stream operations.
 */
std::ostream& operator<<(std::ostream& ostream, const Model& model) {
    const auto& model_data = *model._data;

    // Count named semantic containers through their public iterator interface.
    // This avoids depending on the associative storage selected by Dict.
    auto entry_count = [](const auto& entries) {
        return static_cast<Index>(std::distance(entries.begin(), entries.end()));
    };

    // Use the compiled position field as the authoritative dense node count
    const Index node_count = model_data.positions ? model_data.positions->rows : Index(0);

    // Write a compact and side-effect-free summary with aligned labels
    ostream << "compiled  = " << (model_data.compiled ? "true" : "false") << '\n';
    ostream << "parts     = " << entry_count(model_data.parts)              << '\n';
    ostream << "instances = " << entry_count(model_data.instances)          << '\n';
    ostream << "nodes     = " << node_count                                 << '\n';
    ostream << "elements  = " << model_data.elements.size()                 << '\n';
    ostream << "surfaces  = " << model_data.surfaces.size()                 << '\n';
    ostream << "lines     = " << model_data.lines.size()                    << '\n';
    return ostream;
}

} // namespace fem::model
