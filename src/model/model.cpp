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

/**
 * Initializes persistent element state required during one analysis step.
 *
 * Every compiled structural element receives a lifecycle callback. Elements
 * that do not require persistent step data use the default no-op implementation,
 * while formulations such as finite-rotation shells construct their cached
 * reference geometry here.
 *
 * Initialization is transactional at model level: if one element throws, the
 * method invokes `step_end()` for the complete model to release every cache that
 * may already have been created, then propagates the original exception.
 */
void Model::step_begin() {
    // Initialize analysis-local state for every compiled structural element
    try {
        for (auto& elem : _data->elements) {
            if (elem) {
                if (auto* structural = elem->as<StructuralElement>()) {
                    structural->step_begin();
                }
            }
        }
    } catch (...) {
        // Release partially initialized state before propagating the failure
        step_end();
        throw;
    }
}

/**
 * Releases persistent element state associated with the current analysis step.
 *
 * The callback is sent to every compiled structural element and is safe after a
 * partially completed `step_begin()`. Persistent model topology, material
 * history and reusable thread-local evaluation workspaces remain unchanged.
 */
void Model::step_end() {
    // Release only the analysis-local caches owned by structural elements
    for (auto& elem : _data->elements) {
        if (elem) {
            if (auto* structural = elem->as<StructuralElement>()) {
                structural->step_end();
            }
        }
    }
}

/**
 * Transfers one load into the currently active load collector.
 *
 * Loads reference compiled assembly regions, so registration is permitted only
 * after the semantic topology has been flattened. Ownership of the polymorphic
 * load is moved into the active collector.
 *
 * @param load Load definition to register.
 */
void Model::add_load(bc::Load::Ptr load) {
    // Validate the compiled model state, load ownership and collector target
    logging::error(_data != nullptr && _data->compiled,
        "Model: loads require a compiled model");
    logging::error(load != nullptr,
        "Model: cannot add a null load");
    logging::error(_data->load_cols.has_any() && _data->load_cols.get() != nullptr,
        "Model: no load collector is active");

    // Transfer the validated load into the selected collector
    _data->load_cols.get()->add(std::move(load));
}

/**
 * Registers one named amplitude as a shared model definition.
 *
 * Amplitudes are independent of part/instance identifier flattening and may be
 * added before or after compilation. Names are unique within the model-wide
 * amplitude dictionary, which takes shared ownership of the supplied object.
 *
 * @param amplitude Named time-history definition to register.
 */
void Model::add_amplitude(bc::Amplitude::Ptr amplitude) {
    // Validate ownership and the model-wide amplitude namespace
    logging::error(amplitude != nullptr,
        "Model: cannot add a null amplitude");
    logging::error(!_data->amplitudes.has(amplitude->name),
        "Model: amplitude ", amplitude->name, " is already defined");

    // Transfer the unique named definition into persistent model storage
    _data->amplitudes.add(std::move(amplitude));
}

/**
 * Transfers one support into the currently active support collector.
 *
 * Supports resolve compiled assembly regions when constraint equations are
 * collected. Registration therefore requires a compiled model and an active
 * collector. The support value is moved into that collector.
 *
 * @param support Support definition to register.
 */
void Model::add_support(bc::Support support) {
    // Validate the compiled model state and collector target
    logging::error(_data != nullptr && _data->compiled,
        "Model: supports require a compiled model");
    logging::error(_data->supp_cols.has_any() && _data->supp_cols.get() != nullptr,
        "Model: no support collector is active");

    // Transfer the support into the selected collector
    _data->supp_cols.get()->add(std::move(support));
}

/**
 * Adds a concentrated mass, inertia and spring contribution to a node region.
 *
 * The target node set must already exist in compiled assembly space. A
 * `PointMass` feature stores the resolved region together with translational
 * mass, rotary inertia and the translational and rotational spring coefficients.
 * Later stiffness and mass assembly operations evaluate the registered feature
 * on the active global DOFs.
 *
 * @param nset Name of the compiled target node region.
 * @param mass Concentrated translational mass.
 * @param rotary_inertia Rotary inertia about the three global axes.
 * @param spring Translational spring coefficients along the global axes.
 * @param rotary_spring Rotational spring coefficients about the global axes.
 */
void Model::add_point_mass_feature(const std::string& nset,
                                   Precision mass,
                                   Vec3 rotary_inertia,
                                   Vec3 spring,
                                   Vec3 rotary_spring) {
    // Validate the compiled assembly and resolve the target node region
    logging::error(_data != nullptr && _data->compiled,
        "Model: point-mass features require a compiled model");
    logging::error(_data->node_sets.has(nset),
        "Node set ", nset, " is not a defined node set");

    // Construct the non-element contribution and retain it for global assembly
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
    ostream << "parts     = " << entry_count(model_data.parts)            << '\n';
    ostream << "instances = " << entry_count(model_data.instances)        << '\n';
    ostream << "nodes     = " << node_count                               << '\n';
    ostream << "elements  = " << model_data.elements.size()               << '\n';
    ostream << "surfaces  = " << model_data.surfaces.size()               << '\n';
    ostream << "lines     = " << model_data.lines.size()                  << '\n';
    return ostream;
}

} // namespace fem::model
