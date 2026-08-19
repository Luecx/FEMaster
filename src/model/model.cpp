/**
 * @file model.cpp
 * @brief Implements compact model-level object registration and lifecycle logic.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "model.h"

#include "../bc/load_collector.h"
#include "../bc/support_collector.h"
#include "../feature/point_mass.h"

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

std::ostream& operator<<(std::ostream& ostream, const Model& model) {
    const Index nodes = model._data->positions ? model._data->positions->rows : Index(0);
    ostream << "nodes = "    << nodes                         << '\n';
    ostream << "elements = " << model._data->elements.size() << '\n';
    ostream << "surfaces = " << model._data->surfaces.size() << '\n';

    logging::info(true, "Materials");
    logging::up();
    for (const auto& material : model._data->materials) {
        material.second->info();
    }
    logging::down();

    logging::info(true, "Sections");
    logging::up();
    for (const auto& section : model._data->sections) {
        section->info();
    }
    logging::down();

    logging::info(true, "Profiles");
    logging::up();
    for (const auto& profile : model._data->profiles) {
        profile.second->info();
    }
    logging::down();

    logging::info(true, "Element sets");
    logging::up();
    for (const auto& elem_set : model._data->elem_sets) {
        elem_set.second->info();
    }
    logging::down();

    logging::info(true, "Node sets");
    logging::up();
    for (const auto& node_set : model._data->node_sets) {
        node_set.second->info();
    }
    logging::down();

    return ostream;
}

} // namespace fem::model
