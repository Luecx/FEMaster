#include "model.h"

namespace fem {
namespace model {

void Model::add_connector(const std::string& set1,
                          const std::string& set2,
                          const std::string& coordinate_system,
                          constraint::ConnectorType type) {

    logging::error(_data->node_sets.has(set1), "Node set ", set1, " does not exist");
    logging::error(_data->node_sets.has(set2), "Node set ", set2, " does not exist");

    logging::error(_data->coordinate_systems.has(coordinate_system), "Coordinate system ", coordinate_system, " does not exist");

    logging::error(_data->node_sets.get(set1)->size() == 1, "Set 1 must contain exactly one node");
    logging::error(_data->node_sets.get(set2)->size() == 1, "Set 2 must contain exactly one node");

    ID id1 = _data->node_sets.get(set1)->first();
    ID id2 = _data->node_sets.get(set2)->first();

    _data->connectors.emplace_back(id1, id2, _data->coordinate_systems.get(coordinate_system), type);
}

void Model::add_coupling(const std::string &master_set, const std::string &slave_set, Dofs coupled_dofs, constraint::CouplingType type, bool is_surface) {
    logging::error(_data->node_sets.get(master_set)->size() == 1, "Master set must contain exactly one node");

    if (is_surface) {
        logging::error(_data->surface_sets.has(slave_set), "Slave set ", slave_set, " is not a defined surface set");
    } else {
        logging::error(_data->node_sets.has(slave_set), "Slave set ", slave_set, " is not a defined node set");
    }

    ID master_node = _data->node_sets.get(master_set)->first();

    if (is_surface) {
        _data->couplings.emplace_back(master_node, _data->surface_sets.get(slave_set), coupled_dofs, type);
    } else {
        _data->couplings.emplace_back(master_node, _data->node_sets.get(slave_set), coupled_dofs, type);
    }
}

void Model::add_tie(const std::string& master_set, const std::string& slave_set, Precision distance, bool adjust) {
    logging::error(_data->surface_sets.has(master_set), "Master set ", master_set, " is not a defined surface set");
    logging::error(_data->node_sets.has(slave_set), "Slave set ", slave_set, " is not a defined node set");

    NodeRegion::Ptr slave_ptr = _data->node_sets.get(slave_set);
    SurfaceRegion::Ptr master_ptr = _data->surface_sets.get(master_set);
    _data->ties.emplace_back(master_ptr, slave_ptr, distance, adjust);
}

void Model::add_cload(const std::string& nset, Vec6 load){
    logging::error(_data->node_sets.has(nset), "Node set ", nset, " does not exist");
    auto region_ptr = _data->node_sets.get(nset);
    _data->load_cols.get()->add_cload(region_ptr, load);
}

void Model::add_cload(const ID id, Vec6 load) {
    auto region_ptr = std::make_shared<NodeRegion>("INTERNAL");
    region_ptr->add(id);
    _data->load_cols.get()->add_cload(region_ptr, load);
}

void Model::add_dload(const std::string& sfset, Vec3 load) {
    logging::error(_data->surface_sets.has(sfset), "Surface set ", sfset, " does not exist");
    auto region_ptr = _data->surface_sets.get(sfset);
    _data->load_cols.get()->add_dload(region_ptr, load);
}

void Model::add_dload(ID id, Vec3 load) {
    auto region_ptr = std::make_shared<SurfaceRegion>("INTERNAL");
    region_ptr->add(id);
    _data->load_cols.get()->add_dload(region_ptr, load);
}

void Model::add_pload(const std::string& sfset, Precision load) {
    logging::error(_data->surface_sets.has(sfset), "Surface set ", sfset, " does not exist");
    auto region_ptr = _data->surface_sets.get(sfset);
    _data->load_cols.get()->add_pload(region_ptr, load);
}

void Model::add_pload(ID id, Precision load) {
    auto region_ptr = std::make_shared<SurfaceRegion>("INTERNAL");
    region_ptr->add(id);
    _data->load_cols.get()->add_pload(region_ptr, load);
}

void Model::add_vload(const std::string& elset, Vec3 load) {
    logging::error(_data->elem_sets.has(elset), "Element set ", elset, " does not exist");
    auto region_ptr = _data->elem_sets.get(elset);
    _data->load_cols.get()->add_vload(region_ptr, load);
}

void Model::add_vload(const ID id, Vec3 load) {
    auto region_ptr = std::make_shared<ElementRegion>("INTERNAL");
    region_ptr->add(id);
    _data->load_cols.get()->add_vload(region_ptr, load);
}

void Model::add_tload(std::string& temp_field, Precision ref_temp) {
    if (!_data->node_fields.has(TEMPERATURE)) {
        _data->node_fields.activate(TEMPERATURE);
    }
    auto temp_fields = _data->node_fields.get(TEMPERATURE);

    logging::error(temp_fields->has(temp_field), "Temperature field ", temp_field, " does not exist");
    _data->load_cols.get()->add_tload(temp_fields->get(temp_field), ref_temp);

    // // TODO
    // logging::error(_fields_temperature.has(temp_field), "Temperature field ", temp_field, " does not exist");
    //
    // auto temp_ptr = _fields_temperature.get(temp_field);
    // for (ElementPtr& elem : _data->elements) {
    //     if (elem == nullptr) continue;
    //     if (auto sel = elem->as<StructuralElement>())
    //         sel->apply_tload((*_data->load_cols.get()), *temp_ptr, ref_temp);
    // }
}

void Model::add_support(const std::string& nset, const StaticVector<6> constraint, const std::string& orientation) {
    logging::error(_data->supp_cols.has_any(), "No support collectors have been defined");
    logging::error(_data->supp_cols.get() != nullptr, "No support collectors is currently active");
    if (!orientation.empty())
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");

    SupportCollector::Ptr supp_col = _data->supp_cols.get();
    supp_col->add_supp(_data->node_sets.get(nset), constraint, _data->coordinate_systems.get(orientation));
}

void Model::add_support(const ID id, const StaticVector<6> constraint, const std::string& orientation) {
    logging::error(_data->supp_cols.has_any(), "No support collectors have been defined");
    logging::error(_data->supp_cols.get() != nullptr, "No support collectors is currently active");
    if (!orientation.empty())
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");

    // create a new NodeRegion
    NodeRegion::Ptr region = std::make_shared<NodeRegion>("INTERNAL");
    region->add(id);

    SupportCollector::Ptr supp_col = _data->supp_cols.get();
    supp_col->add_supp(region, constraint, _data->coordinate_systems.get(orientation));
}

void Model::set_field_temperature(const std::string& name, ID id, Precision value) {
    if (!_data->node_fields.has(TEMPERATURE)) {
        _data->node_fields.activate(TEMPERATURE);
    }

    auto temp_fields = _data->node_fields.get(TEMPERATURE);
    if (temp_fields->has(name)) {
        temp_fields->get(name)->operator()(id) = value;
    } else {
        temp_fields->activate(name, _data->max_nodes, 1);
        temp_fields->get(name)->fill(std::numeric_limits<Precision>::quiet_NaN());
        temp_fields->get(name)->operator()(id) = value;
    }
}

void Model::solid_section(const std::string& set, const std::string& material) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material), "Material ", material, " is not a defined material");

    Section::Ptr section = std::make_shared<SolidSection>();
    section->material = _data->materials.get(material);
    section->region   = _data->elem_sets.get(set);
    this->_data->sections.push_back(section);
}

void Model::beam_section(const std::string& set, const std::string& material,  const std::string& profile, Vec3 orientation) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material), "Material ", material, " is not a defined material");
    logging::error(_data->profiles.has(profile), "Profile ", profile, " is not a defined profile");
    BeamSection::Ptr sec = std::make_shared<BeamSection>();
    sec->material = _data->materials.get(material);
    sec->region   = _data->elem_sets.get(set);
    sec->profile  = _data->profiles.get(profile);
    sec->n1       = orientation;
    this->_data->sections.push_back(sec);
}

void Model::shell_section(const std::string& set, const std::string& material, Precision thickness) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material), "Material ", material, " is not a defined material");
    ShellSection::Ptr sec = std::make_shared<ShellSection>();
    sec->material = _data->materials.get(material);
    sec->region = _data->elem_sets.get(set);
    sec->thickness = thickness;
    this->_data->sections.push_back(sec);
}

void Model::point_mass_section(const std::string& set, Precision mass, Vec3 rotary_inertia, Vec3 spring, Vec3 rotary_spring) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    PointMassSection::Ptr sec = std::make_shared<PointMassSection>();
    sec->mass = mass;
    sec->rotary_inertia = rotary_inertia;
    sec->spring_constants = spring;
    sec->rotary_spring_constants = rotary_spring;
    sec->region = _data->elem_sets.get(set);
    this->_data->sections.push_back(sec);
}

std::ostream& operator<<(std::ostream& ostream, const model::Model& model) {
    ostream << "max nodes = " << model._data->max_nodes << '\n';
    ostream << "max elements = " << model._data->max_elems << '\n';
    ostream << "max surfaces = " << model._data->max_surfaces << "\n";

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
    for(const auto &profile : model._data->profiles){
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

    // // print materials
    // ostream << "Materials:\n";
    // for (const auto& material : model._data->materials) {
    //     ostream << *material.second;
    // }
    // return ostream;
}

} // namespace model
} // namespace fem
