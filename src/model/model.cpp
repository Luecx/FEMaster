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
    for (ID id : *_data->node_sets.get(nset)) {
        add_cload(id, load);
    }
}

void Model::add_cload(const ID id, Vec6 load) {
    for (int i = 0; i < 6; i++) {
        (*(_load_sets.get()))(id, i) += load(i);
    }
}

void Model::add_dload(const std::string& sfset, Vec3 load) {
    for (ID id : *_data->surface_sets.get(sfset)) {
        add_dload(id, load);
    }
}

void Model::add_dload(ID id, Vec3 load) {
    _data->surfaces[id]->apply_dload(_data->get(NodeDataEntries::POSITION), (*_load_sets.get()), load);
}

void Model::add_vload(const std::string& elset, Vec3 load) {
    for (ID id : *_data->elem_sets.get(elset)) {
        add_vload(id, load);
    }
}

void Model::add_vload(const ID id, Vec3 load) {
    if (_data->elements[id] == nullptr) return;
    if (auto sel = _data->elements[id]->as<StructuralElement>())
        sel->apply_vload((*_load_sets.get()), load);
}

void Model::add_tload(std::string& temp_field, Precision ref_temp) {
    // TODO
    logging::error(_fields_temperature.has(temp_field), "Temperature field ", temp_field, " does not exist");

    auto temp_ptr = _fields_temperature.get(temp_field);
    for (ElementPtr& elem : _data->elements) {
        if (elem == nullptr) continue;
        if (auto sel = elem->as<StructuralElement>())
            sel->apply_tload((*_load_sets.get()), *temp_ptr, ref_temp);
    }
}

void Model::add_support(const std::string& nset, const StaticVector<6> constraint) {
    for (ID id : *_data->node_sets.get(nset)) {
        add_support(id, constraint);
    }
}

void Model::add_support(const ID id, const StaticVector<6> constraint) {
    for (int i = 0; i < 6; i++) {
        (*_support_sets.get())(id, i) = constraint(i);
    }
}

void Model::set_field_temperature(const std::string& name, ID id, Precision value) {
    if (!_fields_temperature.has(name)) {
        _fields_temperature.activate(name, _data->max_nodes, 1);
        _fields_temperature.get()->fill(std::numeric_limits<Precision>::quiet_NaN());
    }
    _fields_temperature.activate(name, _data->max_nodes, 1);
    _fields_temperature.get()->operator()(id) = value;
}

void Model::solid_section(const std::string& set, const std::string& material) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material), "Material ", material, " is not a defined material");

    Section::Ptr section = std::make_shared<Section>();
    section->material = _data->materials.get(material);
    section->region   = _data->elem_sets.get(set);
    this->_data->sections.push_back(section);
}

std::ostream& operator<<(std::ostream& ostream, const model::Model& model) {
    ostream << "\tmax nodes = " << model._data->max_nodes << '\n';
    ostream << "\tmax elements = " << model._data->max_elems << '\n';
    ostream << "\tmax surfaces = " << model._data->max_surfaces << "\n";

    // print materials
    ostream << "\tMaterials:\n";
    for (const auto& [name, mat] : model._data->materials._data) {
        ostream << "\t\t" << name << ":\n";
    }

    return ostream;
}

} // namespace model
} // namespace fem
