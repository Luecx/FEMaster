
#include "model.h"
namespace fem {

namespace model {

void Model::add_connector(const std::string& set1,
                          const std::string& set2,
                          const std::string& coordinate_system,
                          constraint::ConnectorType type) {
    logging::error(node_sets.has(set1), "Node set ", set1, " does not exist");
    logging::error(node_sets.has(set2), "Node set ", set2, " does not exist");

    logging::error(coordinate_systems.has(coordinate_system), "Coordinate system ", coordinate_system, " does not exist");

    // can only contain exactly one node
    logging::error(node_sets.get(set1)->size() == 1, "Set 1 must contain exactly one node");
    logging::error(node_sets.get(set2)->size() == 1, "Set 2 must contain exactly one node");

    ID id1 = node_sets.get(set1)->first();
    ID id2 = node_sets.get(set2)->first();

    _connectors.push_back(constraint::Connector{id1, id2, coordinate_systems.get(coordinate_system), type});
}

void Model::add_coupling(const std::string &master_set, const std::string &slave_set, Dofs coupled_dofs, constraint::CouplingType type) {
    logging::error(node_sets.get(master_set)->size() == 1, "Master set must contain exactly one node");
    ID master_node = node_sets.get(master_set)->first();

    _couplings.push_back({master_node, node_sets.get(slave_set), coupled_dofs, type});
}

void Model::add_tie(const std::string& master_set, const std::string& slave_set, Precision distance, bool adjust) {
    // check that the sets exist for the surfaces
    logging::error(surface_sets.has(master_set), "Master set ", master_set, " is not a defined surface set");
    logging::error(node_sets   .has(slave_set) , "Slave set " , slave_set , " is not a defined node set");

    NodeRegion::Ptr    slave_ptr  = node_sets.get(slave_set);
    SurfaceRegion::Ptr master_ptr = surface_sets.get(master_set);
    _ties.push_back(constraint::Tie{master_ptr, slave_ptr, distance, adjust});
}


void Model::add_cload(const std::string& nset, Vec6 load){
    for(ID id:*node_sets.get(nset)){
        add_cload(id, load);
    }
}
void Model::add_cload(const ID id, Vec6 load) {
    for(int i = 0; i < 6; i++){
        (*(_load_sets.get()))(id,i) += load(i);
    }
}

void Model::add_dload(const std::string& sfset, Vec3 load) {
    for(ID id:*surface_sets.get(sfset)){
        add_dload(id, load);
    }
}
void Model::add_dload(ID id, Vec3 load){
    surfaces[id]->apply_dload(node_coords, (*_load_sets.get()), load);
}

void Model::add_vload(const std::string& elset, Vec3 load){
    for(ID id:*elem_sets.get(elset)){
        add_vload(id, load);
    }
}
void Model::add_vload(const ID id, Vec3 load){
    if (elements[id] == nullptr) return;
    if (!elements[id] ->is_type(StructuralType))
        return;
    auto sel = elements[id]->as<StructuralElement>();

    sel->apply_vload(node_coords, (*_load_sets.get()), load);
}

void Model::add_tload(std::string& temp_field, Precision ref_temp) {
    logging::error(_fields_temperature.has(temp_field), "Temperature field ", temp_field, " does not exist");

    auto temp_ptr = _fields_temperature.get(temp_field);

    for(ElementPtr& elem:elements){
        if (elem == nullptr)
            continue;
        if(!elem->is_type(StructuralType))
            continue;
        auto sel = elem->as<StructuralElement>();

        sel->apply_tload(node_coords, (*_load_sets.get()), *temp_ptr, ref_temp);
    }
}

void Model::add_support(const std::string& nset, const StaticVector<6> constraint){
    for(ID id:*node_sets.get(nset)){
        add_support(id, constraint);
    }
}
void Model::add_support(const ID id, const StaticVector<6> constraint){
    for(int i = 0; i < 6; i++){
        (*_support_sets.get())(id,i) = constraint(i);
    }
}

void Model::set_field_temperature(const std::string& name, const ID id, Precision value) {
    if (!_fields_temperature.has(name)) {
        _fields_temperature.activate(name, max_nodes, 1);
        _fields_temperature.get()->fill(std::numeric_limits<Precision>::quiet_NaN());
    }
    _fields_temperature.activate(name, max_nodes, 1);
    _fields_temperature.get()->operator ()(id) = value;
}

void Model::solid_section(const std::string& set, const std::string& material){
    material::Material::Ptr ptr = _materials.get(material);
    for(ID id:*elem_sets.get(set)){
        elements[id]->set_material(ptr);
    }
}

std::ostream& operator<<(std::ostream& ostream, const model::Model& model) {
    ostream << "Model (dim = " << model.element_dims << ")\n";
    ostream << "\tmax nodes = " << model.max_nodes << '\n';
    ostream << "\tmax elements = " << model.max_elements << '\n';
    ostream << "\tmax surfaces = " << model.max_surfaces << '\n';

//    ostream << "\tNode sets:\n";
//    for (const auto& set : model.node_sets.m_sets) {
//        ostream << "\t\t" << set.first << ": " << set.second.size() << '\n';
//    }

//    ostream << "\tElement sets:\n";
//    for (const auto& set : model.elem_sets.m_sets) {
//        ostream << "\t\t" << set.first << ": " << set.second.size() << '\n';
//    }
//
//    ostream << "\tSurface sets:\n";
//    for (const auto& set : model.surface_sets.m_sets) {
//        ostream << "\t\t" << set.first << ": " << set.second.size() << '\n';
//    }
//
//    ostream << "\tLoad sets:\n";
//    for (const auto& set : model._load_sets.m_sets) {
//        ostream << "\t\t" << set.first << '\n';
//        ostream << "\t\t\t" << "count=" << set.second.count() << '\n';
//    }
//
//    ostream << "\tSupport sets:\n";
//    for (const auto& set : model._support_sets.m_sets) {
//        ostream << "\t\t" << set.first << '\n';
//        ostream << "\t\t\t" << "count=" << set.second.count() << '\n';
//    }

//    ostream << "\tMaterial sets:\n";
//    for (const auto& set : model._materials.m_sets) {
//        ostream << "\t\t" << set.first << '\n';
//        if(set.second.has_density()) {
//            ostream << "\t\t\tdensity=" << set.second.get_density() << '\n';
//        }
//        if(set.second.has_elasticity()) {
//            ostream << "\t\t\telastic" << '\n';
//        }
//    }

    return ostream;
}


}    // namespace model
}    // namespace fem
