
#include "model.h"
namespace fem {

namespace model {

void Model::activate_node_set(const std::string &name) {
    node_sets.activate(name);
}
void Model::activate_element_set(const std::string &name) {
    elem_sets.activate(name);
}
void Model::activate_load_set(const std::string &name) {
    load_sets.activate(name, max_nodes, 6);
    load_sets.current().setZero();
}
void Model::activate_support_set(const std::string &name) {
    support_sets.activate(name, max_nodes, 6);
    support_sets.current().setZero();
    support_sets.current().fill(std::numeric_limits<Precision>::quiet_NaN());
}
void Model::activate_material(const std::string &name){
    materials.activate(name);
}

void Model::add_cload(const std::string& nset, StaticVector<3> load){
    for(ID id:node_sets.get(nset)){
        add_cload(id, load);
    }
}
void Model::add_cload(const ID id, StaticVector<3> load){
    for(int i = 0; i < 3; i++){
        load_sets.current()(id,i) += load(i);
        if(!load_sets.is_default_set())
            load_sets.all()(id,i) += load(i);
    }
}

void Model::add_vload(const std::string& elset, StaticVector<3> load){
    for(ID id:elem_sets.get(elset)){
        add_vload(id, load);
    }
}
void Model::add_vload(const ID id, StaticVector<3> load){
    for(int i = 0; i < 3; i++){
        elements[i]->apply_vload(node_coords, load_sets.current(), load);
        if(!load_sets.is_default_set())
            elements[i]->apply_vload(node_coords, load_sets.all(), load);
    }
}

void Model::add_support(const std::string& nset, const StaticVector<6> constraint){
    for(ID id:node_sets.get(nset)){
        add_support(id, constraint);
    }
}
void Model::add_support(const std::string& nset, const StaticVector<3> displacement){
    for(ID id:node_sets.get(nset)){
        add_support(id, displacement);
    }
}
void Model::add_support(const ID id, const StaticVector<6> constraint){
    for(int i = 0; i < 6; i++){
        support_sets.current()(id,i) = constraint(i);
        if(!support_sets.is_default_set())
            support_sets.all()(id,i) = constraint(i);
    }
}
void Model::add_support(const ID id, const StaticVector<3> displacement){
    for(int i = 0; i < 3; i++){
        support_sets.current()(id,i) = displacement(i);
        if(!support_sets.is_default_set())
            support_sets.all()(id,i) = displacement(i);
    }
}
void Model::add_support_rot(const std::string& nset, const StaticVector<3> rotation){
    for(ID id:node_sets.get(nset)){
        add_support_rot(id, rotation);
    }
}
void Model::add_support_rot(const ID id, const StaticVector<3> rotation){
    for(int i = 3; i < 6; i++){
        support_sets.current()(id,i) = rotation(i);
        if(!support_sets.is_default_set())
            support_sets.all()(id,i) = rotation(i);
    }
}
void Model::add_support(const ID id, const Dim dim, const Precision displacement){
    support_sets.current()(id,dim) = displacement;
}

material::Material& Model::active_material(){
    return materials.current();
}
NodeData& Model::active_loads(){
    return load_sets.current();
}
NodeData& Model::active_supports(){
    return support_sets.current();
}
std::vector<ID>& Model::active_nodeset(){
    return node_sets.current();
}
std::vector<ID>& Model::active_elemset(){
    return elem_sets.current();
}

Sets<std::vector<ID>>&  Model::nodesets(){
    return node_sets;
}

Sets<std::vector<ID>>&  Model::elemsets(){
    return elem_sets;
}

void Model::solid_section(const std::string& set, const std::string& material){
    material::Material* mat_ptr = &materials.get(material);
    for(ID id:elem_sets.get(set)){
        elements[id]->set_material(mat_ptr);
    }
}

std::ostream& operator<<(std::ostream& ostream, const model::Model& model) {
    ostream << "Model (dim = " << model.element_dims << ")\n";
    ostream << "\tmax_nodes = " << model.max_nodes << '\n';
    ostream << "\tmax_elements = " << model.max_elements << '\n';

    ostream << "\tNode sets:\n";
    for (const auto& set : model.node_sets.m_sets) {
        ostream << "\t\t" << set.first << ": " << set.second.size() << '\n';
    }

    ostream << "\tElement sets:\n";
    for (const auto& set : model.elem_sets.m_sets) {
        ostream << "\t\t" << set.first << ": " << set.second.size() << '\n';
    }

    ostream << "\tLoad sets:\n";
    for (const auto& set : model.load_sets.m_sets) {
        ostream << "\t\t" << set.first << '\n';
        ostream << "\t\t\t" << "count=" << set.second.count() << '\n';
    }

    ostream << "\tSupport sets:\n";
    for (const auto& set : model.support_sets.m_sets) {
        ostream << "\t\t" << set.first << '\n';
        ostream << "\t\t\t" << "count=" << set.second.count() << '\n';
    }

    ostream << "\tMaterial sets:\n";
    for (const auto& set : model.materials.m_sets) {
        ostream << "\t\t" << set.first << '\n';
        if(set.second.has_density()) {
            ostream << "\t\t\tdensity=" << set.second.density() << '\n';
        }
        if(set.second.has_elasticity()) {
            ostream << "\t\t\telastic" << '\n';
        }
    }

    return ostream;
}


}    // namespace model
}    // namespace fem
