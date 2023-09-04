
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
    load_sets.activate(name);
}
void Model::activate_support_set(const std::string &name) {
    support_sets.activate(name);
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
    load_sets.current()(id,0) = load(0);
    load_sets.current()(id,1) = load(1);
    load_sets.current()(id,2) = load(2);
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
    support_sets.current()(id,0) = constraint(0);
    support_sets.current()(id,1) = constraint(1);
    support_sets.current()(id,2) = constraint(2);
    support_sets.current()(id,3) = constraint(3);
    support_sets.current()(id,4) = constraint(4);
    support_sets.current()(id,5) = constraint(5);
}
void Model::add_support(const ID id, const StaticVector<3> displacement){
    support_sets.current()(id,0) = displacement(0);
    support_sets.current()(id,1) = displacement(1);
    support_sets.current()(id,2) = displacement(2);
}
void Model::add_support_rot(const std::string& nset, const StaticVector<3> rotation){
    for(ID id:node_sets.get(nset)){
        add_support_rot(id, rotation);
    }
}
void Model::add_support_rot(const ID id, const StaticVector<3> rotation){
    support_sets.current()(id,3) = rotation(0);
    support_sets.current()(id,4) = rotation(1);
    support_sets.current()(id,5) = rotation(2);
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
std::vector<ID> Model::active_nodeset(){
    return node_sets.current();
}
std::vector<ID> Model::active_elemset(){
    return elem_sets.current();
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
    }

    ostream << "\tSupport sets:\n";
    for (const auto& set : model.support_sets.m_sets) {
        ostream << "\t\t" << set.first << '\n';
    }

    ostream << "\tMaterial sets:\n";
    for (const auto& set : model.materials.m_sets) {
        ostream << "\t\t" << set.first << '\n';
    }

    return ostream;
}


}    // namespace model
}    // namespace fem
