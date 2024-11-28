#include "load_collector.h"

namespace fem {

LoadCollector::LoadCollector(const std::string& p_name)
    : model::Collection<Load::Ptr>(p_name) {}

void LoadCollector::apply(model::ModelData& model_data, NodeData& bc) {
    for (const auto& load : this->_data) {
        load->apply(model_data, bc);
    }
}

void LoadCollector::add_cload(model::NodeRegion::Ptr region, Vec6 values) {
    auto cload = std::make_shared<CLoad>();
    cload->region = region;
    cload->values = values;
    this->add(cload);
}

void LoadCollector::add_dload(model::SurfaceRegion::Ptr region, Vec3 values) {
    auto dload = std::make_shared<DLoad>();
    dload->region = region;
    dload->values = values;
    this->add(dload);
}

void LoadCollector::add_pload(model::SurfaceRegion::Ptr region, Precision pressure) {
    auto pload = std::make_shared<PLoad>();
    pload->region = region;
    pload->pressure = pressure;
    this->add(pload);
}

void LoadCollector::add_vload(model::ElementRegion::Ptr region, Vec3 values) {
    auto vload = std::make_shared<VLoad>();
    vload->region = region;
    vload->values = values;
    this->add(vload);
}

} // namespace fem
