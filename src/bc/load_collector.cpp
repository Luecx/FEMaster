/**
 * @file load_collector.cpp
 * @brief Implements aggregation of load boundary conditions.
 *
 * The `LoadCollector` centralizes creation and application of heterogeneous
 * loads, ensuring that all registered loads are assembled consistently during
 * load-case execution.
 *
 * @see src/bc/load_collector.h
 * @see src/bc/load.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_collector.h"

#include "../data/field.h"

namespace fem {
namespace bc {

/**
 * @copydoc LoadCollector::LoadCollector
 */
LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Load::Ptr>(name) {}

/**
 * @copydoc LoadCollector::apply
 */
void LoadCollector::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    for (const auto& load : this->_data) {
        if (!load) {
            continue;
        }
        load->apply(model_data, bc, time);
    }
}

/**
 * @copydoc LoadCollector::add_cload
 */
void LoadCollector::add_cload(model::NodeRegion::Ptr region,
                              Vec6 values,
                              cos::CoordinateSystem::Ptr orientation,
                              Amplitude::Ptr amplitude) {
    auto cload = std::make_shared<CLoad>();
    cload->region = std::move(region);
    cload->values = values;
    cload->orientation = std::move(orientation);
    cload->amplitude = std::move(amplitude);
    this->add(cload);
}

/**
 * @copydoc LoadCollector::add_dload
 */
void LoadCollector::add_dload(model::SurfaceRegion::Ptr region,
                              Vec3 values,
                              cos::CoordinateSystem::Ptr orientation,
                              Amplitude::Ptr amplitude) {
    auto dload = std::make_shared<DLoad>();
    dload->region = std::move(region);
    dload->values = values;
    dload->orientation = std::move(orientation);
    dload->amplitude = std::move(amplitude);
    this->add(dload);
}

/**
 * @copydoc LoadCollector::add_pload
 */
void LoadCollector::add_pload(model::SurfaceRegion::Ptr region,
                              Precision pressure,
                              Amplitude::Ptr amplitude) {
    auto pload = std::make_shared<PLoad>();
    pload->region = std::move(region);
    pload->pressure = pressure;
    pload->amplitude = std::move(amplitude);
    this->add(pload);
}

/**
 * @copydoc LoadCollector::add_vload
 */
void LoadCollector::add_vload(model::ElementRegion::Ptr region,
                              Vec3 values,
                              cos::CoordinateSystem::Ptr orientation,
                              Amplitude::Ptr amplitude) {
    auto vload = std::make_shared<VLoad>();
    vload->region = std::move(region);
    vload->values = values;
    vload->orientation = std::move(orientation);
    vload->amplitude = std::move(amplitude);
    this->add(vload);
}

/**
 * @copydoc LoadCollector::add_tload
 */
void LoadCollector::add_tload(model::Field::Ptr temp_field, Precision ref_temp) {
    auto tload = std::make_shared<TLoad>();
    tload->temp_field = std::move(temp_field);
    tload->ref_temp = ref_temp;
    this->add(tload);
}

void LoadCollector::add_inertialload(model::ElementRegion::Ptr region,
                                     Vec3 center,
                                     Vec3 center_acc,
                                     Vec3 omega,
                                     Vec3 alpha) {
    auto il = std::make_shared<InertialLoad>();
    il->region      = std::move(region);
    il->center      = center;
    il->center_acc  = center_acc;
    il->omega       = omega;
    il->alpha       = alpha;
    this->add(il);
}

} // namespace bc
} // namespace fem
