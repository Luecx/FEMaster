/**
 * @file load_collector.cpp
 * @brief Implements aggregation of load boundary conditions.
 *
 * The `LoadCollector` centralizes creation and application of heterogeneous
 * loads, ensuring that all registered loads are assembled consistently during
 * load-case execution.
 *
 * @see src/bc/load_collector.h
 * @see src/bc/load_c.cpp
 * @see src/bc/load_d.cpp
 * @see src/bc/load_p.cpp
 * @see src/bc/load_v.cpp
 * @see src/bc/load_t.cpp
 * @see src/bc/load_inertial.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_collector.h"

#include "load_c.h"
#include "load_d.h"
#include "load_inertial.h"
#include "load_p.h"
#include "load_t.h"
#include "load_v.h"

#include "../data/field.h"

#include <memory>
#include <utility>

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
    cload->region_      = std::move(region);
    cload->values_      = values;
    cload->orientation_ = std::move(orientation);
    cload->amplitude_   = std::move(amplitude);
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
    dload->region_      = std::move(region);
    dload->values_      = values;
    dload->orientation_ = std::move(orientation);
    dload->amplitude_   = std::move(amplitude);
    this->add(dload);
}

/**
 * @copydoc LoadCollector::add_pload
 */
void LoadCollector::add_pload(model::SurfaceRegion::Ptr region,
                              Precision pressure,
                              Amplitude::Ptr amplitude) {
    auto pload = std::make_shared<PLoad>();
    pload->region_    = std::move(region);
    pload->pressure_  = pressure;
    pload->amplitude_ = std::move(amplitude);
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
    vload->region_      = std::move(region);
    vload->values_      = values;
    vload->orientation_ = std::move(orientation);
    vload->amplitude_   = std::move(amplitude);
    this->add(vload);
}

/**
 * @copydoc LoadCollector::add_tload
 */
void LoadCollector::add_tload(model::Field::Ptr temp_field, Precision ref_temp) {
    auto tload = std::make_shared<TLoad>();
    tload->temp_field_ = std::move(temp_field);
    tload->ref_temp_   = ref_temp;
    this->add(tload);
}

/**
 * @copydoc LoadCollector::add_inertialload
 */
void LoadCollector::add_inertialload(model::ElementRegion::Ptr region,
                                     Vec3 center,
                                     Vec3 center_acc,
                                     Vec3 omega,
                                     Vec3 alpha,
                                     bool consider_point_masses) {
    auto inertial_load = std::make_shared<InertialLoad>();
    inertial_load->region_                = std::move(region);
    inertial_load->center_                = center;
    inertial_load->center_acc_            = center_acc;
    inertial_load->omega_                 = omega;
    inertial_load->alpha_                 = alpha;
    inertial_load->consider_point_masses_ = consider_point_masses;
    this->add(inertial_load);
}
} // namespace bc
} // namespace fem
