/**
 * @file load_collector.cpp
 * @brief Implements construction and collective application of load objects.
 *
 * The collector converts strongly typed parser or API parameters into concrete
 * load instances, transfers optional shared resources into them and stores the
 * resulting polymorphic pointers. During assembly it dispatches each valid load
 * in insertion order so all contributions accumulate in the same target field.
 *
 * @see load_collector.h
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
 * Constructs a named polymorphic collection for boundary-condition loads.
 *
 * Concrete concentrated, distributed, pressure, thermal, body and inertial
 * loads are stored through the common Load interface and applied together.
 *
 * @param name Name used to identify the load collection.
 */
LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Load::Ptr>(name) {}

/**
 * Applies every valid stored load in insertion order.
 *
 * All loads receive the same model data, target field and analysis time, so
 * their contributions accumulate in one generalized nodal result. Null entries
 * are skipped.
 *
 * @param model_data Model fields and topology required by the loads.
 * @param bc Generalized nodal field receiving all contributions.
 * @param time Analysis time forwarded to amplitude-aware loads.
 */
void LoadCollector::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    for (const auto& load : this->_data) {
        // Null pointers are tolerated so incomplete or conditionally generated
        // collections do not interrupt assembly of the remaining loads.
        if (!load) {
            continue;
        }

        // Dynamic dispatch selects the concrete region traversal and
        // integration strategy. Each load adds its contribution to `bc`.
        load->apply(model_data, bc, time);
    }
}

/**
 * Creates and stores a concentrated nodal force or moment load.
 *
 * Construction parameters are transferred to a concrete CLoad and appended to
 * the collection for later assembly.
 *
 * @param region Target node region.
 * @param values Generalized force and moment components.
 * @param orientation Optional local coordinate system.
 * @param amplitude Optional temporal multiplier.
 */
void LoadCollector::add_cload(model::NodeRegion::Ptr region,
                              Vec6 values,
                              cos::CoordinateSystem::Ptr orientation,
                              Amplitude::Ptr amplitude) {
    // Allocate the concrete type and transfer all construction parameters into
    // its public data members before exposing it through the base pointer.
    auto cload = std::make_shared<CLoad>();
    cload->region_      = std::move(region);
    cload->values_      = values;
    cload->orientation_ = std::move(orientation);
    cload->amplitude_   = std::move(amplitude);
    this->add(cload);
}

/**
 * Creates and stores a distributed surface-traction load.
 *
 * Geometric integration and nodal distribution are deferred until collection
 * application.
 *
 * @param region Target surface region.
 * @param values Traction vector.
 * @param orientation Optional local coordinate system.
 * @param amplitude Optional temporal multiplier.
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
 * Creates and stores a scalar pressure load on a surface region.
 *
 * Surface normals determine the global traction direction when the collection
 * is applied.
 *
 * @param region Target surface region.
 * @param pressure Scalar pressure value.
 * @param amplitude Optional temporal multiplier.
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
 * Creates and stores a distributed body-force load on an element region.
 *
 * Element-level integration and density scaling are deferred until application.
 *
 * @param region Target structural element region.
 * @param values Body-force vector.
 * @param orientation Optional local coordinate system.
 * @param amplitude Optional temporal multiplier.
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
 * Creates and stores a thermal load driven by a nodal temperature field.
 *
 * The field and stress-free reference temperature are retained for structural
 * element-specific thermal-force evaluation.
 *
 * @param temp_field Nodal temperature field.
 * @param ref_temp Stress-free reference temperature.
 */
void LoadCollector::add_tload(model::Field::Ptr temp_field, Precision ref_temp) {
    auto tload = std::make_shared<TLoad>();
    tload->temp_field_ = std::move(temp_field);
    tload->ref_temp_   = ref_temp;
    this->add(tload);
}

/**
 * Creates and stores a rigid-body inertial load.
 *
 * Reference-point accelerations and angular kinematics are retained until the
 * collection is applied; point-mass participation is optional.
 *
 * @param region Target structural element region.
 * @param center Reference point of the rigid-body motion.
 * @param center_acc Translational acceleration at center.
 * @param omega Angular velocity.
 * @param alpha Angular acceleration.
 * @param consider_point_masses Whether point-mass features are included.
 */
void LoadCollector::add_inertialload(model::ElementRegion::Ptr region,
                                     Vec3 center,
                                     Vec3 center_acc,
                                     Vec3 omega,
                                     Vec3 alpha,
                                     bool consider_point_masses) {
    auto inertial_load = std::make_shared<InertialLoad>();

    // Transfer the selected region and copy the kinematic parameters that define
    // acceleration at every spatial point.
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
