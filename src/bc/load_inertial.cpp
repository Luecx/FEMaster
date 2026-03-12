/**
 * @file load_inertial.cpp
 * @brief Implements the rigid-body inertial load boundary condition.
 *
 * @see src/bc/load_inertial.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_inertial.h"

#include "../core/logging.h"
#include "../feature/point_mass.h"
#include "../model/element/element_structural.h"
#include "../model/model_data.h"

#include <sstream>

namespace fem {
namespace bc {

namespace {

void apply_point_mass_inertial_contribution(model::ModelData& model_data,
                                            model::Field& bc,
                                            const Vec3& center,
                                            const Vec3& center_acc,
                                            const Vec3& omega,
                                            const Vec3& alpha) {
    logging::error(model_data.positions != nullptr,
                   "InertialLoad (point masses): positions field not set in model data");
    const auto& positions = *model_data.positions;

    for (const auto& feature_ptr : model_data.features) {
        const auto* pm = dynamic_cast<const feature::PointMass*>(feature_ptr.get());
        if (!pm || !pm->region) {
            continue;
        }

        for (const ID node_id : *pm->region) {
            logging::error(node_id >= 0 && static_cast<Index>(node_id) < positions.rows,
                           "InertialLoad (point masses): node id ", node_id,
                           " out of bounds for positions field with ", positions.rows, " rows");

            const Index node = static_cast<Index>(node_id);
            const Vec3 x = positions.row_vec3(node);
            const Vec3 r = x - center;
            const Vec3 a_rot = alpha.cross(r) + omega.cross(omega.cross(r));

            const Vec3 dF = -pm->mass * (center_acc + a_rot);
            bc(node, 0) += dF(0);
            bc(node, 1) += dF(1);
            bc(node, 2) += dF(2);

            if (bc.components >= 6) {
                const Vec3 Jalpha = pm->rotary_inertia.cwiseProduct(alpha);
                const Vec3 Jomega = pm->rotary_inertia.cwiseProduct(omega);
                const Vec3 dM = -(Jalpha + omega.cross(Jomega));
                bc(node, 3) += dM(0);
                bc(node, 4) += dM(1);
                bc(node, 5) += dM(2);
            }
        }
    }
}

} // namespace

/**
 * @copydoc InertialLoad::apply
 */
void InertialLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    (void)time;
    logging::error(region != nullptr, "InertialLoad: region not set");

    for (auto& el_id : *region) {
        auto& el_ptr = model_data.elements[el_id];
        if (!el_ptr) {
            continue;
        }

        auto structural = el_ptr->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        auto f = [c = center, a0 = center_acc, w = omega, al = alpha](const Vec3& x) -> Vec3 {
            const Vec3 r = x - c;
            const Vec3 a_rot = al.cross(r) + w.cross(w.cross(r));
            return -(a0 + a_rot);
        };
        structural->integrate_vec_field(bc, /*scale_by_density=*/true, f);
    }

    if (consider_point_masses) {
        apply_point_mass_inertial_contribution(model_data, bc, center, center_acc, omega, alpha);
    }
}

/**
 * @copydoc InertialLoad::str
 */
std::string InertialLoad::str() const {
    std::ostringstream os;
    os << "INERTIAL: target=ELSET " << (region ? region->name : std::string("?"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")"
       << ", center=[" << center(0) << ", " << center(1) << ", " << center(2) << "]"
       << ", a0=[" << center_acc(0) << ", " << center_acc(1) << ", " << center_acc(2) << "]"
       << ", omega=[" << omega(0) << ", " << omega(1) << ", " << omega(2) << "]"
       << ", alpha=[" << alpha(0) << ", " << alpha(1) << ", " << alpha(2) << "]"
       << ", consider_point_masses=" << (consider_point_masses ? "true" : "false");
    return os.str();
}

} // namespace bc
} // namespace fem
