/**
 * @file load_inertial.cpp
 * @brief Implements equivalent forces for prescribed rigid-body acceleration.
 */

#include "load_inertial.h"

#include "../../core/logging.h"
#include "../../model/element/element_structural.h"
#include "../../model/element/point.h"
#include "../../model/model_data.h"
#include "../../section/section_point_mass.h"

#include <Eigen/Geometry>
#include <sstream>

namespace fem::bc {
namespace {

void apply_point_element_inertial_contribution(const model::PointElement& point,
                                               const model::Field&        positions,
                                               model::Field&              rhs,
                                               const Vec3&                center,
                                               const Vec3&                center_acc,
                                               const Vec3&                omega,
                                               const Vec3&                alpha) {
    if (!point._section) return;

    const auto* section = point._section->as<PointMassSection>();
    logging::error(section != nullptr,
        "InertialLoad: point element ", point.elem_id, " has a non-point-mass section");

    const ID node_id = point.nodes()[0];
    logging::error(node_id >= 0 && static_cast<Index>(node_id) < positions.rows,
        "InertialLoad: point-element node ", node_id,
        " is outside positions field with ", positions.rows, " rows");

    const Index node  = static_cast<Index>(node_id);
    const Vec3  x     = positions.row_vec3(node);
    const Vec3  r     = x - center;
    const Vec3  a_rot = alpha.cross(r) + omega.cross(omega.cross(r));
    const Vec3  dF    = -section->mass_ * (center_acc + a_rot);

    rhs(node, 0) += dF(0);
    rhs(node, 1) += dF(1);
    rhs(node, 2) += dF(2);

    if (rhs.components >= 6) {
        const Vec3 Jalpha = section->rotary_inertia_.cwiseProduct(alpha);
        const Vec3 Jomega = section->rotary_inertia_.cwiseProduct(omega);
        const Vec3 dM     = -(Jalpha + omega.cross(Jomega));

        rhs(node, 3) += dM(0);
        rhs(node, 4) += dM(1);
        rhs(node, 5) += dM(2);
    }
}

} // namespace

void InertialLoad::apply(model::ModelData& model_data,
                         model::Field&     rhs,
                         Precision         time,
                         bool              ignore_amplitude) {
    (void) time;
    (void) ignore_amplitude;

    logging::error(region_ != nullptr,
        "InertialLoad: region not set");
    logging::error(model_data.positions != nullptr,
        "InertialLoad: positions field not set in model data");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "InertialLoad: target field must be a NODE field with at least six components");

    const auto& positions = *model_data.positions;

    for (const ID element_id : *region_) {
        logging::error(element_id >= 0
                    && static_cast<std::size_t>(element_id) < model_data.elements.size(),
            "InertialLoad: element ", element_id, " is outside compiled element storage");

        auto& element = model_data.elements[static_cast<std::size_t>(element_id)];
        if (!element) continue;

        if (auto* point = element->as<model::PointElement>()) {
            apply_point_element_inertial_contribution(
                *point, positions, rhs, center_, center_acc_, omega_, alpha_);
            continue;
        }

        auto* structural = element->as<model::StructuralElement>();
        if (!structural) continue;

        auto acceleration = [c = center_, a0 = center_acc_, w = omega_, al = alpha_](const Vec3& x) -> Vec3 {
            const Vec3 r = x - c;
            return -(a0 + al.cross(r) + w.cross(w.cross(r)));
        };
        structural->integrate_vector_field(rhs, true, acceleration);
    }

    if (consider_point_masses_) {
        for (const auto& element : model_data.point_elements) {
            if (!element) continue;
            const auto* point = element->as<model::PointElement>();
            logging::error(point != nullptr,
                "InertialLoad: auxiliary point-element storage contains a non-point element");
            apply_point_element_inertial_contribution(
                *point, positions, rhs, center_, center_acc_, omega_, alpha_);
        }
    }
}

std::string InertialLoad::str() const {
    std::ostringstream os;
    os << "INERTIAL: target=ELSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", center=[" << center_(0) << ", " << center_(1) << ", " << center_(2) << "]"
       << ", a0=[" << center_acc_(0) << ", " << center_acc_(1) << ", " << center_acc_(2) << "]"
       << ", omega=[" << omega_(0) << ", " << omega_(1) << ", " << omega_(2) << "]"
       << ", alpha=[" << alpha_(0) << ", " << alpha_(1) << ", " << alpha_(2) << "]"
       << ", consider_point_masses=" << (consider_point_masses_ ? "true" : "false");
    return os.str();
}

} // namespace fem::bc
