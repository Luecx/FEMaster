/**
 * @file rigid_body.h
 * @brief Abaqus-style rigid body constraint with element, PIN and TIE regions.
 */

#pragma once

#include "equation.h"
#include "../../data/region.h"

#include <utility>
#include <vector>

namespace fem::model { struct ModelData; }

namespace fem::constraint {

/**
 * A rigid body is one reference node and up to three independent regions.
 *
 * The input regions are retained as supplied. prepare() classifies the
 * participating elements and builds the effective PIN and TIE node regions.
 * Explicit PIN/TIE regions override element-derived assignments.
 *
 * The mechanical constraint owns its own equations; no Coupling objects are
 * inserted into ModelData by the input parser.
 */
class RigidBody {
public:
    enum class Position { Input, CenterOfMass };

    ID reference_node = -1;
    model::ElementRegion::Ptr element_region = nullptr;
    model::NodeRegion::Ptr pin_region = nullptr;
    model::NodeRegion::Ptr tie_region = nullptr;
    Position position = Position::Input;

    RigidBody() = default;

    RigidBody(ID reference,
              model::ElementRegion::Ptr elements,
              model::NodeRegion::Ptr pin,
              model::NodeRegion::Ptr tie,
              Position location = Position::Input)
        : reference_node(reference),
          element_region(std::move(elements)),
          pin_region(std::move(pin)),
          tie_region(std::move(tie)),
          position(location) {}

    // Resolve element classifications, validate memberships and position RP.
    // Called once by Model::add_rigid_body() after assembly compilation.
    void prepare(model::ModelData& data);

    Equations get_equations(SystemDofIds& indices, model::ModelData& data) const;
    Dofs master_dofs(SystemDofs& mask, model::ModelData& data) const;

    const std::vector<ID>& rigid_elements() const { return elements_; }
    const model::NodeRegion::Ptr& pin_nodes() const { return resolved_pin_; }
    const model::NodeRegion::Ptr& tie_nodes() const { return resolved_tie_; }

private:
    std::vector<ID> elements_;
    model::NodeRegion::Ptr resolved_pin_ = nullptr;
    model::NodeRegion::Ptr resolved_tie_ = nullptr;
};

} // namespace fem::constraint
