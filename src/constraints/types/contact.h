/**
 * @file contact.h
 * @brief Declares a simple frictionless node-to-surface penalty contact.
 *
 * Contact uses the same surface/node set resolution style as tie constraints,
 * but contributes nonlinear internal force and tangent stiffness instead of
 * fixed constraint equations.
 */

#pragma once

#include "../../data/field.h"
#include "../../data/region.h"
#include "../../model/geometry/surface/surface.h"

namespace fem {
namespace model {
struct ModelData;
}

namespace constraint {

class Contact {
    model::SurfaceRegion::Ptr master_surfaces;
    model::NodeRegion::Ptr    slave_nodes;
    model::SurfaceRegion::Ptr slave_surfaces;

    Precision distance;
    Precision penalty;
    Precision clearance;
    bool flip_normal;

public:
    Contact(model::SurfaceRegion::Ptr master,
            model::NodeRegion::Ptr slave,
            Precision search_distance,
            Precision penalty_stiffness,
            Precision contact_clearance,
            bool flip_master_normal);

    Contact(model::SurfaceRegion::Ptr master,
            model::SurfaceRegion::Ptr slave,
            Precision search_distance,
            Precision penalty_stiffness,
            Precision contact_clearance,
            bool flip_master_normal);

    void assemble(SystemDofIds& system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData& nodal_forces,
                  TripletList& triplets) const;
};

} // namespace constraint
} // namespace fem
