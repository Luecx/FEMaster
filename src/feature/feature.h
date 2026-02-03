/**
 * @file feature.h
 * @brief Base interface for non-element features contributing to global operators.
 */

#pragma once

#include "../core/types_eig.h"

namespace fem { namespace feature {

struct Feature {
    using Ptr = std::shared_ptr<Feature>;
    virtual ~Feature() = default;

    // Minimal-invasive: only assemble onto existing active DOFs (indices >= 0).
    virtual void assemble_stiffness(const SystemDofIds& indices, TripletList& out) const = 0;
    virtual void assemble_mass     (const SystemDofIds& indices, TripletList& out) const = 0;
};

} } // namespace fem::feature

