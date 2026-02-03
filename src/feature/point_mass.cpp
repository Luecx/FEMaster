#include "point_mass.h"

namespace fem { namespace feature {

static inline void add_diag_triplet(int gid, Precision v, TripletList& out) {
    if (gid >= 0 && v != Precision(0)) {
        out.emplace_back(gid, gid, v);
    }
}

void PointMass::assemble_stiffness(const SystemDofIds& indices, TripletList& out) const {
    if (!region) return;
    for (const auto node : *region) {
        // translational springs kx, ky, kz
        add_diag_triplet(indices(node, 0), spring_constants(0), out);
        add_diag_triplet(indices(node, 1), spring_constants(1), out);
        add_diag_triplet(indices(node, 2), spring_constants(2), out);
        // rotational springs krx, kry, krz
        add_diag_triplet(indices(node, 3), rotary_spring_constants(0), out);
        add_diag_triplet(indices(node, 4), rotary_spring_constants(1), out);
        add_diag_triplet(indices(node, 5), rotary_spring_constants(2), out);
    }
}

void PointMass::assemble_mass(const SystemDofIds& indices, TripletList& out) const {
    if (!region) return;
    for (const auto node : *region) {
        // translational mass (same m on x,y,z)
        add_diag_triplet(indices(node, 0), mass, out);
        add_diag_triplet(indices(node, 1), mass, out);
        add_diag_triplet(indices(node, 2), mass, out);
        // rotary inertia on rotational DOFs
        add_diag_triplet(indices(node, 3), rotary_inertia(0), out);
        add_diag_triplet(indices(node, 4), rotary_inertia(1), out);
        add_diag_triplet(indices(node, 5), rotary_inertia(2), out);
    }
}

} } // namespace fem::feature

