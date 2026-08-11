/**
 * @file contact.cpp
 * @brief Implements frictionless node-to-surface penalty contact.
 *
 * Slave surface nodes are assigned representative tributary areas. For each
 * slave node every master facet is tested directly using the bounded
 * closest-point projection provided by the surface geometry classes. The nearest
 * master facet owns the current contact evaluation. No spatial acceleration,
 * mortar integration, augmented-Lagrange state or formulation switching is
 * present in this implementation.
 *
 * The normal contact law follows the regularized linear node-to-face law used by
 * CalculiX. Around zero gap the linear penalty is multiplied by an arctangent
 * transition, giving a smooth force response and a very small tensile branch on
 * the open side. Contact evaluation is retained up to a small positive cutoff
 * proportional to the slave nodal characteristic length.
 *
 * The tangent follows the analytic node-to-face linearization used by CalculiX.
 * The selected master face is fixed during one tangent evaluation, while the
 * closest-point coordinates, shape functions, surface normal, gap and force law
 * are differentiated consistently from the closest-point orthogonality
 * equations. The frictionless tangent is symmetrized after local assembly.
 *
 * Compact diagnostics report the current active contact set, gap range, force
 * range and local closest-point tangent quality. Pair-level diagnostics can be
 * disabled independently from the one-line assembly summary below.
 *
 * @see Contact
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 11.08.2026
 */

#include "contact.h"

#include "../../core/logging.h"
#include "../../model/geometry/surface/surface_interface.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fem::constraint {
namespace {

constexpr Precision geometry_tolerance     = Precision(1e-12);
constexpr Precision tangent_tolerance      = Precision(1e-14);
constexpr Precision smoothing_length_ratio = Precision(1e-6);
constexpr Precision contact_cutoff_ratio   = Precision(1e-3);
constexpr Precision pi                     = Precision(3.141592653589793238462643383279502884L);

// Temporary diagnostics for the reduced node-to-surface contact formulation.
// The normal global logging switch still suppresses these messages.
constexpr bool contact_diagnostics      = true;
constexpr bool contact_pair_diagnostics = true;

/**
 * @brief Closest bounded projection of one slave node onto one master face.
 */
struct MasterProjection {
    ID                           surface_id = ID(-1);
    model::SurfaceInterface::Ptr surface;
    Vec2                         local       = Vec2::Zero();
    Vec3                         point       = Vec3::Zero();
    Precision                    distance_sq = std::numeric_limits<Precision>::max();
};

/**
 * @brief Current force state of one slave-node/master-face contact pair.
 */
struct PairEvaluation {
    bool            valid        = false;
    bool            active       = false;
    Precision       gap          = Precision(0);
    Precision       force_scalar = Precision(0);
    Vec2            local        = Vec2::Zero();
    Vec3            relative     = Vec3::Zero();
    Vec3            normal       = Vec3::Zero();
    Vec3            slave_force  = Vec3::Zero();
    DynamicVector   shape;
    std::vector<ID> node_ids;
};

bool valid_surface_id(
    ID                                                surface_id,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    return surface_id >= 0 &&
           static_cast<std::size_t>(surface_id) < surfaces.size() &&
           static_cast<bool>(surfaces[static_cast<std::size_t>(surface_id)]);
}

bool surface_contains_node(const model::SurfaceInterface& surface, ID node_id) {
    for (Index local_node = 0; local_node < surface.n_nodes; ++local_node) {
        if (surface.nodes()[local_node] == node_id) {
            return true;
        }
    }
    return false;
}

/**
 * @brief CalculiX-like positive tributary-area weight of one slave facet node.
 */
Precision slave_area_weight(Index n_nodes, Index local_node) {
    switch (n_nodes) {
        case 3:
            return Precision(1) / Precision(3);
        case 4:
            return Precision(1) / Precision(4);
        case 6:
            return local_node < 3 ? Precision(0) : Precision(1) / Precision(3);
        case 8:
            return local_node < 4 ? Precision(0.01) : Precision(0.24);
        default:
            return Precision(0);
    }
}

Precision slave_characteristic_length(Precision slave_area) {
    return std::sqrt(std::max(slave_area, Precision(0)));
}

Precision smoothing_length(Precision slave_area) {
    return smoothing_length_ratio * slave_characteristic_length(slave_area);
}

Precision contact_cutoff(Precision slave_area) {
    return contact_cutoff_ratio * slave_characteristic_length(slave_area);
}

/**
 * @brief CalculiX-style arctangent transition multiplying the linear penalty.
 */
Precision smooth_contact_factor(Precision gap, Precision slave_area) {
    const Precision eps = smoothing_length(slave_area);
    if (!(eps > geometry_tolerance)) {
        return gap < Precision(0) ? Precision(1) : Precision(0);
    }

    return Precision(0.5) + std::atan(-gap / eps) / pi;
}

/**
 * @brief Derivative of the regularized scalar slave force with respect to gap.
 *
 * For
 *
 *     q(g) = K A g [1/2 + atan(-g/eps)/pi]
 *
 * this returns dq/dg with the representative slave area held fixed, matching
 * the CalculiX contact-spring linearization.
 */
Precision smooth_force_derivative(
    Precision gap,
    Precision slave_area,
    Precision penalty
) {
    const Precision eps = smoothing_length(slave_area);
    if (!(eps > geometry_tolerance)) {
        return gap < Precision(0) ? penalty * slave_area : Precision(0);
    }

    const Precision factor = smooth_contact_factor(gap, slave_area);
    const Precision ratio  = gap / eps;
    const Precision df_dg  = -Precision(1) / (pi * eps * (Precision(1) + ratio * ratio));

    return penalty * slave_area * (factor + gap * df_dg);
}

std::unordered_map<ID, Precision> build_slave_nodal_areas(
    const model::SurfaceRegion&                       slave_region,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                               node_coords
) {
    std::unordered_map<ID, Precision> areas;

    for (ID surface_id : slave_region) {
        if (!valid_surface_id(surface_id, surfaces)) {
            continue;
        }

        const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];
        const Precision area = surface->area(node_coords);

        if (!std::isfinite(area) || area <= geometry_tolerance) {
            continue;
        }

        for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
            const Precision weight = slave_area_weight(surface->n_nodes, local_node);
            if (weight > Precision(0)) {
                areas[surface->nodes()[local_node]] += weight * area;
            }
        }
    }

    return areas;
}

std::optional<MasterProjection> find_master_projection(
    ID                                                slave_node,
    const Vec3&                                       slave_position,
    const model::SurfaceRegion&                       master_region,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                               node_coords
) {
    std::optional<MasterProjection> best;

    for (ID surface_id : master_region) {
        if (!valid_surface_id(surface_id, surfaces)) {
            continue;
        }

        const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];
        if (surface_contains_node(*surface, slave_node)) {
            continue;
        }

        const Vec2 local = surface->global_to_local(slave_position, node_coords, true);
        if (!local.allFinite() || !surface->in_bounds(local)) {
            continue;
        }

        const Vec3 point = surface->local_to_global(local, node_coords);
        if (!point.allFinite()) {
            continue;
        }

        const Precision distance_sq = (slave_position - point).squaredNorm();
        if (!std::isfinite(distance_sq)) {
            continue;
        }

        const bool better =
            !best.has_value() ||
            distance_sq < best->distance_sq - geometry_tolerance ||
            (std::abs(distance_sq - best->distance_sq) <= geometry_tolerance &&
             surface_id < best->surface_id);

        if (better) {
            MasterProjection projection;
            projection.surface_id  = surface_id;
            projection.surface     = surface;
            projection.local       = local;
            projection.point       = point;
            projection.distance_sq = distance_sq;
            best = std::move(projection);
        }
    }

    return best;
}

PairEvaluation evaluate_pair(
    ID                                  slave_node,
    Precision                           slave_area,
    const model::SurfaceInterface::Ptr& master,
    model::Field&                       node_coords,
    Precision                           penalty,
    Precision                           clearance,
    bool                                flip_normal
) {
    PairEvaluation result;

    if (!master || slave_area <= Precision(0)) {
        return result;
    }

    const Vec3 slave_position = node_coords.row_vec3(slave_node);
    const Vec2 local = master->global_to_local(slave_position, node_coords, true);
    if (!local.allFinite() || !master->in_bounds(local)) {
        return result;
    }

    const Vec3 master_position = master->local_to_global(local, node_coords);
    Vec3 normal = master->normal(node_coords, local);
    if (flip_normal) {
        normal = -normal;
    }

    const DynamicVector shape = master->shape_function(local);
    if (!slave_position.allFinite() || !master_position.allFinite() ||
        !normal.allFinite() || normal.norm() <= geometry_tolerance ||
        shape.size() != master->n_nodes || !shape.allFinite()) {
        return result;
    }
    normal.normalize();

    const Vec3 relative = slave_position - master_position;
    const Precision gap = relative.dot(normal) - clearance;
    if (!std::isfinite(gap)) {
        return result;
    }

    result.valid        = true;
    result.active       = gap <= contact_cutoff(slave_area);
    result.gap          = gap;
    result.local        = local;
    result.relative     = relative;
    result.normal       = normal;
    result.shape        = shape;
    result.force_scalar = penalty * slave_area * gap * smooth_contact_factor(gap, slave_area);
    result.slave_force  = result.force_scalar * normal;

    result.node_ids.reserve(static_cast<std::size_t>(master->n_nodes + 1));
    result.node_ids.push_back(slave_node);
    for (Index local_node = 0; local_node < master->n_nodes; ++local_node) {
        result.node_ids.push_back(master->nodes()[local_node]);
    }

    return result;
}

/**
 * @brief Builds the CalculiX-style analytic tangent for one fixed master face.
 *
 * The closest point satisfies
 *
 *     r . x_r = 0,    r . x_s = 0.
 *
 * Differentiating these two equations gives d(r,s)/du from a 2x2 system. The
 * resulting parameter derivatives are then propagated through the relative
 * vector, both surface tangents, the normalized normal, the normal gap and the
 * master shape functions. The representative slave area and discrete master
 * face choice are held fixed during this linearization, as in CalculiX.
 *
 * `determinant_out` exposes the determinant of the closest-point 2x2 system for
 * diagnostics. `valid_out` is false whenever the local linearization cannot be
 * formed safely.
 */
DynamicMatrix analytic_pair_tangent(
    const PairEvaluation&               pair,
    Precision                           slave_area,
    const model::SurfaceInterface::Ptr& master,
    const model::Field&                 node_coords,
    Precision                           penalty,
    bool                                flip_normal,
    Precision&                          determinant_out,
    bool&                               valid_out
) {
    const Index master_nodes = master->n_nodes;
    const Index local_nodes  = master_nodes + 1;
    const Index local_dofs   = Index(3) * local_nodes;

    determinant_out = Precision(0);
    valid_out       = false;

    DynamicMatrix tangent = DynamicMatrix::Zero(local_dofs, local_dofs);

    const DynamicMatrix dshape  = master->shape_derivative(pair.local);
    const DynamicMatrix ddshape = master->shape_second_derivative(pair.local);
    if (dshape.rows() != master_nodes || dshape.cols() != 2 ||
        ddshape.rows() != master_nodes || ddshape.cols() != 3 ||
        !dshape.allFinite() || !ddshape.allFinite()) {
        return tangent;
    }

    DynamicMatrix coordinates(master_nodes, 3);
    for (Index node = 0; node < master_nodes; ++node) {
        coordinates.row(node) = node_coords.row_vec3(master->nodes()[node]).transpose();
    }

    const Vec3 xr  = coordinates.transpose() * dshape.col(0);
    const Vec3 xs  = coordinates.transpose() * dshape.col(1);
    const Vec3 xrr = coordinates.transpose() * ddshape.col(0);
    const Vec3 xss = coordinates.transpose() * ddshape.col(1);
    const Vec3 xrs = coordinates.transpose() * ddshape.col(2);

    const Precision orientation = flip_normal ? Precision(-1) : Precision(1);
    const Vec3 normal_vector = orientation * xr.cross(xs);
    const Precision normal_length = normal_vector.norm();
    if (!(normal_length > geometry_tolerance) || !normal_vector.allFinite()) {
        return tangent;
    }
    const Vec3 normal = normal_vector / normal_length;

    // Hessian of the closest-point orthogonality equations. This is the same
    // 2x2 system denoted a11/a12/a22 in CalculiX springstiff_n2f.f.
    const Precision a11 = -xr.dot(xr) + pair.relative.dot(xrr);
    const Precision a12 = -xr.dot(xs) + pair.relative.dot(xrs);
    const Precision a22 = -xs.dot(xs) + pair.relative.dot(xss);
    const Precision determinant = a11 * a22 - a12 * a12;
    determinant_out = determinant;

    const Precision determinant_scale = std::max({
        Precision(1), std::abs(a11 * a22), a12 * a12
    });
    if (!std::isfinite(determinant) ||
        std::abs(determinant) <= geometry_tolerance * determinant_scale) {
        return tangent;
    }

    const Precision c11 =  a22 / determinant;
    const Precision c12 = -a12 / determinant;
    const Precision c22 =  a11 / determinant;
    const Precision dq_dg = smooth_force_derivative(pair.gap, slave_area, penalty);

    for (Index col_node = 0; col_node < local_nodes; ++col_node) {
        const bool slave_column = col_node == 0;
        const Index master_col  = col_node - 1;

        for (Dim component = 0; component < 3; ++component) {
            Vec3 direction = Vec3::Zero();
            direction(component) = Precision(1);

            Precision b1;
            Precision b2;

            if (slave_column) {
                b1 = -xr(component);
                b2 = -xs(component);
            } else {
                b1 = pair.shape(master_col) * xr(component)
                   - dshape(master_col, 0) * pair.relative(component);
                b2 = pair.shape(master_col) * xs(component)
                   - dshape(master_col, 1) * pair.relative(component);
            }

            const Precision dr_local = c11 * b1 + c12 * b2;
            const Precision ds_local = c12 * b1 + c22 * b2;

            Vec3 drelative = -xr * dr_local - xs * ds_local;
            if (slave_column) {
                drelative += direction;
            } else {
                drelative -= pair.shape(master_col) * direction;
            }

            Vec3 dxr = xrr * dr_local + xrs * ds_local;
            Vec3 dxs = xrs * dr_local + xss * ds_local;
            if (!slave_column) {
                dxr += dshape(master_col, 0) * direction;
                dxs += dshape(master_col, 1) * direction;
            }

            const Vec3 dnormal_vector = orientation * (dxr.cross(xs) + xr.cross(dxs));
            const Vec3 dnormal =
                (dnormal_vector - normal * normal.dot(dnormal_vector)) / normal_length;

            const Precision dgap = drelative.dot(normal) + pair.relative.dot(dnormal);
            const Vec3 dslave_force =
                dq_dg * dgap * normal + pair.force_scalar * dnormal;

            const Index col = Index(3) * col_node + component;

            for (Dim row_component = 0; row_component < 3; ++row_component) {
                tangent(row_component, col) = dslave_force(row_component);
            }

            for (Index row_master = 0; row_master < master_nodes; ++row_master) {
                const Precision dN =
                    dshape(row_master, 0) * dr_local +
                    dshape(row_master, 1) * ds_local;

                const Vec3 dmaster_force =
                    -dN * pair.slave_force - pair.shape(row_master) * dslave_force;

                const Index row_node = row_master + 1;
                for (Dim row_component = 0; row_component < 3; ++row_component) {
                    tangent(Index(3) * row_node + row_component, col) = dmaster_force(row_component);
                }
            }
        }
    }

    tangent = Precision(0.5) * (tangent + tangent.transpose()).eval();
    valid_out = tangent.allFinite();
    return tangent;
}

void add_force(model::NodeData& nodal_forces, ID node_id, const Vec3& force) {
    for (Dim component = 0; component < 3; ++component) {
        nodal_forces(node_id, component) += force(component);
    }
}

void add_tangent_entry(
    ID            row_node,
    Dim           row_component,
    ID            col_node,
    Dim           col_component,
    Precision     value,
    SystemDofIds& system_nodal_dofs,
    TripletList&  triplets
) {
    if (std::abs(value) <= tangent_tolerance) {
        return;
    }

    const int row = system_nodal_dofs(row_node, row_component);
    const int col = system_nodal_dofs(col_node, col_component);
    if (row >= 0 && col >= 0) {
        triplets.emplace_back(row, col, value);
    }
}

} // namespace

Contact::Contact(
    model::SurfaceRegion::Ptr master,
    model::SurfaceRegion::Ptr slave,
    Precision                 penalty_stiffness,
    Precision                 contact_clearance,
    bool                      flip_master_normal
)
    : master_surfaces(std::move(master)),
      slave_surfaces(std::move(slave)),
      penalty(penalty_stiffness),
      clearance(contact_clearance),
      flip_normal(flip_master_normal) {
    logging::error(master_surfaces != nullptr && master_surfaces->size() > 0,
        "CONTACT: master surface region is empty or undefined");
    logging::error(slave_surfaces != nullptr && slave_surfaces->size() > 0,
        "CONTACT: slave surface region is empty or undefined");
    logging::error(std::isfinite(penalty) && penalty > Precision(0),
        "CONTACT: PENALTY must be finite and positive");
    logging::error(std::isfinite(clearance),
        "CONTACT: CLEARANCE must be finite");
}

void Contact::assemble(
    SystemDofIds&     system_nodal_dofs,
    model::ModelData& model_data,
    model::NodeData&  nodal_forces,
    TripletList*      triplets
) const {
    logging::error(master_surfaces != nullptr && slave_surfaces != nullptr,
        "CONTACT: surface regions are not defined");
    logging::error(model_data.positions != nullptr,
        "CONTACT: POSITION field is not set");

    model::Field& node_coords = *model_data.positions;
    auto& surfaces = model_data.surfaces;

    const auto slave_areas = build_slave_nodal_areas(*slave_surfaces, surfaces, node_coords);

    Index projected_count       = 0;
    Index valid_count           = 0;
    Index active_count          = 0;
    Index tangent_failure_count = 0;

    Precision min_gap       = std::numeric_limits<Precision>::max();
    Precision max_gap       = -std::numeric_limits<Precision>::max();
    Precision max_penetration = Precision(0);
    Precision max_force       = Precision(0);
    Vec3 slave_resultant       = Vec3::Zero();

    ID worst_slave = ID(-1);
    ID worst_face  = ID(-1);
    Vec2 worst_local = Vec2::Zero();
    Precision worst_area  = Precision(0);
    Precision worst_gap   = Precision(0);
    Precision worst_force = Precision(0);

    for (const auto& [slave_node, slave_area] : slave_areas) {
        if (slave_area <= Precision(0)) {
            continue;
        }

        const Vec3 slave_position = node_coords.row_vec3(slave_node);
        const auto projection = find_master_projection(
            slave_node,
            slave_position,
            *master_surfaces,
            surfaces,
            node_coords
        );
        if (!projection.has_value()) {
            continue;
        }
        ++projected_count;

        const PairEvaluation pair = evaluate_pair(
            slave_node,
            slave_area,
            projection->surface,
            node_coords,
            penalty,
            clearance,
            flip_normal
        );
        if (!pair.valid) {
            continue;
        }
        ++valid_count;

        min_gap = std::min(min_gap, pair.gap);
        max_gap = std::max(max_gap, pair.gap);

        if (!pair.active) {
            continue;
        }
        ++active_count;

        const Precision penetration = std::max(Precision(0), -pair.gap);
        const Precision force_abs   = std::abs(pair.force_scalar);
        max_penetration = std::max(max_penetration, penetration);
        max_force       = std::max(max_force, force_abs);
        slave_resultant += pair.slave_force;

        if (penetration >= std::max(Precision(0), -worst_gap)) {
            worst_slave = slave_node;
            worst_face  = projection->surface_id;
            worst_local = pair.local;
            worst_area  = slave_area;
            worst_gap   = pair.gap;
            worst_force = pair.force_scalar;
        }

        add_force(nodal_forces, pair.node_ids[0], pair.slave_force);
        for (Index master_node = 0; master_node < projection->surface->n_nodes; ++master_node) {
            add_force(
                nodal_forces,
                pair.node_ids[static_cast<std::size_t>(master_node + 1)],
                -pair.shape(master_node) * pair.slave_force
            );
        }

        Precision tangent_determinant = Precision(0);
        bool tangent_valid = true;
        DynamicMatrix tangent;

        if (triplets != nullptr) {
            tangent = analytic_pair_tangent(
                pair,
                slave_area,
                projection->surface,
                node_coords,
                penalty,
                flip_normal,
                tangent_determinant,
                tangent_valid
            );

            if (!tangent_valid) {
                ++tangent_failure_count;
            }

            logging::warning(
                tangent_valid,
                "CONTACT tangent invalid: slave=", slave_node,
                " face=", projection->surface_id,
                " local=(`", pair.local(0), ", ", pair.local(1), ")",
                " gap=", pair.gap,
                " area=", slave_area,
                " detH=", tangent_determinant
            );

            logging::info(
                contact_pair_diagnostics,
                "CONTACT pair: slave=", slave_node,
                " face=", projection->surface_id,
                " local=(`", pair.local(0), ", ", pair.local(1), ")",
                " gap=", pair.gap,
                " area=", slave_area,
                " fn=", pair.force_scalar,
                " detH=", tangent_determinant,
                " |Kt|=", tangent.norm()
            );

            const Index local_nodes = static_cast<Index>(pair.node_ids.size());
            for (Index row_node = 0; row_node < local_nodes; ++row_node) {
                for (Index col_node = 0; col_node < local_nodes; ++col_node) {
                    for (Dim row_component = 0; row_component < 3; ++row_component) {
                        for (Dim col_component = 0; col_component < 3; ++col_component) {
                            add_tangent_entry(
                                pair.node_ids[static_cast<std::size_t>(row_node)],
                                row_component,
                                pair.node_ids[static_cast<std::size_t>(col_node)],
                                col_component,
                                tangent(Index(3) * row_node + row_component,
                                        Index(3) * col_node + col_component),
                                system_nodal_dofs,
                                *triplets
                            );
                        }
                    }
                }
            }
        } else {
            logging::info(
                contact_pair_diagnostics,
                "CONTACT pair: slave=", slave_node,
                " face=", projection->surface_id,
                " local=(`", pair.local(0), ", ", pair.local(1), ")",
                " gap=", pair.gap,
                " area=", slave_area,
                " fn=", pair.force_scalar
            );
        }
    }

    const bool have_valid_gap = valid_count > 0;
    logging::info(
        contact_diagnostics,
        "CONTACT summary: slaves=", slave_areas.size(),
        " projected=", projected_count,
        " valid=", valid_count,
        " active=", active_count,
        " min_gap=", have_valid_gap ? min_gap : Precision(0),
        " max_gap=", have_valid_gap ? max_gap : Precision(0),
        " max_pen=", max_penetration,
        " max|fn|=", max_force,
        " slave_resultant=(`", slave_resultant(0), ", ", slave_resultant(1), ", ", slave_resultant(2), ")",
        " tangent_failures=", tangent_failure_count,
        " tangent=", triplets != nullptr ? "yes" : "no"
    );

    if (active_count > 0) {
        logging::info(
            contact_diagnostics,
            "CONTACT worst: slave=", worst_slave,
            " face=", worst_face,
            " local=(`", worst_local(0), ", ", worst_local(1), ")",
            " gap=", worst_gap,
            " area=", worst_area,
            " fn=", worst_force
        );
    }
}

} // namespace fem::constraint
