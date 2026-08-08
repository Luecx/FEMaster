/**
 * @file contact.cpp
 * @brief Implements frictionless contact against a reconstructed Nagata master surface.
 *
 * Every assembly reconstructs the selected master region from the current nodal
 * coordinate field. New slave points use global Nagata projection, while
 * established points continue from their transactional locations and walk
 * continuously across internal triangular charts.
 *
 * The augmented-Lagrange law operates on the signed Nagata normal gap. Its
 * slave gradient is the analytical closest-point normal. Master derivatives
 * are evaluated from coordinate differences of the reconstructed position at
 * the stationary surface location; the closest-point envelope property avoids
 * differentiating or repeating the projection. Residual assembly uses this
 * variational gap gradient and optional tangent assembly retains the symmetric
 * positive-semidefinite Gauss-Newton term.
 *
 * @see Contact
 * @see model::NagataSurface
 *
 * @author Finn Eggers
 * @date 08.08.2026
 */

#include "contact.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fem::constraint {
namespace {

// The augmentation tolerance is tied to the originating master-surface size.
// Multipliers inside this geometric band remain unchanged so the outer loop can
// terminate without oscillating around a zero gap.
constexpr Precision augmentation_gap_relative_tolerance        = Precision(1e-4);
constexpr Precision augmentation_gap_absolute_tolerance        = Precision(1e-10);
constexpr Precision augmentation_multiplier_relative_tolerance = Precision(1e-6);
constexpr Precision augmentation_multiplier_absolute_tolerance = Precision(1e-10);

constexpr bool print_contact_summary = false;

/**
 * @brief Diagnostic quantities collected during one Nagata contact assembly.
 *
 * The counters distinguish global acquisition from tracked continuation and
 * describe changes of physical connected master components. Patch changes are
 * intentionally absent because patches are merely local coordinate charts.
 */
struct ContactDiagnostics {
    Index slave_nodes             = 0;
    Index global_projections      = 0;
    Index tracked_projections     = 0;
    Index no_partner              = 0;
    Index open_closest_partner    = 0;
    Index active_contacts         = 0;
    Index activations             = 0;
    Index deactivations           = 0;
    Index component_switches      = 0;
    Index self_contact_rejections = 0;
    Index active_partner_losses   = 0;
    Index approximate_tangents    = 0;

    Precision maximum_penetration      = Precision(0);
    Precision maximum_closest_distance = Precision(0);
    Precision maximum_slave_force      = Precision(0);
    Precision slave_force_squared_sum  = Precision(0);
    Precision maximum_tangent_norm     = Precision(0);

    ID                         worst_slave     = ID(-1);
    model::nagata::ComponentID worst_component = std::numeric_limits<model::nagata::ComponentID>::max();
    model::nagata::PatchID     worst_patch     = std::numeric_limits<model::nagata::PatchID>::max();
    Vec2                       worst_local     = Vec2::Zero();
    Precision                  worst_gap       = Precision(0);
    Precision                  worst_distance  = Precision(0);

    std::uint64_t signature = 1469598103934665603ull;
};

/**
 * @brief Local active-set value of the augmented frictionless normal law.
 *
 * `value` is the shifted signed gap used directly in the residual. `pressure`
 * is its non-negative compressive counterpart used for diagnostics. During a
 * frozen active-set continuation, `value` may become positive until the outer
 * contact update releases the slave node.
 */
struct ContactLaw {
    Precision value      = Precision(0);
    Precision derivative = Precision(0);
    Precision pressure   = Precision(0);
    bool      active     = false;
};

/**
 * @brief Appends one integer value to the deterministic contact signature.
 *
 * FNV-1a hashing is applied bytewise so the signature depends only on active
 * slave/component pairs and is independent of hash-container iteration order.
 *
 * @param signature Signature accumulated in sorted slave-node order.
 * @param value Unsigned value to append.
 */
void hash_value(std::uint64_t& signature, std::uint64_t value) {
    constexpr std::uint64_t prime = 1099511628211ull;

    for (Index byte = 0; byte < 8; ++byte) {
        signature ^= value & 0xffull;
        signature *= prime;
        value >>= 8;
    }
}

/**
 * @brief Collects the unique slave nodes represented by the contact definition.
 *
 * A nodal slave region is copied directly. For a surface slave region all
 * surface connectivities are merged. Sorting establishes deterministic assembly
 * and signature order.
 *
 * @param slave_nodes Optional explicit slave-node region.
 * @param slave_surfaces Optional slave-surface region.
 * @param model_data Model topology containing referenced surfaces.
 * @return Sorted unique slave-node identifiers.
 */
std::vector<ID> collect_slave_nodes(const model::NodeRegion::Ptr&    slave_nodes,
                                    const model::SurfaceRegion::Ptr& slave_surfaces,
                                    model::ModelData&                model_data) {
    std::vector<ID> result;

    if (slave_nodes) {
        result.reserve(static_cast<std::size_t>(slave_nodes->size()));

        for (const ID node_id : *slave_nodes) {
            result.push_back(node_id);
        }

        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
    }

    logging::error(slave_surfaces != nullptr,
        "CONTACT requires either a slave node set or a slave surface set");

    std::unordered_set<ID> unique_nodes;
    const auto&            surfaces = model_data.surfaces;

    for (const ID surface_id : *slave_surfaces) {
        logging::error(surface_id >= 0
                    && static_cast<std::size_t>(surface_id) < surfaces.size(),
            "CONTACT slave surface ID ", surface_id,
            " lies outside the surface container");

        const model::SurfaceInterface::Ptr& surface =
            surfaces[static_cast<std::size_t>(surface_id)];

        logging::error(surface != nullptr,
            "CONTACT slave surface ID ", surface_id, " is null");

        for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
            unique_nodes.insert(surface->nodes()[local_node]);
        }
    }

    result.reserve(unique_nodes.size());

    for (const ID node_id : unique_nodes) {
        result.push_back(node_id);
    }

    std::sort(result.begin(), result.end());
    return result;
}

/**
 * @brief Builds positive lumped tributary areas for surface-based slave nodes.
 *
 * Every slave-surface area is divided uniformly among its nodes. Uniform
 * positive lumping is used deliberately because higher-order serendipity shape
 * functions can produce negative consistent nodal area weights.
 *
 * @param slave_surfaces Optional slave-surface region.
 * @param model_data Model topology containing referenced surfaces.
 * @param node_coords Current nodal geometry.
 * @return Accumulated positive area per slave node; empty for nodal slaves.
 */
std::unordered_map<ID, Precision> build_slave_tributary_areas(
    const model::SurfaceRegion::Ptr& slave_surfaces,
    model::ModelData&                model_data,
    const model::Field&              node_coords) {
    std::unordered_map<ID, Precision> areas;

    if (!slave_surfaces) {
        return areas;
    }

    const auto& surfaces = model_data.surfaces;

    for (const ID surface_id : *slave_surfaces) {
        logging::error(surface_id >= 0
                    && static_cast<std::size_t>(surface_id) < surfaces.size(),
            "CONTACT slave surface ID ", surface_id,
            " lies outside the surface container");

        const model::SurfaceInterface::Ptr& surface =
            surfaces[static_cast<std::size_t>(surface_id)];

        logging::error(surface != nullptr && surface->n_nodes > 0,
            "CONTACT slave surface ID ", surface_id,
            " has no valid surface topology");

        const Precision surface_area = surface->area(node_coords);

        logging::error(std::isfinite(surface_area) && surface_area > Precision(0),
            "CONTACT slave surface ID ", surface_id,
            " has invalid area ", surface_area);

        const Precision nodal_area =
            surface_area / static_cast<Precision>(surface->n_nodes);

        for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
            areas[surface->nodes()[local_node]] += nodal_area;
        }
    }

    return areas;
}

/**
 * @brief Evaluates the active branch of the augmented normal contact law.
 *
 * With signed normal gap `g`, multiplier `lambda >= 0` and penalty `k > 0`,
 * the shifted gap is `g_bar = g - lambda / k`. Contact is active for
 * `g_bar <= 0`, producing compressive pressure `-k g_bar`. At equality the
 * pressure vanishes while the normal tangent remains available for an
 * initially closed contact configuration.
 *
 * @param gap Signed Nagata normal gap including clearance.
 * @param normal_multiplier Non-negative augmented normal multiplier.
 * @param penalty Positive penalty stiffness.
 * @return Active-set value, derivative and compressive pressure.
 */
ContactLaw evaluate_augmented_lagrange_law(Precision gap,
                                           Precision normal_multiplier,
                                           Precision penalty) {
    const Precision shifted_gap = gap - normal_multiplier / penalty;

    if (shifted_gap > Precision(0)) {
        return {};
    }

    return {
        shifted_gap,
        Precision(1),
        -penalty * shifted_gap,
        true
    };
}

/**
 * @brief Continues one frozen branch of the augmented normal contact law.
 *
 * Newton and line-search trials must evaluate one differentiable residual. An
 * inherited inactive slave therefore remains inactive even if a trial step
 * penetrates the master surface. An inherited active slave remains on
 *
 *     g_bar = g - lambda / k
 *
 * with unit derivative, even if a trial step temporarily opens the gap. The
 * enclosing contact update subsequently reclassifies the active set from the
 * accepted geometry and restarts Newton when that discrete state changes.
 *
 * @param gap Signed Nagata normal gap including clearance.
 * @param normal_multiplier Non-negative augmented normal multiplier.
 * @param penalty Positive penalty stiffness.
 * @param active Frozen active-set membership inherited by the trial.
 * @return Continued active or inactive branch of the normal contact law.
 */
ContactLaw continue_augmented_lagrange_branch(Precision gap,
                                              Precision normal_multiplier,
                                              Precision penalty,
                                              bool      active) {
    if (!active) {
        return {};
    }

    const Precision shifted_gap = gap - normal_multiplier / penalty;

    return {
        shifted_gap,
        Precision(1),
        std::max(Precision(0), -penalty * shifted_gap),
        true
    };
}

/**
 * @brief Tests whether a source finite-element surface contains one slave node.
 *
 * This prevents a node shared by master and slave definitions from contacting
 * its own source surface. The current global Nagata projection does not search
 * for a second-best source surface after this rejection.
 *
 * @param surface Originating finite-element surface of a Nagata evaluation.
 * @param node_id Slave node to test.
 * @return True if the node occurs in the source connectivity.
 */
bool surface_contains_node(const model::SurfaceInterface& surface, ID node_id) {
    for (Index local_node = 0; local_node < surface.n_nodes; ++local_node) {
        if (surface.nodes()[local_node] == node_id) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Assembles one active variational Nagata contact pair.
 *
 * For a stationary closest point, the derivative of signed distance with
 * respect to the slave position is exactly the oriented master normal. The
 * closest-point envelope property similarly removes the derivative of the
 * projected surface coordinates from every master derivative:
 *
 *     dg / dx_s = n,
 *     dg / du_a = -n^T dx_N / du_a.
 *
 * Only `dx_N / du_a` is currently evaluated by central differences. Both
 * perturbed reconstructions are evaluated at the unperturbed Nagata location;
 * no perturbed projection is necessary. The master stencil includes the source
 * surfaces used to form all three averaged patch-vertex normals.
 *
 * With the discrete gap gradient `B = dg/dq`, the augmented contact potential
 * contributes
 *
 *     r_c = A_s k g_bar B
 *
 * and the symmetric Gauss-Newton tangent
 *
 *     K_c = A_s k (d g_bar / d g) B^T B.
 *
 * This retains the complete first-order Nagata geometry and closest-point
 * sensitivity while omitting only the gap-Hessian term proportional to
 * `g_bar`. The omitted term vanishes at exact contact and the retained matrix
 * remains positive semidefinite for the LDLT solution path.
 *
 * @param node_coords Current global nodal coordinate field.
 * @param master_geometry Nagata reconstruction at the unperturbed state.
 * @param location Tracked closest-point location at the unperturbed state.
 * @param normal Oriented unit normal at the unperturbed closest point.
 * @param slave_node_id Global slave node identifier.
 * @param contact_law Active augmented contact-law state.
 * @param penalty Positive penalty stiffness.
 * @param slave_weight Positive tributary slave area or one for nodal contact.
 * @param characteristic_length Positive local length controlling perturbations.
 * @param system_nodal_dofs Global translational degree-of-freedom mapping.
 * @param nodal_forces Global nodal force accumulator.
 * @param triplets Global tangent triplet accumulator.
 * @param assemble_tangent Whether to assemble the sparse Gauss-Newton block.
 * @param diagnostics Assembly diagnostics updated by this pair.
 */
void assemble_contact_pair(const model::Field&                   node_coords,
                           const model::NagataSurface&           master_geometry,
                           const model::NagataSurface::Location& location,
                           const Vec3&                           normal,
                           ID                                    slave_node_id,
                           const ContactLaw&                     contact_law,
                           Precision                             penalty,
                           Precision                             slave_weight,
                           Precision                             characteristic_length,
                           SystemDofIds&                         system_nodal_dofs,
                           model::NodeData&                      nodal_forces,
                           TripletList&                          triplets,
                           bool                                  assemble_tangent,
                           ContactDiagnostics&                   diagnostics) {
    const Precision contact_scale = slave_weight * penalty;

    // Collect the exact local master stencil before adding the slave coordinate.
    std::vector<ID> dependency_nodes = master_geometry.dependency_nodes(location);
    const std::unordered_set<ID> master_dependency_nodes(
        dependency_nodes.begin(), dependency_nodes.end());

    if (master_dependency_nodes.find(slave_node_id) == master_dependency_nodes.end()) {
        dependency_nodes.push_back(slave_node_id);
    }

    std::sort(dependency_nodes.begin(), dependency_nodes.end());

    const Index local_dofs = static_cast<Index>(dependency_nodes.size()) * 3;

    std::vector<ID>        dof_nodes(static_cast<std::size_t>(local_dofs));
    std::vector<Dim>       dof_components(static_cast<std::size_t>(local_dofs));
    std::vector<int>       global_dofs(static_cast<std::size_t>(local_dofs));
    std::vector<Precision> gap_gradient(static_cast<std::size_t>(local_dofs));

    const Precision difference_step_scale =
        std::cbrt(std::numeric_limits<Precision>::epsilon());

    // Reuse one coordinate field for all perturbations. Every reconstruction
    // observes exactly one changed master coordinate and the baseline value is
    // restored before continuing with the next degree of freedom.
    model::Field perturbed_coords = node_coords;

    // Assemble the analytical slave gradient and the central-difference master
    // position sensitivities. If one node participates on both sides, both
    // contributions are accumulated in the same coordinate derivative.
    Index local_dof = 0;

    for (const ID node_id : dependency_nodes) {
        for (Dim component = 0; component < 3; ++component) {
            const Precision coordinate =
                node_coords(static_cast<Index>(node_id), component);

            const Precision step = difference_step_scale * std::max({
                Precision(1),
                characteristic_length,
                std::abs(coordinate)
            });

            Precision derivative =
                node_id == slave_node_id
                    ? normal(component)
                    : Precision(0);

            if (master_dependency_nodes.find(node_id) != master_dependency_nodes.end()) {
                perturbed_coords(static_cast<Index>(node_id), component) = coordinate + step;

                const Vec3 position_plus =
                    master_geometry.evaluate_position(location, perturbed_coords);

                perturbed_coords(static_cast<Index>(node_id), component) = coordinate - step;

                const Vec3 position_minus =
                    master_geometry.evaluate_position(location, perturbed_coords);

                perturbed_coords(static_cast<Index>(node_id), component) = coordinate;

                derivative -= normal.dot(position_plus - position_minus)
                            / (Precision(2) * step);
            }

            logging::error(std::isfinite(derivative),
                "CONTACT produced a non-finite numerical Nagata gap sensitivity");

            const std::size_t entry = static_cast<std::size_t>(local_dof);

            dof_nodes[entry]      = node_id;
            dof_components[entry] = component;
            global_dofs[entry]    = system_nodal_dofs(node_id, component);
            gap_gradient[entry]   = derivative;

            ++local_dof;
        }
    }

    // Assemble the variational residual of the shifted quadratic contact
    // potential, including reactions at constrained dependency coordinates.
    const Precision residual_scale = contact_scale * contact_law.value;

    for (Index row = 0; row < local_dofs; ++row) {
        const std::size_t entry = static_cast<std::size_t>(row);
        nodal_forces(
            static_cast<Index>(dof_nodes[entry]),
            dof_components[entry]) += residual_scale * gap_gradient[entry];
    }

    if (!assemble_tangent) {
        return;
    }

    // Assemble the symmetric positive-semidefinite Gauss-Newton contact block.
    const Precision tangent_scale = contact_scale * contact_law.derivative;
    Precision       free_gradient_squared_norm = Precision(0);

    for (Index row = 0; row < local_dofs; ++row) {
        const std::size_t row_entry = static_cast<std::size_t>(row);
        const int         global_row = global_dofs[row_entry];

        if (global_row < 0) {
            continue;
        }

        free_gradient_squared_norm += gap_gradient[row_entry] * gap_gradient[row_entry];

        for (Index col = 0; col < local_dofs; ++col) {
            const std::size_t col_entry = static_cast<std::size_t>(col);
            const int         global_col = global_dofs[col_entry];

            if (global_col < 0) {
                continue;
            }

            const Precision value =
                tangent_scale * gap_gradient[row_entry] * gap_gradient[col_entry];

            if (std::abs(value) > Precision(1e-14)) {
                triplets.emplace_back(global_row, global_col, value);
            }
        }
    }

    diagnostics.maximum_tangent_norm = std::max(
        diagnostics.maximum_tangent_norm,
        tangent_scale * free_gradient_squared_norm);

    ++diagnostics.approximate_tangents;
}

} // namespace

/**
 * @brief Constructs contact with an explicit nodal slave region.
 *
 * @param master Master finite-element surfaces reconstructed as one Nagata surface.
 * @param slave Slave node region.
 * @param penalty_stiffness Positive normal penalty stiffness.
 * @param contact_clearance Signed normal clearance subtracted from the gap.
 * @param flip_master_normal Whether to reverse reconstructed master normals.
 */
Contact::Contact(model::SurfaceRegion::Ptr master,
                 model::NodeRegion::Ptr    slave,
                 Precision                 penalty_stiffness,
                 Precision                 contact_clearance,
                 bool                      flip_master_normal)
    : master_surfaces(std::move(master)),
      slave_nodes(std::move(slave)),
      slave_surfaces(nullptr),
      penalty(penalty_stiffness),
      clearance(contact_clearance),
      flip_normal(flip_master_normal) {
    logging::error(master_surfaces != nullptr,
        "CONTACT requires a master surface set");
    logging::error(slave_nodes != nullptr,
        "CONTACT requires a slave node set");
    logging::error(penalty > Precision(0),
        "CONTACT PENALTY must be positive");
}

/**
 * @brief Constructs contact with a surface-based slave region.
 *
 * Slave nodes are collected from the selected surfaces and receive positive
 * lumped tributary areas during force and tangent assembly.
 *
 * @param master Master finite-element surfaces reconstructed as one Nagata surface.
 * @param slave Slave surfaces defining nodes and tributary weights.
 * @param penalty_stiffness Positive normal penalty stiffness.
 * @param contact_clearance Signed normal clearance subtracted from the gap.
 * @param flip_master_normal Whether to reverse reconstructed master normals.
 */
Contact::Contact(model::SurfaceRegion::Ptr master,
                 model::SurfaceRegion::Ptr slave,
                 Precision                 penalty_stiffness,
                 Precision                 contact_clearance,
                 bool                      flip_master_normal)
    : master_surfaces(std::move(master)),
      slave_nodes(nullptr),
      slave_surfaces(std::move(slave)),
      penalty(penalty_stiffness),
      clearance(contact_clearance),
      flip_normal(flip_master_normal) {
    logging::error(master_surfaces != nullptr,
        "CONTACT requires a master surface set");
    logging::error(slave_surfaces != nullptr,
        "CONTACT requires a slave surface set");
    logging::error(penalty > Precision(0),
        "CONTACT PENALTY must be positive");
}

/**
 * @brief Pushes a complete child contact state onto the transactional stack.
 *
 * Frozen trials preserve the connected master component selected by the parent
 * while allowing tracked locations to cross internal charts. If no parent
 * location exists yet, the first frozen evaluation acquires an initial chart
 * without persisting it outside that transaction. Update trials may acquire a
 * new component once and freeze it after their first assembly.
 *
 * @param freeze_components Whether current master components are immediately frozen.
 * @param freeze_after_update Whether to freeze after the first assembly update.
 */
void Contact::begin_trial(bool freeze_components, bool freeze_after_update) const {
    const AssemblyState& parent_state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    runtime_state.trials.push_back({
        parent_state,
        freeze_components,
        freeze_after_update
    });
}

/**
 * @brief Commits the current child state into its parent or persistent state.
 */
void Contact::commit_trial() const {
    logging::error(!runtime_state.trials.empty(),
        "CONTACT has no active trial state to commit");

    AssemblyState accepted_state = std::move(runtime_state.trials.back().state);
    runtime_state.trials.pop_back();

    if (runtime_state.trials.empty()) {
        runtime_state.committed = std::move(accepted_state);
    } else {
        runtime_state.trials.back().state = std::move(accepted_state);
    }
}

/**
 * @brief Discards the current child contact state without modifying its parent.
 */
void Contact::rollback_trial() const {
    logging::error(!runtime_state.trials.empty(),
        "CONTACT has no active trial state to roll back");

    runtime_state.trials.pop_back();
}

/**
 * @brief Reports a change of active slave/component pairs.
 *
 * Nagata patch and originating FE-surface changes are deliberately excluded
 * because both are internal charts of one physical connected master component.
 *
 * @return True if the most recent assembly changed the discrete contact signature.
 */
bool Contact::partner_signature_changed() const {
    const AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    return state.last_signature_changed;
}

/**
 * @brief Performs the outer augmented-Lagrange multiplier update.
 *
 * A changed active/component signature is first equilibrated with existing
 * multipliers. Otherwise each projected slave updates
 * `lambda <- max(0, lambda - penalty * gap)` outside a length-scaled geometric
 * tolerance. The transactional state owns the resulting multipliers.
 *
 * @return True if at least one multiplier changed beyond its tolerance.
 */
bool Contact::update_augmented_lagrange() const {
    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    // A discrete component or active-set change defines a new inner problem.
    // Re-equilibrate it before changing multipliers to avoid a two-cycle.
    if (state.last_signature_changed) {
        state.last_augmentation_changed = false;

        if (print_contact_summary) {
            std::cout
                << std::scientific
                << std::setprecision(3)
                << "[CONTACT_AUGMENT]"
                << " deferred=1"
                << " updated=0"
                << " max_dlambda=" << Precision(0)
                << " max_pen="     << Precision(0)
                << '\n';
        }

        return false;
    }

    Index     updated_points      = 0;
    Precision maximum_penetration = Precision(0);
    Precision maximum_change      = Precision(0);

    for (const auto& [slave_node_id, location] : state.locations) {
        (void) location;

        const auto gap_iterator = state.gaps.find(slave_node_id);

        if (gap_iterator == state.gaps.end() || !std::isfinite(gap_iterator->second)) {
            continue;
        }

        const Precision gap = gap_iterator->second;

        const auto length_iterator = state.characteristic_lengths.find(slave_node_id);

        const Precision characteristic_length =
            length_iterator != state.characteristic_lengths.end()
                ? std::max(length_iterator->second, Precision(0))
                : Precision(0);

        const Precision gap_tolerance =
            std::max(
                augmentation_gap_absolute_tolerance,
                augmentation_gap_relative_tolerance * characteristic_length);

        Precision& normal_multiplier = state.normal_multipliers[slave_node_id];

        const Precision old_multiplier = std::max(normal_multiplier, Precision(0));
        Precision       new_multiplier = old_multiplier;

        // Retain the multiplier in the geometric tolerance band. Outside it,
        // penetrations increase and openings reduce the non-negative value.
        if (gap < -gap_tolerance
            || (gap > gap_tolerance && old_multiplier > Precision(0))) {
            new_multiplier =
                std::max(Precision(0), old_multiplier - penalty * gap);
        }

        const Precision multiplier_change = std::abs(new_multiplier - old_multiplier);

        const Precision multiplier_scale =
            std::max({
                old_multiplier,
                new_multiplier,
                penalty * gap_tolerance,
                Precision(1)
            });

        const Precision multiplier_tolerance =
            std::max(
                augmentation_multiplier_absolute_tolerance,
                augmentation_multiplier_relative_tolerance * multiplier_scale);

        normal_multiplier = new_multiplier;

        maximum_penetration =
            std::max(maximum_penetration, std::max(Precision(0), -gap));
        maximum_change = std::max(maximum_change, multiplier_change);

        if (multiplier_change > multiplier_tolerance) {
            ++updated_points;
        }
    }

    state.last_augmentation_changed = updated_points > 0;

    if (print_contact_summary) {
        std::cout
            << std::scientific
            << std::setprecision(3)
            << "[CONTACT_AUGMENT]"
            << " updated="     << updated_points
            << " max_dlambda=" << maximum_change
            << " max_pen="     << maximum_penetration
            << '\n';
    }

    return state.last_augmentation_changed;
}

/**
 * @brief Reconstructs the Nagata master and assembles contact contributions.
 *
 * The algorithm performs the following phases:
 *
 * 1. Reconstruct the complete smooth master from current nodal coordinates.
 * 2. Initialize deterministic slave membership and positive slave weights.
 * 3. Continue established locations or globally acquire released slave points.
 * 4. Evaluate Nagata position, normal and signed clearance gap.
 * 5. Determine the local Nagata coordinate dependency stencil.
 * 6. Differentiate the projected gap and assemble variational contributions.
 * 7. Replace the transactional locations, components and active-set geometry.
 *
 * A tracked projection remains on its connected component and may cross any
 * internal Nagata patch or source FE-surface boundary. Global acquisition
 * returns the closest point without an additional geometric search cutoff.
 *
 * @param system_nodal_dofs Global nodal degree-of-freedom equation numbers.
 * @param model_data Current model topology and position field.
 * @param nodal_forces Global nodal residual-force accumulator.
 * @param triplets Global sparse Gauss-Newton tangent triplet accumulator.
 * @param assemble_tangent Whether to assemble the sparse Gauss-Newton block.
 */
void Contact::assemble(SystemDofIds&     system_nodal_dofs,
                       model::ModelData& model_data,
                       model::NodeData&  nodal_forces,
                       TripletList&      triplets,
                       bool              assemble_tangent) const {
    // ---------------------------------------------------------------------
    // Validate current geometry and reconstruct the smooth master
    // ---------------------------------------------------------------------

    logging::error(model_data.positions != nullptr,
        "CONTACT requires the current model position field");
    logging::error(master_surfaces != nullptr && !master_surfaces->data().empty(),
        "CONTACT requires a non-empty master surface set");

    const model::Field& node_coords = *model_data.positions;
    const auto&         surfaces    = model_data.surfaces;

    logging::error(node_coords.domain == model::FieldDomain::NODE
                && node_coords.components >= 3,
        "CONTACT requires a nodal coordinate field with three components");

    const auto assembly_start       = std::chrono::steady_clock::now();
    const std::size_t initial_triplets = triplets.size();

    model::NagataSurface master_geometry(*master_surfaces, surfaces, node_coords);

    // Map retained source pointers back to stable model surface IDs for
    // diagnostics and self-contact validation without exposing IDs in Nagata.
    std::unordered_map<const model::SurfaceInterface*, ID> source_surface_ids;
    source_surface_ids.reserve(static_cast<std::size_t>(master_surfaces->size()));

    for (const ID surface_id : *master_surfaces) {
        source_surface_ids.emplace(
            surfaces[static_cast<std::size_t>(surface_id)].get(),
            surface_id);
    }

    // ---------------------------------------------------------------------
    // Select transactional state and initialize fixed slave information
    // ---------------------------------------------------------------------

    ++runtime_state.call;

    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    const bool freeze_master_components =
        !runtime_state.trials.empty()
        && runtime_state.trials.back().freeze_master_components;

    if (!runtime_state.slave_nodes_initialized) {
        runtime_state.slave_node_ids =
            collect_slave_nodes(slave_nodes, slave_surfaces, model_data);
        runtime_state.slave_nodes_initialized = true;
    }

    if (!runtime_state.slave_weights_initialized) {
        runtime_state.slave_tributary_areas =
            build_slave_tributary_areas(slave_surfaces, model_data, node_coords);
        runtime_state.slave_weights_initialized = true;
    }

    const std::vector<ID>& slave_node_ids = runtime_state.slave_node_ids;
    const auto& slave_tributary_areas = runtime_state.slave_tributary_areas;

    ContactDiagnostics diagnostics;
    diagnostics.slave_nodes = static_cast<Index>(slave_node_ids.size());

    std::unordered_map<ID, model::NagataSurface::Location> current_locations;
    std::unordered_map<ID, model::nagata::ComponentID>      current_components;
    std::unordered_set<ID>                                 current_active_slaves;
    std::unordered_map<ID, Precision>                      current_multipliers;
    std::unordered_map<ID, Precision>                      current_gaps;
    std::unordered_map<ID, Precision>                      current_characteristic_lengths;

    current_locations.reserve(slave_node_ids.size());
    current_components.reserve(slave_node_ids.size());
    current_active_slaves.reserve(slave_node_ids.size());
    current_multipliers.reserve(slave_node_ids.size());
    current_gaps.reserve(slave_node_ids.size());
    current_characteristic_lengths.reserve(slave_node_ids.size());

    // ---------------------------------------------------------------------
    // Project every slave and assemble active contact pairs
    // ---------------------------------------------------------------------

    for (const ID slave_node_id : slave_node_ids) {
        logging::error(slave_node_id >= 0
                    && static_cast<Index>(slave_node_id) < node_coords.rows,
            "CONTACT slave node ID ", slave_node_id,
            " lies outside the coordinate field");

        const Vec3 slave_position =
            node_coords.row_vec3(static_cast<Index>(slave_node_id));

        logging::error(slave_position.allFinite(),
            "CONTACT slave node ID ", slave_node_id,
            " has invalid current coordinates");

        const auto previous_location_iterator = state.locations.find(slave_node_id);
        const auto previous_component_iterator = state.components.find(slave_node_id);
        const auto previous_multiplier_iterator = state.normal_multipliers.find(slave_node_id);

        const Precision normal_multiplier =
            previous_multiplier_iterator != state.normal_multipliers.end()
                ? std::max(previous_multiplier_iterator->second, Precision(0))
                : Precision(0);

        const bool previously_active =
            state.active_slaves.find(slave_node_id) != state.active_slaves.end();

        const bool has_previous_location =
            previous_location_iterator != state.locations.end()
            && master_geometry.valid(previous_location_iterator->second);

        const bool track_established_location =
            has_previous_location
            && (freeze_master_components || previously_active
                || normal_multiplier > Precision(0));

        model::NagataSurface::Location location;
        bool                           has_location = false;

        // Established locations are continued on their connected component.
        // A frozen trial without inherited geometry must still acquire an
        // initial chart. This occurs for the state-neutral predictor before the
        // outer update trial has performed its first assembly.
        if (track_established_location) {
            location = master_geometry.project(
                slave_position,
                previous_location_iterator->second);
            has_location = true;
            ++diagnostics.tracked_projections;
        } else if (!freeze_master_components || !has_previous_location) {
            location = master_geometry.project(slave_position);
            has_location = true;
            ++diagnostics.global_projections;
        }

        if (!has_location) {
            ++diagnostics.no_partner;

            if (previously_active || normal_multiplier > Precision(0)) {
                ++diagnostics.active_partner_losses;
            }

            continue;
        }

        model::NagataSurface::Evaluation evaluation =
            master_geometry.evaluate(location);

        logging::error(evaluation.surface != nullptr,
            "CONTACT Nagata evaluation has no originating FE surface");
        logging::error(evaluation.position.allFinite()
                    && evaluation.normal.allFinite(),
            "CONTACT Nagata evaluation produced invalid geometry");

        const auto source_id_iterator = source_surface_ids.find(evaluation.surface);

        logging::error(source_id_iterator != source_surface_ids.end(),
            "CONTACT Nagata evaluation references an unknown source surface");

        // A shared slave/master node must not contact the FE surface containing
        // itself. Global Nagata projection currently exposes no second-best
        // source candidate, so the rejected point remains without a partner.
        if (surface_contains_node(*evaluation.surface, slave_node_id)) {
            ++diagnostics.self_contact_rejections;
            ++diagnostics.no_partner;

            if (previously_active || normal_multiplier > Precision(0)) {
                ++diagnostics.active_partner_losses;
            }

            continue;
        }

        Vec3 normal = evaluation.normal;

        if (flip_normal) {
            normal = -normal;
        }

        const Vec3      difference = slave_position - evaluation.position;
        const Precision closest_distance = difference.norm();
        const Precision gap = difference.dot(normal) - clearance;

        logging::error(std::isfinite(closest_distance) && std::isfinite(gap),
            "CONTACT produced an invalid Nagata distance or normal gap");

        const model::nagata::ComponentID component =
            master_geometry.component(location);

        const auto slave_area_iterator = slave_tributary_areas.find(slave_node_id);

        const Precision slave_weight =
            slave_area_iterator != slave_tributary_areas.end()
                ? slave_area_iterator->second
                : Precision(1);

        logging::error(std::isfinite(slave_weight) && slave_weight > Precision(0),
            "CONTACT slave node ID ", slave_node_id,
            " has invalid tributary weight ", slave_weight);

        const Precision characteristic_length =
            std::sqrt(std::max(evaluation.surface->area(node_coords), Precision(0)));

        logging::error(std::isfinite(characteristic_length)
                    && characteristic_length > Precision(0),
            "CONTACT source surface has an invalid characteristic length");

        // Component and active-set decisions remain fixed inside Newton and
        // line-search trials. The outer update performs the discrete contact
        // classification and requests another Newton solve if it changes.
        const ContactLaw contact_law =
            freeze_master_components
                ? continue_augmented_lagrange_branch(
                    gap, normal_multiplier, penalty, previously_active)
                : evaluate_augmented_lagrange_law(
                    gap, normal_multiplier, penalty);

        current_locations[slave_node_id]              = location;
        current_components[slave_node_id]             = component;
        current_multipliers[slave_node_id]            = normal_multiplier;
        current_gaps[slave_node_id]                   = gap;
        current_characteristic_lengths[slave_node_id] = characteristic_length;

        diagnostics.maximum_closest_distance =
            std::max(diagnostics.maximum_closest_distance, closest_distance);

        if (previous_component_iterator != state.components.end()
            && previous_component_iterator->second != component) {
            ++diagnostics.component_switches;
        }

        if (!contact_law.active) {
            ++diagnostics.open_closest_partner;
            continue;
        }

        ++diagnostics.active_contacts;
        current_active_slaves.insert(slave_node_id);

        if (!previously_active) {
            ++diagnostics.activations;
        }

        const Precision penetration = std::max(Precision(0), -gap);
        const Precision slave_force = slave_weight * contact_law.pressure;

        diagnostics.slave_force_squared_sum += slave_force * slave_force;
        diagnostics.maximum_slave_force =
            std::max(diagnostics.maximum_slave_force, slave_force);

        if (penetration > diagnostics.maximum_penetration) {
            diagnostics.maximum_penetration = penetration;
            diagnostics.worst_slave         = slave_node_id;
            diagnostics.worst_component     = component;
            diagnostics.worst_patch         = location.patch;
            diagnostics.worst_local         = location.local;
            diagnostics.worst_gap           = gap;
            diagnostics.worst_distance      = closest_distance;
        }

        // Only the physical connected component participates in the discrete
        // signature. Internal Nagata and FE chart transitions remain smooth.
        hash_value(
            diagnostics.signature,
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(slave_node_id)));
        hash_value(
            diagnostics.signature,
            static_cast<std::uint64_t>(component));

        assemble_contact_pair(
            node_coords,
            master_geometry,
            location,
            normal,
            slave_node_id,
            contact_law,
            penalty,
            slave_weight,
            characteristic_length,
            system_nodal_dofs,
            nodal_forces,
            triplets,
            assemble_tangent,
            diagnostics);
    }

    // ---------------------------------------------------------------------
    // Finalize transactional state and discrete change detection
    // ---------------------------------------------------------------------

    for (const ID slave_node_id : state.active_slaves) {
        if (current_active_slaves.find(slave_node_id) == current_active_slaves.end()) {
            ++diagnostics.deactivations;
        }
    }

    const bool signature_changed =
        state.previous_signature != 0
        && state.previous_signature != diagnostics.signature;

    const std::int64_t active_change =
        static_cast<std::int64_t>(diagnostics.active_contacts)
        - static_cast<std::int64_t>(state.previous_active);

    state.locations                 = std::move(current_locations);
    state.components                = std::move(current_components);
    state.active_slaves             = std::move(current_active_slaves);
    state.normal_multipliers        = std::move(current_multipliers);
    state.gaps                      = std::move(current_gaps);
    state.characteristic_lengths    = std::move(current_characteristic_lengths);
    state.previous_signature        = diagnostics.signature;
    state.previous_active           = diagnostics.active_contacts;
    state.last_signature_changed    = signature_changed;

    if (!freeze_master_components
        && !runtime_state.trials.empty()
        && runtime_state.trials.back().freeze_after_update) {
        runtime_state.trials.back().freeze_master_components = true;
        runtime_state.trials.back().freeze_after_update      = false;
    }

    // ---------------------------------------------------------------------
    // Report compact geometry and assembly diagnostics
    // ---------------------------------------------------------------------

    const Precision contact_force_norm =
        std::sqrt(diagnostics.slave_force_squared_sum);

    const auto assembly_end = std::chrono::steady_clock::now();
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(assembly_end - assembly_start).count();

    const std::size_t added_triplets = triplets.size() - initial_triplets;

    if (print_contact_summary) {
        const auto previous_flags     = std::cout.flags();
        const auto previous_precision = std::cout.precision();

        std::cout
            << std::scientific
            << std::setprecision(3)
            << "[CONTACT]"
            << " call="             << runtime_state.call
            << " depth="            << runtime_state.trials.size()
            << " frozen="           << (freeze_master_components ? 1 : 0)
            << " active="           << diagnostics.active_contacts
            << " d_active="         << active_change
            << " activated="        << diagnostics.activations
            << " deactivated="      << diagnostics.deactivations
            << " component_switch=" << diagnostics.component_switches
            << " changed="          << (signature_changed ? 1 : 0)
            << " lost_active="      << diagnostics.active_partner_losses
            << " self_rej="         << diagnostics.self_contact_rejections
            << " global_proj="      << diagnostics.global_projections
            << " tracked_proj="     << diagnostics.tracked_projections
            << " no_partner="       << diagnostics.no_partner
            << " open="             << diagnostics.open_closest_partner
            << " K_approx="         << diagnostics.approximate_tangents
            << " max_dist="         << diagnostics.maximum_closest_distance
            << " max_pen="          << diagnostics.maximum_penetration
            << " max_force="        << diagnostics.maximum_slave_force
            << " force_norm="       << contact_force_norm
            << " max_Kn="           << diagnostics.maximum_tangent_norm
            << " triplets="         << added_triplets
            << " ms="               << elapsed_ms
            << " signature="        << diagnostics.signature
            << '\n';

        if (diagnostics.active_contacts > 0) {
            std::cout
                << std::scientific
                << std::setprecision(6)
                << "[CONTACT_MAX]"
                << " call="      << runtime_state.call
                << " slave="     << diagnostics.worst_slave
                << " component=" << diagnostics.worst_component
                << " patch="     << diagnostics.worst_patch
                << " local=("    << diagnostics.worst_local(0)
                << ','            << diagnostics.worst_local(1)
                << ")"
                << " gap="       << diagnostics.worst_gap
                << " pen="       << diagnostics.maximum_penetration
                << " distance="  << diagnostics.worst_distance
                << '\n';
        }

        std::cout.flags(previous_flags);
        std::cout.precision(previous_precision);
    }
}

} // namespace fem::constraint
