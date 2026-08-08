/**
 * @file nagata_surface.cpp
 * @brief Implements preprocessing of connected G1 Nagata surfaces.
 *
 * Construction resolves the finite-element surfaces referenced by one model
 * surface region, merges their nodes into smooth reconstructed vertices,
 * averages incident normals, triangulates supported FE natural domains, builds
 * patch adjacency and computes the quadratic and rational Nagata coefficients.
 *
 * @see NagataSurface
 * @see nagata_surface_evaluate.cpp
 * @see nagata_surface_project.cpp
 *
 * @author Finn Eggers
 * @date 08.08.2026
 */

#include "nagata_surface.h"

#include "../../../core/logging.h"

#include <Eigen/QR>

#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include <utility>

namespace fem::model {

/**
 * @brief Reconstructs all finite-element surfaces referenced by a surface set.
 *
 * Surface IDs are first resolved against the supplied container and retained as
 * shared pointers. Each global FE node then becomes one smooth reconstructed
 * vertex whose normal is the normalized sum of the unit normals evaluated by
 * all incident selected surfaces.
 *
 * Supported natural domains are subdivided into one triangle for S3, two for
 * S4, four for S6 and six for S8. Every patch stores affine source-natural
 * coordinates, explicit shared-edge adjacency, the quadratic Nagata edge
 * coefficients and the rational G1 correction coefficients. A point BVH over
 * reconstructed vertices supplies seeds for later global projection.
 *
 * @param surface_set Surface IDs included in the reconstruction.
 * @param surfaces Global surface container indexed by `surface_set` IDs.
 * @param node_coords Finite global nodal coordinate field for this snapshot.
 */
NagataSurface::NagataSurface(
    const SurfaceRegion&                      surface_set,
    const std::vector<SurfaceInterface::Ptr>& surfaces,
    const Field&                              node_coords) {
    // ---------------------------------------------------------------------
    // Resolve and validate the selected source surfaces
    // ---------------------------------------------------------------------

    logging::error(!surface_set.data().empty(),
        "NagataSurface requires a non-empty surface set");
    logging::error(node_coords.domain == FieldDomain::NODE,
        "NagataSurface requires a NODE coordinate field");
    logging::error(node_coords.components >= 3,
        "NagataSurface requires at least three coordinate components");

    source_surfaces_.reserve(surface_set.data().size());

    for (const ID surface_id : surface_set) {
        logging::error(surface_id >= 0
                    && static_cast<std::size_t>(surface_id) < surfaces.size(),
            "NagataSurface surface ID ", surface_id,
            " lies outside the supplied surface container");

        const SurfaceInterface::Ptr& surface =
            surfaces[static_cast<std::size_t>(surface_id)];

        logging::error(surface != nullptr,
            "NagataSurface surface ID ", surface_id, " is null");
        logging::error(surface->n_nodes == 3
                    || surface->n_nodes == 4
                    || surface->n_nodes == 6
                    || surface->n_nodes == 8,
            "NagataSurface supports only 3-, 4-, 6- and 8-node surfaces");

        source_surfaces_.push_back(surface);
    }

    constexpr Precision normal_tolerance = Precision(1e-12);

    // Nearly parallel normal constraints represent one geometrically resolved
    // direction. Truncating their small singular mode prevents infinitesimal
    // coordinate changes from producing finite jumps in Nagata coefficients.
    constexpr Precision coefficient_rank_tolerance = Precision(1e-4);

    std::unordered_map<ID, Index> vertex_by_node;

    // ---------------------------------------------------------------------
    // Gather reconstructed vertices and accumulate incident FE normals
    // ---------------------------------------------------------------------

    for (Index surface_index = 0;
         surface_index < static_cast<Index>(source_surfaces_.size());
         ++surface_index) {
        const SurfaceInterface::Ptr& surface =
            source_surfaces_[static_cast<std::size_t>(surface_index)];
        const DynamicMatrix natural_coords = surface->node_coords_natural();

        logging::error(static_cast<Index>(natural_coords.rows()) == surface->n_nodes
                    && natural_coords.cols() == 2,
            "NagataSurface source natural coordinates do not match surface topology");

        for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
            const ID node_id = surface->nodes()[local_node];

            logging::error(node_id >= 0
                        && static_cast<Index>(node_id) < node_coords.rows,
                "NagataSurface source node ", node_id,
                " lies outside the supplied coordinate field");

            auto iterator = vertex_by_node.find(node_id);

            if (iterator == vertex_by_node.end()) {
                const Index vertex_index = static_cast<Index>(vertices_.size());

                Vertex vertex;
                vertex.node_id  = node_id;
                vertex.position = node_coords.row_vec3(static_cast<Index>(node_id));

                logging::error(vertex.position.allFinite(),
                    "NagataSurface source node ", node_id,
                    " has invalid coordinates");

                vertices_.push_back(vertex);
                iterator = vertex_by_node.emplace(node_id, vertex_index).first;
            }

            const Vec2 local =
                natural_coords.row(static_cast<Eigen::Index>(local_node)).transpose();

            Vec3 normal = surface->normal(node_coords, local);

            logging::error(normal.allFinite() && normal.norm() > normal_tolerance,
                "NagataSurface encountered an invalid normal at source node ", node_id);

            normal.normalize();

            // All incident surface normals contribute equally to the current
            // smooth-only vertex normal. Hard-edge grouping is intentionally
            // absent from this MVP.
            Vertex& vertex = vertices_[static_cast<std::size_t>(iterator->second)];
            vertex.normal += normal;
            vertex.source_surfaces.push_back(surface_index);
            vertex.source_local_nodes.push_back(local_node);
        }
    }

    // Normalize the accumulated smooth normal at every reconstructed vertex.
    for (Vertex& vertex : vertices_) {
        logging::error(vertex.normal.allFinite() && vertex.normal.norm() > normal_tolerance,
            "NagataSurface cannot construct a unique averaged normal at node ",
            vertex.node_id);

        vertex.normal.normalize();
    }

    // ---------------------------------------------------------------------
    // Triangulate each selected FE surface in its natural domain
    // ---------------------------------------------------------------------

    for (Index surface_index = 0;
         surface_index < static_cast<Index>(source_surfaces_.size());
         ++surface_index) {
        const SurfaceInterface::Ptr& surface =
            source_surfaces_[static_cast<std::size_t>(surface_index)];

        const DynamicMatrix natural_coords = surface->node_coords_natural();
        std::vector<Index> surface_vertices(static_cast<std::size_t>(surface->n_nodes));

        for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
            const ID node_id = surface->nodes()[local_node];
            const auto iterator = vertex_by_node.find(node_id);

            logging::error(iterator != vertex_by_node.end(),
                "NagataSurface internal vertex lookup failed for node ", node_id);

            surface_vertices[static_cast<std::size_t>(local_node)] = iterator->second;
        }

        // Append one counter-clockwise triangular patch and retain the affine
        // mapping from patch barycentric coordinates to source FE coordinates.
        const auto append_patch = [&](Index a, Index b, Index c) {
            Patch patch;

            patch.surface = surface_index;
            patch.vertices = {
                surface_vertices[static_cast<std::size_t>(a)],
                surface_vertices[static_cast<std::size_t>(b)],
                surface_vertices[static_cast<std::size_t>(c)]
            };

            patch.element_local = {
                natural_coords.row(static_cast<Eigen::Index>(a)).transpose(),
                natural_coords.row(static_cast<Eigen::Index>(b)).transpose(),
                natural_coords.row(static_cast<Eigen::Index>(c)).transpose()
            };

            patches_.push_back(patch);
        };

        if (surface->n_nodes == 3) {
            append_patch(0, 1, 2);
        } else if (surface->n_nodes == 4) {
            append_patch(0, 1, 2);
            append_patch(0, 2, 3);
        } else if (surface->n_nodes == 6) {
            // Corner-first S6 ordering followed by edge nodes 0-1, 1-2, 2-0.
            append_patch(0, 3, 5);
            append_patch(3, 1, 4);
            append_patch(5, 4, 2);
            append_patch(3, 4, 5);
        } else {
            // Corner-first S8 ordering followed by edge nodes 0-1, 1-2, 2-3, 3-0.
            append_patch(0, 4, 7);
            append_patch(1, 5, 4);
            append_patch(2, 6, 5);
            append_patch(3, 7, 6);
            append_patch(4, 5, 6);
            append_patch(4, 6, 7);
        }
    }

    logging::error(!patches_.empty(),
        "NagataSurface construction produced no triangular patches");

    // ---------------------------------------------------------------------
    // Build vertex incidence and shared-edge patch adjacency
    // ---------------------------------------------------------------------

    for (nagata::PatchID patch_id = 0;
         patch_id < static_cast<nagata::PatchID>(patches_.size());
         ++patch_id) {
        const Patch& patch = patches_[static_cast<std::size_t>(patch_id)];

        for (Index local_vertex = 0; local_vertex < 3; ++local_vertex) {
            vertices_[static_cast<std::size_t>(
                patch.vertices[static_cast<std::size_t>(local_vertex)])]
                .patches.push_back(patch_id);
        }
    }

    using EdgeOccurrence = std::pair<nagata::PatchID, Index>;
    using EdgeKey        = std::pair<Index, Index>;

    std::map<EdgeKey, std::vector<EdgeOccurrence>> edge_occurrences;

    for (nagata::PatchID patch_id = 0;
         patch_id < static_cast<nagata::PatchID>(patches_.size());
         ++patch_id) {
        const Patch& patch = patches_[static_cast<std::size_t>(patch_id)];

        for (Index edge = 0; edge < 3; ++edge) {
            const Index j = (edge + 1) % 3;
            const Index k = (edge + 2) % 3;

            const Index vertex_j = patch.vertices[static_cast<std::size_t>(j)];
            const Index vertex_k = patch.vertices[static_cast<std::size_t>(k)];
            const EdgeKey key = vertex_j < vertex_k
                ? EdgeKey{vertex_j, vertex_k}
                : EdgeKey{vertex_k, vertex_j};

            edge_occurrences[key].push_back({patch_id, edge});
        }
    }

    for (const auto& entry : edge_occurrences) {
        const std::vector<EdgeOccurrence>& occurrences = entry.second;

        logging::error(occurrences.size() <= 2,
            "NagataSurface does not support non-manifold patch edges");

        if (occurrences.size() == 1) {
            continue;
        }

        const nagata::PatchID patch_a = occurrences[0].first;
        const Index            edge_a  = occurrences[0].second;
        const nagata::PatchID patch_b = occurrences[1].first;
        const Index            edge_b  = occurrences[1].second;

        patches_[static_cast<std::size_t>(patch_a)]
            .neighbors[static_cast<std::size_t>(edge_a)] = {patch_b, edge_b};
        patches_[static_cast<std::size_t>(patch_b)]
            .neighbors[static_cast<std::size_t>(edge_b)] = {patch_a, edge_a};
    }

    // ---------------------------------------------------------------------
    // Label connected components independently of internal patch charts
    // ---------------------------------------------------------------------

    nagata::ComponentID component_id = 0;

    for (nagata::PatchID seed = 0;
         seed < static_cast<nagata::PatchID>(patches_.size());
         ++seed) {
        Patch& seed_patch = patches_[static_cast<std::size_t>(seed)];

        if (seed_patch.component != std::numeric_limits<nagata::ComponentID>::max()) {
            continue;
        }

        std::vector<nagata::PatchID> pending{seed};
        seed_patch.component = component_id;

        while (!pending.empty()) {
            const nagata::PatchID patch_id = pending.back();
            pending.pop_back();

            const Patch& patch = patches_[static_cast<std::size_t>(patch_id)];

            for (const Neighbor& neighbor : patch.neighbors) {
                if (neighbor.patch == invalid_patch) {
                    continue;
                }

                Patch& next = patches_[static_cast<std::size_t>(neighbor.patch)];

                if (next.component != std::numeric_limits<nagata::ComponentID>::max()) {
                    continue;
                }

                next.component = component_id;
                pending.push_back(neighbor.patch);
            }
        }

        ++component_id;
    }

    // ---------------------------------------------------------------------
    // Construct quadratic C0 and rational G1 Nagata coefficients
    // ---------------------------------------------------------------------

    for (nagata::PatchID patch_id = 0;
         patch_id < static_cast<nagata::PatchID>(patches_.size());
         ++patch_id) {
        Patch& patch = patches_[static_cast<std::size_t>(patch_id)];

        std::array<Vec3, 3> position;
        std::array<Vec3, 3> normal;

        for (Index i = 0; i < 3; ++i) {
            const Vertex& vertex = vertices_[static_cast<std::size_t>(
                patch.vertices[static_cast<std::size_t>(i)])];

            position[static_cast<std::size_t>(i)] = vertex.position;
            normal  [static_cast<std::size_t>(i)] = vertex.normal;
        }

        // Edge i joins vertices j and k. Its curvature vector is the
        // minimum-norm solution of both endpoint-normal constraints.
        for (Index i = 0; i < 3; ++i) {
            const Index j = (i + 1) % 3;
            const Index k = (i + 2) % 3;
            const Vec3 edge = position[static_cast<std::size_t>(k)]
                            - position[static_cast<std::size_t>(j)];

            logging::error(edge.norm() > normal_tolerance,
                "NagataSurface patch ", patch_id, " contains a zero-length edge");

            StaticMatrix<2, 3> system;
            system.row(0) = normal[static_cast<std::size_t>(j)].transpose();
            system.row(1) = normal[static_cast<std::size_t>(k)].transpose();

            Vec2 rhs;
            rhs(0) =  normal[static_cast<std::size_t>(j)].dot(edge);
            rhs(1) = -normal[static_cast<std::size_t>(k)].dot(edge);

            Eigen::CompleteOrthogonalDecomposition<StaticMatrix<2, 3>> decomposition(system);
            decomposition.setThreshold(coefficient_rank_tolerance);

            const Vec3 curvature = decomposition.solve(rhs);

            logging::error(curvature.allFinite(),
                "NagataSurface failed to construct C0 coefficient for patch ",
                patch_id, " edge ", i);

            patch.q[static_cast<std::size_t>(i)] =
                position[static_cast<std::size_t>(j)]
                + position[static_cast<std::size_t>(k)]
                - curvature;
        }

        for (Index i = 0; i < 3; ++i) {
            const Index j = (i + 1) % 3;
            const Index k = (i + 2) % 3;

            const Vec3& pj = position[static_cast<std::size_t>(j)];
            const Vec3& pk = position[static_cast<std::size_t>(k)];
            const Vec3& nj = normal[static_cast<std::size_t>(j)];
            const Vec3& nk = normal[static_cast<std::size_t>(k)];
            const Vec3& qi = patch.q[static_cast<std::size_t>(i)];
            const Vec3& qj = patch.q[static_cast<std::size_t>(j)];
            const Vec3& qk = patch.q[static_cast<std::size_t>(k)];

            // Nagata edge tangent and cross-boundary vectors determine the
            // normalized edge-normal coefficients and their C0 mismatch.
            const Vec3 t_ij = qi - Precision(2) * pj;
            const Vec3 t_ik = qi - Precision(2) * pk;
            const Vec3 s_ij = qi - qk;
            const Vec3 s_ik = Precision(2) * pk - qj;

            const Precision c = nk.dot(t_ij);
            const Precision d = nj.dot(t_ik);

            Precision c_hat = Precision(1);
            Precision d_hat = Precision(1);

            if (std::abs(c) + std::abs(d) > normal_tolerance) {
                const Precision scale = Precision(1) / (std::abs(c) + std::abs(d));
                c_hat = scale * c;
                d_hat = scale * d;
            }

            const Precision edge_mismatch =
                c_hat * nj.dot(s_ik) + d_hat * nk.dot(s_ij);

            StaticMatrix<2, 3> system;
            system.row(0) = c_hat * nj.transpose();
            system.row(1) = d_hat * nk.transpose();

            const Vec2 rhs = Vec2::Constant(edge_mismatch);
            Vec3& gamma = patch.gamma[static_cast<std::size_t>(i)];

            Eigen::CompleteOrthogonalDecomposition<StaticMatrix<2, 3>> decomposition(system);
            decomposition.setThreshold(coefficient_rank_tolerance);

            gamma = decomposition.solve(rhs);

            logging::error(gamma.allFinite(),
                "NagataSurface failed to construct G1 coefficient for patch ",
                patch_id, " edge ", i);
        }
    }

    // ---------------------------------------------------------------------
    // Build the nearest-vertex acceleration structure
    // ---------------------------------------------------------------------

    logging::error(vertices_.size()
                   <= static_cast<std::size_t>(std::numeric_limits<ID>::max()),
        "NagataSurface contains too many vertices for BvhAabb IDs");

    for (Index vertex = 0; vertex < static_cast<Index>(vertices_.size()); ++vertex) {
        vertex_bvh_.add_point(
            static_cast<ID>(vertex),
            vertices_[static_cast<std::size_t>(vertex)].position);
    }

    vertex_bvh_.finalize();

    logging::error(vertex_bvh_.valid(),
        "NagataSurface failed to construct the vertex BVH");
}

/**
 * @brief Returns all global nodes controlling one Nagata patch.
 *
 * A patch depends directly on its three reconstructed vertex positions. Each
 * vertex normal is additionally the normalized sum of normals from all selected
 * source surfaces incident to that vertex. The complete dependency set is
 * therefore the union of the nodes of those incident source surfaces.
 *
 * The returned IDs are sorted and unique. They form the local coordinate
 * stencil required by contact sensitivity calculations without exposing the
 * private Nagata patch representation. The stencil is constructed lazily and
 * cached on first use so reconstructions used only for perturbed evaluation do
 * not preprocess unused dependency information.
 *
 * @param location Valid location identifying the patch of interest.
 * @return Sorted global node IDs controlling patch geometry and vertex normals.
 */
std::vector<ID> NagataSurface::dependency_nodes(const Location& location) const {
    logging::error(valid(location),
        "NagataSurface::dependency_nodes received an invalid location");

    const Patch& patch = patches_[static_cast<std::size_t>(location.patch)];

    if (!patch.dependency_nodes.empty()) {
        return patch.dependency_nodes;
    }

    std::array<ID, 3> patch_node_ids{};

    for (Index local_vertex = 0; local_vertex < 3; ++local_vertex) {
        patch_node_ids[static_cast<std::size_t>(local_vertex)] =
            vertices_[static_cast<std::size_t>(
                patch.vertices[static_cast<std::size_t>(local_vertex)])].node_id;
    }

    // A vertex normal depends on every coordinate of every selected FE surface
    // incident to that vertex. Merge those connectivities into one patch-local
    // coordinate stencil.
    for (const SurfaceInterface::Ptr& surface : source_surfaces_) {
        bool incident = false;

        for (Index local_node = 0; local_node < surface->n_nodes && !incident; ++local_node) {
            incident = std::find(
                patch_node_ids.begin(),
                patch_node_ids.end(),
                surface->nodes()[local_node]) != patch_node_ids.end();
        }

        if (!incident) {
            continue;
        }

        for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
            patch.dependency_nodes.push_back(surface->nodes()[local_node]);
        }
    }

    std::sort(patch.dependency_nodes.begin(), patch.dependency_nodes.end());
    patch.dependency_nodes.erase(
        std::unique(patch.dependency_nodes.begin(), patch.dependency_nodes.end()),
        patch.dependency_nodes.end());

    logging::error(!patch.dependency_nodes.empty(),
        "NagataSurface produced an empty patch dependency stencil");

    return patch.dependency_nodes;
}

} // namespace fem::model
