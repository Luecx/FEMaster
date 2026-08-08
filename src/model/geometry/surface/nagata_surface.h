/**
 * @file nagata_surface.h
 * @brief Defines a connected G1 Nagata reconstruction of a finite-element surface set.
 *
 * `NagataSurface` converts the surfaces referenced by a model `SurfaceRegion`
 * into connected triangular Nagata patches. The reconstructed geometry is a
 * snapshot of one nodal coordinate field; model ownership, field updates and
 * contact state remain responsibilities of the caller.
 *
 * Shared global nodes are treated as smooth. Their source-surface normals are
 * averaged before the patch coefficients are constructed. The class retains the
 * source surfaces so every reconstructed location can be mapped back to its
 * originating finite-element surface and natural coordinates.
 *
 * @see NagataSurface::Location
 * @see NagataSurface::Evaluation
 * @see BvhAabb
 *
 * @author Finn Eggers
 * @date 08.08.2026
 */

#pragma once

#include "../bvh/bvh_aabb.h"
#include "../../../data/region.h"
#include "surface_interface.h"

#include <array>
#include <limits>
#include <vector>

namespace fem::model {

namespace nagata {

using PatchID     = Index;
using ComponentID = Index;

} // namespace nagata

/**
 * @brief Connected G1 Nagata surface reconstructed from a model surface region.
 *
 * The surface region supplies stable IDs into the caller-owned surface
 * container. Supported 3-, 4-, 6- and 8-node finite-element surfaces are
 * triangulated in their natural domains. Generated triangles share reconstructed
 * vertices, averaged vertex normals and explicit edge adjacency.
 *
 * A `Location` addresses one triangular patch in reduced barycentric
 * coordinates. `evaluate()` returns its differential geometry and maps the
 * location back to the source finite-element surface. `project()` performs
 * closest-point Newton iteration and can walk across internal patch edges.
 *
 * The object owns shared pointers to the selected source surfaces, but it does
 * not own or retain the coordinate field. Geometry changes therefore require a
 * new reconstruction. Every shared node is currently smooth; hard edges require
 * distinct reconstructed vertices and are outside this MVP.
 */
class NagataSurface {
public:
    /**
     * @brief Parametric location on the complete reconstructed surface.
     *
     * `local = (xi, eta)` represents barycentric coordinates
     * `beta_0 = 1 - xi - eta`, `beta_1 = xi`, and `beta_2 = eta` on `patch`.
     */
    struct Location {
        nagata::PatchID patch = std::numeric_limits<nagata::PatchID>::max();
        Vec2             local = Vec2::Zero();
    };

    /**
     * @brief Position, derivatives and source-surface mapping at one location.
     *
     * Jacobian columns are derivatives with respect to `(xi, eta)`. The three
     * second derivatives use the same coordinates. `surface` remains valid for
     * the lifetime of this reconstruction because the source shared pointers are
     * retained internally.
     */
    struct Evaluation {
        // Reconstructed differential geometry
        Vec3               position = Vec3::Zero();
        Vec3               normal   = Vec3::Zero();
        StaticMatrix<3, 2> jacobian = StaticMatrix<3, 2>::Zero();

        Vec3 d2_xixi   = Vec3::Zero();
        Vec3 d2_xieta  = Vec3::Zero();
        Vec3 d2_etaeta = Vec3::Zero();

        // Originating finite-element surface and its natural coordinates
        const SurfaceInterface* surface       = nullptr;
        Vec2                    element_local = Vec2::Zero();
    };

    // Reconstruction from the surfaces referenced by one model surface set
    NagataSurface(const SurfaceRegion&                      surface_set,
                  const std::vector<SurfaceInterface::Ptr>& surfaces,
                  const Field&                              node_coords);

    // Global projection and tracked continuation from an existing location
    [[nodiscard]] Location project(const Vec3& point) const;
    [[nodiscard]] Location project(const Vec3& point, const Location& initial_guess) const;

    // Differential geometry and mapping to the originating FE surface
    [[nodiscard]] Evaluation evaluate(const Location& location) const;

    // Stable topology queries for tracked locations and connected components
    [[nodiscard]] bool valid(const Location& location) const;
    [[nodiscard]] nagata::ComponentID component(const Location& location) const;

private:
    static constexpr nagata::PatchID invalid_patch =
        std::numeric_limits<nagata::PatchID>::max();

    /**
     * @brief Smooth reconstructed vertex associated with one global model node.
     *
     * Position is copied from the construction field. Normal is the normalized
     * sum of unit normals contributed by all selected source surfaces incident
     * to the node. The patch list supplies global projection seeds.
     */
    struct Vertex {
        ID   node_id  = ID(-1);
        Vec3 position = Vec3::Zero();
        Vec3 normal   = Vec3::Zero();

        std::vector<nagata::PatchID> patches;
    };

    /**
     * @brief Topological connection across one triangular patch edge.
     *
     * Edge `i` is opposite local patch vertex `i`. A missing neighbor denotes a
     * true boundary of the reconstructed surface.
     */
    struct Neighbor {
        nagata::PatchID patch = invalid_patch;
        Index            edge  = 0;
    };

    /**
     * @brief Triangular rational G1 Nagata patch and source mapping.
     *
     * `vertices` and `neighbors` define reconstructed topology. `surface` and
     * `element_local` preserve the affine mapping into the originating FE
     * natural domain. `q` contains the mixed coefficients of the quadratic C0
     * patch; `gamma` contains the vector coefficients of the rational quintic
     * G1 correction.
     */
    struct Patch {
        // Reconstructed topology
        std::array<Index, 3>    vertices{};
        std::array<Neighbor, 3> neighbors{};

        // Originating finite-element surface
        Index                surface = 0;
        nagata::ComponentID component = std::numeric_limits<nagata::ComponentID>::max();
        std::array<Vec2, 3> element_local{};

        // Precomputed C0 and rational G1 coefficients
        std::array<Vec3, 3> q{
            Vec3::Zero(),
            Vec3::Zero(),
            Vec3::Zero()
        };

        std::array<Vec3, 3> gamma{
            Vec3::Zero(),
            Vec3::Zero(),
            Vec3::Zero()
        };
    };

    // Selected source FE surfaces retained for evaluation mapping and lifetime
    std::vector<SurfaceInterface::Ptr> source_surfaces_;

    // Reconstructed smooth topology and patch coefficients
    std::vector<Vertex> vertices_;
    std::vector<Patch>  patches_;

    // Nearest reconstructed vertex accelerator for global projection seeds
    BvhAabb vertex_bvh_;
};

} // namespace fem::model
