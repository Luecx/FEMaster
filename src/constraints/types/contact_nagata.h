/**
 * @file contact_nagata.h
 * @brief Declares the Nagata-smoothed master geometry used by contact.
 *
 * The contact geometry tessellates the finite-element master surfaces into
 * three-node triangles without introducing additional geometric vertices.
 * Every patch corner therefore remains an existing model node. Adjacent patch
 * corners share one averaged normal across smooth edges and retain separate
 * normals across sharp edges.
 *
 * @author Finn Eggers
 */

#pragma once

#include "../../core/core.h"
#include "../../core/types_eig.h"
#include "../../data/field.h"
#include "../../data/region.h"
#include "../../model/geometry/surface/surface_interface.h"
#include "bvh.h"

#include <array>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

namespace fem::constraint {

enum class NagataEdgeType : std::uint8_t {
    Internal,
    Smooth,
    Sharp,
    Boundary
};

enum class NagataProjectionFeature : std::uint8_t {
    Interior,
    Edge,
    Corner
};

struct NagataPatch {
    ID id             = ID(-1);
    ID parent_surface = ID(-1);

    std::array<ID, 3> nodes{ID(-1), ID(-1), ID(-1)};
    std::array<Vec2, 3> parent_local{Vec2::Zero(), Vec2::Zero(), Vec2::Zero()};

    std::array<Vec3, 3> positions{Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};
    std::array<Vec3, 3> normals  {Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};

    // Local edge order: (0,1), (1,2), (2,0).
    std::array<Mat3, 3> curvature_maps{Mat3::Zero(), Mat3::Zero(), Mat3::Zero()};
    std::array<Vec3, 3> curvatures    {Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};

    std::array<Index, 3> normal_groups{Index(-1), Index(-1), Index(-1)};

    std::array<NagataPatch*, 3> neighbors{nullptr, nullptr, nullptr};
    std::array<NagataEdgeType, 3> edge_types{
        NagataEdgeType::Boundary,
        NagataEdgeType::Boundary,
        NagataEdgeType::Boundary
    };

    BvhAabb::Aabb bounds;
    Precision     area = Precision(0);
};

class NagataContactSurface {
public:
    struct Evaluation {
        bool valid = false;

        Vec3 position = Vec3::Zero();
        StaticMatrix<3, 2> first  = StaticMatrix<3, 2>::Zero();
        StaticMatrix<3, 3> second = StaticMatrix<3, 3>::Zero();

        Vec3      normal   = Vec3::Zero();
        Precision jacobian = Precision(0);

        // D_i = dx / dx_i with lagged nodal normals.
        std::array<Mat3, 3> position_derivative{
            Mat3::Zero(), Mat3::Zero(), Mat3::Zero()
        };

        // Derivatives of D_i with respect to r and s. These are also the direct
        // derivatives of the two patch tangents with respect to x_i.
        std::array<Mat3, 3> position_derivative_r{
            Mat3::Zero(), Mat3::Zero(), Mat3::Zero()
        };
        std::array<Mat3, 3> position_derivative_s{
            Mat3::Zero(), Mat3::Zero(), Mat3::Zero()
        };
    };

    struct Projection {
        bool valid = false;

        const NagataPatch* patch = nullptr;
        Vec2               local = Vec2::Zero();

        Vec3 position = Vec3::Zero();
        Vec3 normal   = Vec3::Zero();

        Precision distance = std::numeric_limits<Precision>::max();

        NagataProjectionFeature feature = NagataProjectionFeature::Interior;
        Index feature_index = Index(-1);
    };

    struct CornerRef {
        ID    patch  = ID(-1);
        Index corner = Index(-1);
    };

    NagataContactSurface(
        const model::SurfaceRegion&                       region,
        const std::vector<model::SurfaceInterface::Ptr>& surfaces,
        const model::Field&                               reference_positions,
        Precision feature_angle = Precision(0.78539816339744830962)
    );

    NagataContactSurface(const NagataContactSurface&)            = delete;
    NagataContactSurface& operator=(const NagataContactSurface&) = delete;
    NagataContactSurface(NagataContactSurface&&)                 = delete;
    NagataContactSurface& operator=(NagataContactSurface&&)      = delete;

    void update(
        const std::vector<model::SurfaceInterface::Ptr>& surfaces,
        const model::Field&                              positions,
        Precision                                        search_radius
    );

    [[nodiscard]] bool valid_patch(ID patch_id) const;

    [[nodiscard]] NagataPatch&       patch(ID patch_id);
    [[nodiscard]] const NagataPatch& patch(ID patch_id) const;

    [[nodiscard]] const std::vector<NagataPatch>& patches() const;
    [[nodiscard]] Index patch_count() const;
    [[nodiscard]] bool valid() const;

    [[nodiscard]] Evaluation evaluate(
        const NagataPatch& patch,
        const Vec2&        local
    ) const;

    [[nodiscard]] Projection project_on_patch(
        ID          patch_id,
        const Vec3& point,
        bool        clip
    ) const;

    [[nodiscard]] Projection walk(
        ID          start_patch,
        const Vec3& point
    ) const;

    [[nodiscard]] const std::vector<ID>& query_point(
        const Vec3& point,
        std::vector<ID>* buffer = nullptr
    ) const;

    [[nodiscard]] static bool in_bounds(
        const Vec2& local,
        Precision tolerance = Precision(0)
    );

    [[nodiscard]] static ID crossed_edge(const Vec2& local);
    [[nodiscard]] static Vec2 edge_direction(Index edge);

private:
    struct EdgeKey {
        ID first  = ID(-1);
        ID second = ID(-1);

        bool operator==(const EdgeKey& other) const {
            return first == other.first && second == other.second;
        }
    };

    struct EdgeKeyHash {
        std::size_t operator()(const EdgeKey& key) const;
    };

    struct EdgeOwner {
        ID    patch = ID(-1);
        Index edge  = Index(-1);
    };

    std::vector<NagataPatch>            patches_;
    std::vector<std::vector<CornerRef>> normal_groups_;
    BvhAabb                             bvh_;
    Precision                           feature_angle_ = Precision(0);

    void build_patches(
        const model::SurfaceRegion&                       region,
        const std::vector<model::SurfaceInterface::Ptr>& surfaces,
        const model::Field&                               reference_positions
    );

    void build_neighbors();
    void classify_edges();
    void build_normal_groups();

    void add_parent_patch(
        ID                              parent_surface,
        const model::SurfaceInterface& surface,
        const DynamicMatrix&            parent_local,
        const std::array<Index, 3>&     local_nodes
    );

    void update_positions(const model::Field& positions);
    void update_normals(
        const std::vector<model::SurfaceInterface::Ptr>& surfaces,
        const model::Field&                              positions
    );
    void update_patch_geometry();
    void update_bvh(Precision search_radius);

    [[nodiscard]] static std::array<Index, 2> edge_corners(Index edge);
    [[nodiscard]] static Index find_edge(
        const NagataPatch& patch,
        ID                 node_0,
        ID                 node_1
    );

    [[nodiscard]] static Vec3 linear_normal(const NagataPatch& patch);
    [[nodiscard]] static Mat3 curvature_map(
        const Vec3& normal_0,
        const Vec3& normal_1
    );

    [[nodiscard]] Projection project_stationary(
        const NagataPatch& patch,
        const Vec3&        point,
        bool               require_inside
    ) const;

    [[nodiscard]] Projection project_edge(
        const NagataPatch& patch,
        Index              edge,
        const Vec3&        point
    ) const;

    [[nodiscard]] Projection make_projection(
        const NagataPatch& patch,
        const Vec2&        local,
        const Vec3&        point
    ) const;

    [[nodiscard]] static Vec2 edge_local(Index edge, Precision coordinate);
    [[nodiscard]] static void classify_feature(Projection& projection);
};

} // namespace fem::constraint
