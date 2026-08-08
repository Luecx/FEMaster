/**
 * @file bvh_aabb.h
 * @brief Defines a spatial search hierarchy of axis-aligned bounding boxes.
 *
 * `BvhAabb` accelerates broadphase searches over geometric primitives identified
 * exclusively by an external `ID`. The hierarchy owns only copied bounding boxes,
 * centroids and tree indices. It has no knowledge of finite elements, model
 * nodes, contact constraints or the geometry represented by a primitive.
 *
 * Callers construct conservative primitive boxes, insert them with `add_aabb()`
 * or use `add_point()` for zero-volume point primitives, and explicitly complete
 * construction with `finalize()`. Point containment, box intersection and
 * nearest-box queries are available after finalization. Exact projection,
 * clipping and distance to the underlying geometry remain caller responsibilities.
 *
 * A uniform insertion-time inflation can provide a broadphase search margin.
 * Any insertion invalidates the current tree so stale hierarchy bounds can never
 * be reported as valid after the primitive set changes.
 *
 * @see BvhAabb
 * @see BvhAabb::Aabb
 *
 * @author Finn Eggers
 * @date 08.08.2026
 */

#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#include "../../../core/types_eig.h"

namespace fem {
namespace model {

/**
 * @brief Binary bounding-volume hierarchy over externally identified AABBs.
 *
 * Each primitive consists only of an `ID` and a finite, ordered axis-aligned
 * bounding box. The primitive ID is opaque: it may identify a surface patch, a
 * line, a point or any other object owned by the calling subsystem. Primitive
 * geometry is never dereferenced or interpreted by the hierarchy.
 *
 * `finalize()` recursively partitions the primitive permutation by median
 * centroid along the widest axis of each subtree. Primitive storage remains
 * stable while the permutation and node arrays form a balanced binary tree.
 * An empty or not-yet-finalized hierarchy is invalid and returns no query hits.
 *
 * A configured non-negative inflation radius is applied exactly once when a
 * primitive is inserted. Internal boxes are unions of those inflated primitive
 * boxes. `clear()` removes geometry and hierarchy state while retaining this
 * configuration. Adding another primitive after finalization invalidates the
 * tree until the next call to `finalize()`.
 *
 * Construction, insertion and public queries validate finite coordinates,
 * ordered box intervals, non-negative inflation and a positive leaf size through
 * `logging::error()` in the implementation. The low-level AABB predicates assume
 * valid operands because hierarchy traversal invokes them in inner loops.
 *
 * Point and AABB queries may write into caller-provided vectors. Without one,
 * they reuse an internal mutable scratch vector. Consequently, queries without
 * an external result buffer are not safe to execute concurrently on the same
 * instance. Read-only queries using separate external buffers do not share
 * mutable output state.
 */
class BvhAabb {
public:
    /**
     * @brief Three-dimensional axis-aligned bounding box in global coordinates.
     *
     * A valid box satisfies `lo(k) <= hi(k)` for all Cartesian directions and
     * contains only finite coordinates. The default state deliberately reverses
     * every interval and is therefore invalid. This sentinel supports incremental
     * construction: the first valid point passed to `expand_point()` establishes
     * a valid zero-volume box, and subsequent points or boxes enlarge it.
     *
     * Spatial predicates include the boundary. The squared point distance is zero
     * inside the box and otherwise measures the Euclidean distance to its closest
     * point. Callers must not query geometric properties of the invalid sentinel.
     */
    struct Aabb {
        // Cartesian lower and upper bounds; reversed in the invalid default state
        Vec3 lo = Vec3::Constant(std::numeric_limits<Precision>::max());
        Vec3 hi = Vec3::Constant(std::numeric_limits<Precision>::lowest());

        // Invalid construction state and incremental expansion
        [[nodiscard]] static Aabb invalid();

        void expand_point(const Vec3& point);
        void expand_aabb(const Aabb& other);
        void inflate(Precision radius);

        // Center and side lengths of a valid box
        [[nodiscard]] Vec3 centroid() const;
        [[nodiscard]] Vec3 extent() const;

        // Closed-set predicates and squared Euclidean point distance
        [[nodiscard]] bool intersects(const Aabb& other) const;
        [[nodiscard]] bool contains_point(const Vec3& point) const;
        [[nodiscard]] Precision squared_distance(const Vec3& point) const;
    };

    // Empty hierarchy construction and configuration-preserving reset
    BvhAabb() = default;
    explicit BvhAabb(Precision inflation);

    void clear();

    // Primitive insertion independent of the represented geometry
    void add_aabb(ID id, const Aabb& box);
    void add_point(ID id, const Vec3& point);

    // Explicit hierarchy construction over all currently inserted primitives
    void finalize(int leaf_size = 4);

    // Broadphase queries with optional caller-owned reusable result storage
    [[nodiscard]] const std::vector<ID>& query_point(const Vec3& point,
                                                     std::vector<ID>* out = nullptr) const;
    [[nodiscard]] const std::vector<ID>& query_aabb(const Aabb& box,
                                                    std::vector<ID>* out = nullptr) const;

    // Primitive whose AABB has minimum squared distance to the query point
    [[nodiscard]] ID query_nearest(const Vec3& point) const;

    // True only for a finalized, non-empty hierarchy not invalidated by insertion
    [[nodiscard]] bool valid() const;

private:
    /**
     * @brief Copied primitive data stored independently of caller geometry.
     *
     * `box` already includes the configured insertion-time inflation. `centroid`
     * is derived from that box and is used only to partition the hierarchy. The
     * opaque `id` is returned unchanged by successful queries.
     */
    struct Primitive {
        ID   id       = ID(-1);
        Aabb box;
        Vec3 centroid = Vec3::Zero();
    };

    /**
     * @brief One internal node or leaf of the binary hierarchy.
     *
     * Every node box is the union of all primitive boxes in its subtree. Internal
     * nodes reference two non-negative child indices and do not use their range.
     * Leaves have no children and reference the half-open permutation interval
     * `[begin, end)`. Node indices remain stable because recursive construction
     * stores indices rather than references across vector growth.
     */
    struct Node {
        Aabb box;

        int left  = -1;
        int right = -1;

        int begin = 0;
        int end   = 0;

        [[nodiscard]] bool is_leaf() const;
    };

    // Recursive median construction over a non-empty half-open permutation range
    [[nodiscard]] int build_recursive(int begin, int end);

    // Persistent construction configuration
    Precision inflation_ = Precision(0);
    int       leaf_size_ = 4;

    // Primitive data and hierarchy state rebuilt by finalize()
    std::vector<Primitive> primitives_;
    std::vector<int>       permutation_;
    std::vector<Node>      nodes_;

    int root_ = -1;

    // Shared only by queries that do not receive caller-owned result storage
    mutable std::vector<ID> scratch_;
};

} // namespace model
} // namespace fem
