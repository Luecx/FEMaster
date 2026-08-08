/**
 * @file bvh_aabb.cpp
 * @brief Implements the spatial axis-aligned bounding-box search hierarchy.
 *
 * The implementation provides validated AABB construction, model-independent
 * primitive insertion, balanced median hierarchy construction and iterative
 * spatial queries. No finite-element or model data is accessed in this file.
 *
 * Point and box queries first prune complete subtrees using node bounds and then
 * test each primitive box in a surviving leaf. Nearest lookup uses squared
 * point-to-AABB distance as a lower bound and visits the closer child first so a
 * useful upper bound is established early.
 *
 * @see BvhAabb
 * @see BvhAabb::Aabb
 *
 * @author Finn Eggers
 * @date 08.08.2026
 */

#include "bvh_aabb.h"

#include "../../../core/logging.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>

namespace fem {
namespace model {

// ============================================================================
// AABB construction and expansion
// ============================================================================

/**
 * @brief Creates the reversed-interval sentinel used for incremental expansion.
 *
 * Every lower coordinate starts at the largest representable value and every
 * upper coordinate at the lowest. The first finite point therefore replaces all
 * six bounds and produces a valid zero-volume box.
 *
 * @return Invalid axis-aligned box with reversed Cartesian intervals.
 */
BvhAabb::Aabb BvhAabb::Aabb::invalid() {
    return Aabb{};
}

/**
 * @brief Expands this box to include one global Cartesian point.
 *
 * Each lower bound is replaced by the smaller coordinate and each upper bound by
 * the larger coordinate. Applied to the invalid sentinel, the operation creates
 * a valid box at the supplied point.
 *
 * @param point Finite point in global Cartesian coordinates.
 */
void BvhAabb::Aabb::expand_point(const Vec3& point) {
    logging::error(point.allFinite(),
        "BvhAabb::Aabb::expand_point requires finite coordinates");

    for (int direction = 0; direction < 3; ++direction) {
        lo(direction) = std::min(lo(direction), point(direction));
        hi(direction) = std::max(hi(direction), point(direction));
    }
}

/**
 * @brief Expands this box by the complete extent of another valid box.
 *
 * Component-wise minima and maxima form the smallest axis-aligned union box.
 * The receiving box may still be the invalid construction sentinel; the supplied
 * box must already contain finite and ordered coordinate intervals.
 *
 * @param other Valid finite bounding box to include.
 */
void BvhAabb::Aabb::expand_aabb(const Aabb& other) {
    const bool finite = other.lo.allFinite() && other.hi.allFinite();
    const bool ordered =
           other.lo(0) <= other.hi(0)
        && other.lo(1) <= other.hi(1)
        && other.lo(2) <= other.hi(2);

    logging::error(finite,
        "BvhAabb::Aabb::expand_aabb requires finite box coordinates");
    logging::error(ordered,
        "BvhAabb::Aabb::expand_aabb requires ordered box coordinates");

    for (int direction = 0; direction < 3; ++direction) {
        lo(direction) = std::min(lo(direction), other.lo(direction));
        hi(direction) = std::max(hi(direction), other.hi(direction));
    }
}

/**
 * @brief Inflates a valid box uniformly in all Cartesian directions.
 *
 * The same non-negative radius is subtracted from every lower coordinate and
 * added to every upper coordinate. Inflation is deliberately performed on
 * primitive boxes before hierarchy construction so every internal node becomes
 * the exact union of already conservative leaf bounds.
 *
 * @param radius Non-negative broadphase margin in physical length units.
 */
void BvhAabb::Aabb::inflate(Precision radius) {
    const bool finite = lo.allFinite() && hi.allFinite() && std::isfinite(radius);
    const bool ordered =
           lo(0) <= hi(0)
        && lo(1) <= hi(1)
        && lo(2) <= hi(2);

    logging::error(finite,
        "BvhAabb::Aabb::inflate requires finite box coordinates and radius");
    logging::error(ordered,
        "BvhAabb::Aabb::inflate requires an ordered bounding box");
    logging::error(radius >= Precision(0),
        "BvhAabb::Aabb requires a non-negative inflation radius");

    for (int direction = 0; direction < 3; ++direction) {
        lo(direction) -= radius;
        hi(direction) += radius;
    }
}

// ============================================================================
// AABB geometric properties
// ============================================================================

/**
 * @brief Computes the midpoint of each Cartesian box interval.
 *
 * @return Global Cartesian centroid `(lo + hi) / 2` of a valid box.
 */
Vec3 BvhAabb::Aabb::centroid() const {
    return (lo + hi) * Precision(0.5);
}

/**
 * @brief Computes the three non-negative side lengths of a valid box.
 *
 * @return Cartesian box extent `hi - lo`.
 */
Vec3 BvhAabb::Aabb::extent() const {
    return hi - lo;
}

// ============================================================================
// AABB predicates and distance
// ============================================================================

/**
 * @brief Tests whether two closed axis-aligned boxes overlap.
 *
 * Separation in any Cartesian direction makes the boxes disjoint. Equality at
 * an interval boundary counts as intersection, which is required for a
 * conservative broadphase at touching primitive bounds.
 *
 * @param other Valid bounding box tested against this box.
 * @return `true` when the boxes overlap or touch in all three directions.
 */
bool BvhAabb::Aabb::intersects(const Aabb& other) const {
    return lo(0) <= other.hi(0) && hi(0) >= other.lo(0)
        && lo(1) <= other.hi(1) && hi(1) >= other.lo(1)
        && lo(2) <= other.hi(2) && hi(2) >= other.lo(2);
}

/**
 * @brief Tests containment in the closed Cartesian box intervals.
 *
 * @param point Point in global Cartesian coordinates.
 * @return `true` when every point coordinate lies between its box bounds.
 */
bool BvhAabb::Aabb::contains_point(const Vec3& point) const {
    return point(0) >= lo(0) && point(0) <= hi(0)
        && point(1) >= lo(1) && point(1) <= hi(1)
        && point(2) >= lo(2) && point(2) <= hi(2);
}

/**
 * @brief Computes squared Euclidean distance from a point to a closed AABB.
 *
 * A coordinate inside its interval contributes zero. An exterior coordinate
 * contributes the squared distance to its nearest interval endpoint. Summing
 * the three independent contributions yields the exact point-to-box distance
 * and a lower bound for every primitive contained by a hierarchy node.
 *
 * @param point Point in global Cartesian coordinates.
 * @return Squared distance to the closest point of this box.
 */
Precision BvhAabb::Aabb::squared_distance(const Vec3& point) const {
    Precision distance = Precision(0);

    for (int direction = 0; direction < 3; ++direction) {
        Precision delta = Precision(0);

        if (point(direction) < lo(direction)) {
            delta = lo(direction) - point(direction);
        } else if (point(direction) > hi(direction)) {
            delta = point(direction) - hi(direction);
        }

        distance += delta * delta;
    }

    return distance;
}

// ============================================================================
// BVH construction and reset
// ============================================================================

/**
 * @brief Constructs an empty hierarchy with uniform primitive inflation.
 *
 * Inflation is stored as construction configuration and applied when primitives
 * are inserted. A hierarchy root does not exist until at least one primitive has
 * been inserted and `finalize()` has completed.
 *
 * @param inflation Non-negative margin applied to every primitive box.
 */
BvhAabb::BvhAabb(Precision inflation)
    : inflation_(inflation) {
    logging::error(std::isfinite(inflation) && inflation >= Precision(0),
        "BvhAabb requires a finite non-negative inflation radius");
}

/**
 * @brief Removes all primitives, hierarchy nodes and reusable query output.
 *
 * The inflation radius and most recently selected leaf size remain unchanged so
 * the instance can be repopulated under the same configuration. The hierarchy
 * becomes invalid until new primitives are inserted and finalized.
 */
void BvhAabb::clear() {
    primitives_.clear();
    permutation_.clear();
    nodes_.clear();
    scratch_.clear();

    root_ = -1;
}

// ============================================================================
// Primitive insertion
// ============================================================================

/**
 * @brief Inserts one externally identified primitive bounding box.
 *
 * The finite ordered box is copied, inflated once by the configured radius and
 * assigned its centroid for later median partitioning. Insertion clears only
 * derived tree state; all previously inserted primitives remain available for
 * the next complete `finalize()` rebuild.
 *
 * @param id Opaque primitive identifier returned unchanged by queries.
 * @param box Primitive box before hierarchy-wide inflation.
 */
void BvhAabb::add_aabb(ID id, const Aabb& box) {
    const bool finite = box.lo.allFinite() && box.hi.allFinite();
    const bool ordered =
           box.lo(0) <= box.hi(0)
        && box.lo(1) <= box.hi(1)
        && box.lo(2) <= box.hi(2);

    logging::error(finite,
        "BvhAabb::add_aabb requires finite box coordinates");
    logging::error(ordered,
        "BvhAabb::add_aabb requires lower bounds not exceeding upper bounds");

    Primitive primitive;
    primitive.id  = id;
    primitive.box = box;

    primitive.box.inflate(inflation_);
    primitive.centroid = primitive.box.centroid();

    primitives_.push_back(primitive);

    // Existing node bounds and permutation do not cover the enlarged primitive
    // set. Invalidate them immediately so valid() cannot expose a stale tree.
    permutation_.clear();
    nodes_.clear();
    root_ = -1;
}

/**
 * @brief Inserts one point as a zero-volume primitive box.
 *
 * With zero configured inflation, nearest-AABB lookup over point primitives is
 * exactly a nearest-point lookup. Positive inflation turns every point into a
 * closed cube for conservative containment queries.
 *
 * @param id Opaque point identifier returned unchanged by queries.
 * @param point Finite point in global Cartesian coordinates.
 */
void BvhAabb::add_point(ID id, const Vec3& point) {
    logging::error(point.allFinite(),
        "BvhAabb::add_point requires finite point coordinates");

    Aabb box = Aabb::invalid();
    box.expand_point(point);
    add_aabb(id, box);
}

/**
 * @brief Rebuilds the binary hierarchy from all inserted primitives.
 *
 * Primitive indices are initialized in insertion order and recursively split at
 * their median centroid along the widest dimension of the current subtree.
 * This keeps the tree height logarithmic without modifying primitive storage.
 * Existing hierarchy data is discarded first. An empty primitive set remains a
 * valid construction state but produces no root and therefore `valid() == false`.
 *
 * @param leaf_size Positive maximum number of primitives stored in one leaf.
 */
void BvhAabb::finalize(int leaf_size) {
    const std::size_t maximum_primitive_count =
        static_cast<std::size_t>(std::numeric_limits<int>::max()) / 2;

    logging::error(leaf_size > 0,
        "BvhAabb requires leaf_size > 0");
    logging::error(primitives_.size() <= maximum_primitive_count,
        "BvhAabb contains too many primitives for integer hierarchy indexing");

    permutation_.clear();
    nodes_.clear();
    root_ = -1;

    leaf_size_ = leaf_size;

    if (primitives_.empty()) {
        return;
    }

    permutation_.resize(primitives_.size());
    std::iota(permutation_.begin(), permutation_.end(), 0);

    // A binary hierarchy over N non-empty leaves contains fewer than 2N nodes.
    nodes_.reserve(2 * primitives_.size());
    root_ = build_recursive(0, static_cast<int>(permutation_.size()));
}

// ============================================================================
// Point containment query
// ============================================================================

/**
 * @brief Returns primitive IDs whose individual AABB contains a query point.
 *
 * Node boxes first reject complete subtrees. Inside each surviving leaf, every
 * primitive is tested again because a leaf box is the union of several boxes and
 * may contain regions belonging to none of them. Results are not sorted and
 * follow traversal order.
 *
 * The supplied output vector is cleared before traversal. Without one, the
 * internal scratch vector is reused and the returned reference remains valid
 * only until the next scratch-backed query or hierarchy mutation.
 *
 * @param point Finite query point in global Cartesian coordinates.
 * @param out Optional caller-owned result vector used as reusable storage.
 * @return Reference to the caller-owned vector or internal scratch storage.
 */
const std::vector<ID>& BvhAabb::query_point(const Vec3& point,
                                            std::vector<ID>* out) const {
    logging::error(point.allFinite(),
        "BvhAabb::query_point requires finite point coordinates");

    std::vector<ID>& result = out != nullptr ? *out : scratch_;
    result.clear();

    if (root_ < 0) {
        return result;
    }

    // Median subdivision bounds the depth by the binary logarithm of the
    // integer-indexed primitive count; 64 entries exceed that possible depth.
    std::array<int, 64> stack;
    int stack_size = 0;

    stack[static_cast<std::size_t>(stack_size++)] = root_;

    while (stack_size > 0) {
        const int node_index = stack[static_cast<std::size_t>(--stack_size)];
        const Node& node = nodes_[static_cast<std::size_t>(node_index)];

        if (!node.box.contains_point(point)) {
            continue;
        }

        if (node.is_leaf()) {
            for (int i = node.begin; i < node.end; ++i) {
                const int primitive_index = permutation_[static_cast<std::size_t>(i)];
                const Primitive& primitive =
                    primitives_[static_cast<std::size_t>(primitive_index)];

                if (primitive.box.contains_point(point)) {
                    result.push_back(primitive.id);
                }
            }

            continue;
        }

        stack[static_cast<std::size_t>(stack_size++)] = node.left;
        stack[static_cast<std::size_t>(stack_size++)] = node.right;
    }

    return result;
}

// ============================================================================
// AABB intersection query
// ============================================================================

/**
 * @brief Returns primitive IDs whose individual AABB intersects a query box.
 *
 * Closed-box intersection at hierarchy nodes provides conservative subtree
 * pruning. Surviving leaves perform the same predicate on individual primitive
 * boxes so union-box false positives are not returned.
 *
 * The supplied output vector is cleared before traversal. Without one, the
 * internal scratch vector is reused and the returned reference remains valid
 * only until the next scratch-backed query or hierarchy mutation.
 *
 * @param box Finite ordered query box in global Cartesian coordinates.
 * @param out Optional caller-owned result vector used as reusable storage.
 * @return Reference to the caller-owned vector or internal scratch storage.
 */
const std::vector<ID>& BvhAabb::query_aabb(const Aabb& box,
                                           std::vector<ID>* out) const {
    const bool finite = box.lo.allFinite() && box.hi.allFinite();
    const bool ordered =
           box.lo(0) <= box.hi(0)
        && box.lo(1) <= box.hi(1)
        && box.lo(2) <= box.hi(2);

    logging::error(finite,
        "BvhAabb::query_aabb requires finite box coordinates");
    logging::error(ordered,
        "BvhAabb::query_aabb requires lower bounds not exceeding upper bounds");

    std::vector<ID>& result = out != nullptr ? *out : scratch_;
    result.clear();

    if (root_ < 0) {
        return result;
    }

    std::array<int, 64> stack;
    int stack_size = 0;

    stack[static_cast<std::size_t>(stack_size++)] = root_;

    while (stack_size > 0) {
        const int node_index = stack[static_cast<std::size_t>(--stack_size)];
        const Node& node = nodes_[static_cast<std::size_t>(node_index)];

        if (!node.box.intersects(box)) {
            continue;
        }

        if (node.is_leaf()) {
            for (int i = node.begin; i < node.end; ++i) {
                const int primitive_index = permutation_[static_cast<std::size_t>(i)];
                const Primitive& primitive =
                    primitives_[static_cast<std::size_t>(primitive_index)];

                if (primitive.box.intersects(box)) {
                    result.push_back(primitive.id);
                }
            }

            continue;
        }

        stack[static_cast<std::size_t>(stack_size++)] = node.left;
        stack[static_cast<std::size_t>(stack_size++)] = node.right;
    }

    return result;
}

// ============================================================================
// Nearest-AABB query
// ============================================================================

/**
 * @brief Finds the primitive box nearest to a query point.
 *
 * Squared point-to-box distance is an exact distance for a primitive box and a
 * lower bound for every primitive below a node union box. A node is discarded
 * when this lower bound cannot improve the current best primitive. At internal
 * nodes, the farther child is pushed first onto the LIFO stack so the closer
 * child is processed immediately and can tighten the upper bound early.
 *
 * The result concerns AABB distance rather than distance to any geometry enclosed
 * by that AABB. For zero-inflation point primitives, the result is an exact
 * nearest-point query. Equidistant primitives have unspecified ownership.
 *
 * @param point Finite query point in global Cartesian coordinates.
 * @return Nearest primitive ID, or `ID(-1)` for an invalid or empty hierarchy.
 */
ID BvhAabb::query_nearest(const Vec3& point) const {
    logging::error(point.allFinite(),
        "BvhAabb::query_nearest requires finite point coordinates");

    if (root_ < 0) {
        return ID(-1);
    }

    Precision best_distance = std::numeric_limits<Precision>::max();
    ID        best_id       = ID(-1);

    std::array<int, 64> stack;
    int stack_size = 0;

    stack[static_cast<std::size_t>(stack_size++)] = root_;

    while (stack_size > 0) {
        const int node_index = stack[static_cast<std::size_t>(--stack_size)];
        const Node& node = nodes_[static_cast<std::size_t>(node_index)];

        // The node box encloses all descendant primitive boxes. Its distance is
        // therefore a lower bound for every candidate in this complete subtree.
        if (node.box.squared_distance(point) >= best_distance) {
            continue;
        }

        if (node.is_leaf()) {
            for (int i = node.begin; i < node.end; ++i) {
                const int primitive_index = permutation_[static_cast<std::size_t>(i)];
                const Primitive& primitive =
                    primitives_[static_cast<std::size_t>(primitive_index)];
                const Precision distance = primitive.box.squared_distance(point);

                if (distance < best_distance) {
                    best_distance = distance;
                    best_id       = primitive.id;
                }
            }

            continue;
        }

        const Precision left_distance =
            nodes_[static_cast<std::size_t>(node.left)].box.squared_distance(point);
        const Precision right_distance =
            nodes_[static_cast<std::size_t>(node.right)].box.squared_distance(point);

        // The stack is LIFO: pushing the farther child first makes the closer
        // subtree execute next and usually improves pruning of the farther one.
        if (left_distance < right_distance) {
            if (right_distance < best_distance) {
                stack[static_cast<std::size_t>(stack_size++)] = node.right;
            }
            if (left_distance < best_distance) {
                stack[static_cast<std::size_t>(stack_size++)] = node.left;
            }
        } else {
            if (left_distance < best_distance) {
                stack[static_cast<std::size_t>(stack_size++)] = node.left;
            }
            if (right_distance < best_distance) {
                stack[static_cast<std::size_t>(stack_size++)] = node.right;
            }
        }
    }

    return best_id;
}

// ============================================================================
// Hierarchy state
// ============================================================================

/**
 * @brief Reports whether the current primitive set has a finalized tree.
 *
 * Empty hierarchies and hierarchies invalidated by later insertion have no root
 * and return `false`.
 *
 * @return `true` exactly when a valid non-empty root node exists.
 */
bool BvhAabb::valid() const {
    return root_ >= 0;
}

// ============================================================================
// Recursive hierarchy construction
// ============================================================================

/**
 * @brief Builds one subtree over a non-empty primitive permutation interval.
 *
 * The subtree box is the component-wise union of all inflated primitive boxes in
 * `[begin, end)`. Ranges no larger than `leaf_size_` terminate as leaves. Larger
 * ranges choose the widest subtree-box direction and use `std::nth_element` to
 * place the median centroid at the split without fully sorting either half.
 *
 * Recursive calls append nodes and may reallocate `nodes_`. The implementation
 * therefore retains only the stable parent index and reacquires the parent after
 * both child subtrees have been constructed.
 *
 * @param begin First permutation entry included in this subtree.
 * @param end One-past-last permutation entry included in this subtree.
 * @return Stable integer index of the newly constructed subtree root.
 */
int BvhAabb::build_recursive(int begin, int end) {
    const int count      = end - begin;
    const int node_index = static_cast<int>(nodes_.size());

    nodes_.push_back(Node{});

    // Construct the exact union bound for the current primitive range before
    // deciding whether it becomes a leaf or requires another median split.
    Aabb box = Aabb::invalid();

    for (int i = begin; i < end; ++i) {
        const int primitive_index = permutation_[static_cast<std::size_t>(i)];
        box.expand_aabb(primitives_[static_cast<std::size_t>(primitive_index)].box);
    }

    nodes_[static_cast<std::size_t>(node_index)].box = box;

    if (count <= leaf_size_) {
        Node& node = nodes_[static_cast<std::size_t>(node_index)];
        node.begin = begin;
        node.end   = end;
        return node_index;
    }

    // The widest box direction provides the largest available centroid spread
    // proxy and avoids introducing assumptions about primitive semantics.
    const Vec3 extent = box.extent();

    int axis = 0;
    if (extent(1) > extent(0)) axis = 1;
    if (extent(2) > extent(axis)) axis = 2;

    const int middle = begin + count / 2;

    std::nth_element(
        permutation_.begin() + begin,
        permutation_.begin() + middle,
        permutation_.begin() + end,
        [&](int left, int right) {
            return primitives_[static_cast<std::size_t>(left)].centroid(axis)
                 < primitives_[static_cast<std::size_t>(right)].centroid(axis);
        }
    );

    const int left  = build_recursive(begin, middle);
    const int right = build_recursive(middle, end);

    Node& node = nodes_[static_cast<std::size_t>(node_index)];
    node.left  = left;
    node.right = right;

    return node_index;
}

// ============================================================================
// Internal node state
// ============================================================================

/**
 * @brief Distinguishes primitive-range leaves from two-child internal nodes.
 *
 * @return `true` when neither child index is assigned.
 */
bool BvhAabb::Node::is_leaf() const {
    return left < 0 && right < 0;
}

} // namespace model
} // namespace fem
