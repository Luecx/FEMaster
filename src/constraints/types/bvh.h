/**
 * @file bvh.h
 * @brief Minimal BVH (AABB tree) for Tie broadphase: add element AABBs and query candidate IDs.
 *
 * Design goals:
 *  - Simple API: add_element(elem_id, node_coords, node_ids_ptr, n_nodes)
 *  - No geometry dependencies (works for surfaces and lines alike)
 *  - Conservative broadphase: element AABBs are inflated by 'inflate' (tie distance)
 *  - Query returns element IDs whose inflated AABB contains the query point
 *
 * Notes:
 *  - This is a broadphase accelerator. You still need your exact projection + clipping as narrowphase.
 *  - If master node coordinates change, rebuild the BVH.
 *
 * Typical usage:
 *  BvhAabb bvh(distance);
 *  for (each master element) bvh.add_element(id, node_coords, elem_nodes, n_nodes);
 *  bvh.finalize();
 *
 *  std::vector<Index> candidates;
 *  candidates.reserve(64);
 *  const auto& cand = bvh.query_point(p, &candidates);
 *
 * @author Finn Eggers
 * @date 16.01.2026
 */

#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

#include "../../core/core.h"
#include "../../core/types_eig.h"
#include "../../core/types_cls.h"

namespace fem {
namespace constraint {

class BvhAabb {
public:
    struct Aabb {
        Vec3 lo = Vec3::Constant(std::numeric_limits<Precision>::max());    // lower corner
        Vec3 hi = Vec3::Constant(std::numeric_limits<Precision>::lowest()); // upper corner

        static Aabb invalid() {
            return Aabb{};
        }

        void expand_point(const Vec3& p) {
            for (Index k = 0; k < 3; ++k) {
                if (p(k) < lo(k)) lo(k) = p(k);
                if (p(k) > hi(k)) hi(k) = p(k);
            }
        }

        void expand_aabb(const Aabb& o) {
            for (Index k = 0; k < 3; ++k) {
                if (o.lo(k) < lo(k)) lo(k) = o.lo(k);
                if (o.hi(k) > hi(k)) hi(k) = o.hi(k);
            }
        }

        Vec3 centroid() const { return (lo + hi) * Precision(0.5); }

        Vec3 extent() const { return (hi - lo); }

        void inflate(Precision r) {
            for (Index k = 0; k < 3; ++k) {
                lo(k) -= r;
                hi(k) += r;
            }
        }

        bool contains_point(const Vec3& p) const {
            return (p(0) >= lo(0) && p(0) <= hi(0) &&
                    p(1) >= lo(1) && p(1) <= hi(1) &&
                    p(2) >= lo(2) && p(2) <= hi(2));
        }

        bool intersects(const Aabb& other) const {
            return (lo(0) <= other.hi(0) && hi(0) >= other.lo(0) &&
                    lo(1) <= other.hi(1) && hi(1) >= other.lo(1) &&
                    lo(2) <= other.hi(2) && hi(2) >= other.lo(2));
        }
    };

    BvhAabb() = default;

    explicit BvhAabb(Precision inflate) : inflate_(inflate) {}

    void clear() {
        elems_.clear();
        perm_.clear();
        nodes_.clear();
        root_ = -1;
        scratch_.clear();
    }

    /**
     * @brief Add an element by listing its node IDs.
     *
     * @param elem_id     Element/geometry ID you want returned from queries.
     * @param node_coords Coordinate matrix-like object: node_coords(node_id, dim) -> Precision
     * @param node_ids    Pointer to node IDs
     * @param n_nodes     Number of nodes in this element
     */
    template <class NodeCoords>
    void add_element(Index elem_id, const NodeCoords& node_coords, const Index* node_ids, Index n_nodes) {
        if (node_ids == nullptr || n_nodes <= 0) {
            return;
        }

        Element e;
        e.id  = elem_id;
        e.box = Aabb::invalid();

        for (Index i = 0; i < n_nodes; ++i) {
            const Index nid = node_ids[i];

            Vec3 p;
            p(0) = node_coords(nid, 0);
            p(1) = node_coords(nid, 1);
            p(2) = node_coords(nid, 2);

            e.box.expand_point(p);
        }

        // Inflate leaf AABB once (broadphase margin)
        e.box.inflate(inflate_);
        e.centroid = e.box.centroid();

        elems_.push_back(e);
    }

    void add_aabb(Index elem_id, const Aabb& box) {
        Element e;
        e.id       = elem_id;
        e.box      = box;
        e.centroid = e.box.centroid();

        e.box.inflate(inflate_);

        elems_.push_back(e);
    }

    /**
     * @brief Build BVH from added elements.
     *
     * @param leaf_size Desired maximum number of primitives per leaf (typical 2..8).
     */
    void finalize(Index leaf_size = 4) {
        nodes_.clear();
        perm_.clear();
        root_ = -1;

        if (elems_.empty()) {
            return;
        }

        leaf_size_ = (leaf_size <= 0) ? 4 : leaf_size;

        perm_.resize(elems_.size());
        for (Index i = 0; i < static_cast<Index>(elems_.size()); ++i) {
            perm_[i] = i;
        }

        nodes_.reserve(2 * elems_.size());
        root_ = build_recursive(0, static_cast<Index>(perm_.size()));
    }

    /**
     * @brief Query candidate element IDs whose inflated AABB contains point p.
     *
     * @param p   Query point in global coordinates
     * @param out Optional output buffer to reuse allocations
     */
    const std::vector<Index>& query_point(const Vec3& p, std::vector<Index>* out = nullptr) const {
        std::vector<Index>& result = out ? *out : scratch_;
        result.clear();

        if (root_ < 0) {
            return result;
        }

        // Iterative stack. Depth is typically ~log2(N). 64 is safe for very large N.
        Index stack[64];
        Index sp = 0;
        stack[sp++] = root_;

        while (sp > 0) {
            const Index ni = stack[--sp];
            const Node& n = nodes_[ni];

            if (!n.box.contains_point(p)) {
                continue;
            }

            if (n.is_leaf()) {
                for (Index i = n.begin; i < n.end; ++i) {
                    const Index ei = perm_[i];
                    result.push_back(elems_[ei].id);
                }
            } else {
                if (n.left >= 0)  stack[sp++] = n.left;
                if (n.right >= 0) stack[sp++] = n.right;
            }
        }

        return result;
    }

    const std::vector<Index>& query_aabb(const Aabb& box, std::vector<Index>* out = nullptr) const {
        std::vector<Index>& result = out ? *out : scratch_;
        result.clear();

        if (root_ < 0) {
            return result;
        }

        Index stack[64];
        Index sp = 0;
        stack[sp++] = root_;

        while (sp > 0) {
            const Index ni = stack[--sp];
            const Node& n = nodes_[ni];

            if (!n.box.intersects(box)) {
                continue;
            }

            if (n.is_leaf()) {
                for (Index i = n.begin; i < n.end; ++i) {
                    const Index ei = perm_[i];
                    const Element& e = elems_[ei];
                    if (e.box.intersects(box)) {
                        result.push_back(e.id);
                    }
                }
            } else {
                if (n.left >= 0)  stack[sp++] = n.left;
                if (n.right >= 0) stack[sp++] = n.right;
            }
        }

        return result;
    }

    bool valid() const { return root_ >= 0; }

private:
    struct Element {
        Index id = Index(-1);
        Aabb  box;
        Vec3  centroid = Vec3::Zero();
    };

    struct Node {
        Aabb box;
        Index left  = -1;
        Index right = -1;
        Index begin = 0;
        Index end   = 0;

        bool is_leaf() const { return left < 0 && right < 0; }
    };

    Index build_recursive(Index begin, Index end) {
        const Index count = end - begin;

        const Index node_idx = static_cast<Index>(nodes_.size());
        nodes_.push_back(Node{});
        Node& node = nodes_.back();

        // Node AABB = union of element AABBs in range
        node.box = Aabb::invalid();
        for (Index i = begin; i < end; ++i) {
            const Index ei = perm_[i];
            node.box.expand_aabb(elems_[ei].box);
        }

        if (count <= leaf_size_) {
            node.begin = begin;
            node.end   = end;
            return node_idx;
        }

        // Split axis by widest extent
        const Vec3 ext = node.box.extent();
        Index axis = 0;
        if (ext(1) > ext(0)) axis = 1;
        if (ext(2) > ext(axis)) axis = 2;

        // Median split by centroid along axis
        const Index mid = begin + count / 2;
        std::nth_element(
            perm_.begin() + begin,
            perm_.begin() + mid,
            perm_.begin() + end,
            [&](Index a, Index b) {
                return elems_[a].centroid(axis) < elems_[b].centroid(axis);
            });

        const Index left  = build_recursive(begin, mid);
        const Index right = build_recursive(mid, end);

        node.left  = left;
        node.right = right;
        return node_idx;
    }

private:
    Precision inflate_   = Precision(0);
    Index leaf_size_     = 4;

    std::vector<Element> elems_;
    std::vector<Index>   perm_;
    std::vector<Node>    nodes_;
    Index root_          = -1;

    mutable std::vector<Index> scratch_;
};
} // namespace constraint
} // namespace fem
