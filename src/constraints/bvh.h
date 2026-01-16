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
 *  std::vector<ID> candidates;
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

#include "../core/core.h"
#include "../core/types_eig.h"
#include "../core/types_cls.h"

namespace fem {
namespace constraint {

class BvhAabb {
public:
    BvhAabb() = default;

    explicit BvhAabb(Precision inflate) : inflate_(inflate) {}

    void set_inflate(Precision inflate) { inflate_ = inflate; }

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
    void add_element(ID elem_id, const NodeCoords& node_coords, const ID* node_ids, int n_nodes) {
        if (node_ids == nullptr || n_nodes <= 0) {
            return;
        }

        Element e;
        e.id  = elem_id;
        e.box = Aabb::invalid();

        for (int i = 0; i < n_nodes; ++i) {
            const ID nid = node_ids[i];

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

    /**
     * @brief Build BVH from added elements.
     *
     * @param leaf_size Desired maximum number of primitives per leaf (typical 2..8).
     */
    void finalize(int leaf_size = 4) {
        nodes_.clear();
        perm_.clear();
        root_ = -1;

        if (elems_.empty()) {
            return;
        }

        leaf_size_ = (leaf_size <= 0) ? 4 : leaf_size;

        perm_.resize(elems_.size());
        for (std::size_t i = 0; i < elems_.size(); ++i) {
            perm_[i] = static_cast<int>(i);
        }

        nodes_.reserve(2 * elems_.size());
        root_ = build_recursive(0, static_cast<int>(perm_.size()));
    }

    /**
     * @brief Query candidate element IDs whose inflated AABB contains point p.
     *
     * @param p   Query point in global coordinates
     * @param out Optional output buffer to reuse allocations
     */
    const std::vector<ID>& query_point(const Vec3& p, std::vector<ID>* out = nullptr) const {
        std::vector<ID>& result = out ? *out : scratch_;
        result.clear();

        if (root_ < 0) {
            return result;
        }

        // Iterative stack. Depth is typically ~log2(N). 64 is safe for very large N.
        int stack[64];
        int sp = 0;
        stack[sp++] = root_;

        while (sp > 0) {
            const int ni = stack[--sp];
            const Node& n = nodes_[static_cast<std::size_t>(ni)];

            if (!n.box.contains_point(p)) {
                continue;
            }

            if (n.is_leaf()) {
                for (int i = n.begin; i < n.end; ++i) {
                    const int ei = perm_[static_cast<std::size_t>(i)];
                    result.push_back(elems_[static_cast<std::size_t>(ei)].id);
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
    struct Aabb {
        Vec3 lo; // lower corner
        Vec3 hi; // upper corner

        static Aabb invalid() {
            Aabb b;
            b.lo.setConstant(std::numeric_limits<Precision>::max());
            b.hi.setConstant(std::numeric_limits<Precision>::lowest());
            return b;
        }

        void expand_point(const Vec3& p) {
            for (int k = 0; k < 3; ++k) {
                if (p(k) < lo(k)) lo(k) = p(k);
                if (p(k) > hi(k)) hi(k) = p(k);
            }
        }

        void expand_aabb(const Aabb& o) {
            for (int k = 0; k < 3; ++k) {
                if (o.lo(k) < lo(k)) lo(k) = o.lo(k);
                if (o.hi(k) > hi(k)) hi(k) = o.hi(k);
            }
        }

        Vec3 centroid() const { return (lo + hi) * Precision(0.5); }

        Vec3 extent() const { return (hi - lo); }

        void inflate(Precision r) {
            for (int k = 0; k < 3; ++k) {
                lo(k) -= r;
                hi(k) += r;
            }
        }

        bool contains_point(const Vec3& p) const {
            return (p(0) >= lo(0) && p(0) <= hi(0) &&
                    p(1) >= lo(1) && p(1) <= hi(1) &&
                    p(2) >= lo(2) && p(2) <= hi(2));
        }
    };

    struct Element {
        ID   id = ID(-1);
        Aabb box;
        Vec3 centroid = Vec3::Zero();
    };

    struct Node {
        Aabb box;
        int left  = -1;
        int right = -1;
        int begin = 0;
        int end   = 0;

        bool is_leaf() const { return left < 0 && right < 0; }
    };

    int build_recursive(int begin, int end) {
        const int count = end - begin;

        const int node_idx = static_cast<int>(nodes_.size());
        nodes_.push_back(Node{});
        Node& node = nodes_.back();

        // Node AABB = union of element AABBs in range
        node.box = Aabb::invalid();
        for (int i = begin; i < end; ++i) {
            const int ei = perm_[static_cast<std::size_t>(i)];
            node.box.expand_aabb(elems_[static_cast<std::size_t>(ei)].box);
        }

        if (count <= leaf_size_) {
            node.begin = begin;
            node.end   = end;
            return node_idx;
        }

        // Split axis by widest extent
        const Vec3 ext = node.box.extent();
        int axis = 0;
        if (ext(1) > ext(0)) axis = 1;
        if (ext(2) > ext(axis)) axis = 2;

        // Median split by centroid along axis
        const int mid = begin + count / 2;
        std::nth_element(
            perm_.begin() + begin,
            perm_.begin() + mid,
            perm_.begin() + end,
            [&](int a, int b) {
                return elems_[static_cast<std::size_t>(a)].centroid(axis) <
                       elems_[static_cast<std::size_t>(b)].centroid(axis);
            });

        const int left  = build_recursive(begin, mid);
        const int right = build_recursive(mid, end);

        node.left  = left;
        node.right = right;
        return node_idx;
    }

private:
    Precision inflate_   = Precision(0);
    int leaf_size_       = 4;

    std::vector<Element> elems_;
    std::vector<int>     perm_;
    std::vector<Node>    nodes_;
    int root_            = -1;

    mutable std::vector<ID> scratch_;
};

} // namespace constraint
} // namespace fem
