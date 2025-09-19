// ----------------------------------------------------------------------------
#pragma once

#include "../../core/types_eig.h"
#include "../../core/logging.h"
#include <Eigen/Sparse>
#include <string>
#include <unordered_map>
#include <vector>
#include <type_traits>

namespace fem { namespace constraint {

struct MemItem {
    std::string label;   // e.g. "C_use", "Rsp", "T_full", "X_cols"
    std::size_t bytes{0};
};

class MemoryTracker {
public:
    void add(const std::string& label, std::size_t bytes) {
        totals_[label] += bytes;
    }

    template<class T>
    void add_vector(const std::string& label, const std::vector<T>& v) {
        add(label, v.capacity() * sizeof(T));
    }

    template<class K, class V>
    void add_unordered_map(const std::string& label, const std::unordered_map<K,V>& m) {
        // very rough approximation
        add(label, m.size() * (sizeof(std::pair<const K,V>) + 2*sizeof(void*)) );
    }

    template<class T>
    void add_nested_vector(const std::string& label, const std::vector<std::vector<T>>& vv) {
        std::size_t bytes = vv.capacity() * sizeof(std::vector<T>);
        for (const auto& v : vv) bytes += v.capacity() * sizeof(T);
        add(label, bytes);
    }

    template<class T1, class T2>
    void add_nested_vector_pair(const std::string& label,
                                const std::vector<std::vector<std::pair<T1,T2>>>& vv) {
        std::size_t bytes = vv.capacity() * sizeof(std::vector<std::pair<T1,T2>>);
        for (const auto& v : vv) bytes += v.capacity() * sizeof(std::pair<T1,T2>);
        add(label, bytes);
    }

    template<class Scalar, int Options, typename StorageIndex>
    void add_sparse(const std::string& label,
                    const Eigen::SparseMatrix<Scalar,Options,StorageIndex>& A) {
        // Compressed column storage estimate
        std::size_t nnz   = (std::size_t)A.nonZeros();
        std::size_t ncols = (std::size_t)A.cols();
        std::size_t bytes = nnz * (sizeof(Scalar) + sizeof(StorageIndex))
                          + (ncols + 1) * sizeof(StorageIndex);
        add(label, bytes);
    }

    std::vector<MemItem> items() const {
        std::vector<MemItem> out; out.reserve(totals_.size());
        for (auto& kv : totals_) out.push_back({kv.first, kv.second});
        return out;
    }

    std::size_t total_bytes() const {
        std::size_t s=0; for (auto& kv : totals_) s += kv.second; return s; }

    static void log_summary(const MemoryTracker& mt) {
        logging::info(true, "");
        logging::info(true, "Memory usage (live objects, approx):");
        for (auto& it : mt.items()) {
            logging::info(true, "  ", it.label, ": ", (unsigned long long)it.bytes, " B");
        }
        logging::info(true, "  Total: ", (unsigned long long)mt.total_bytes(), " B");
    }

private:
    std::unordered_map<std::string,std::size_t> totals_{};
};

}} // namespace fem::constraint

// ============================================================================
