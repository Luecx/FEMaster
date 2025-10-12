
#pragma once
#include <string>
#include <unordered_map>
#include <sstream>
#include <type_traits>
#include <utility>

namespace fem { namespace reader2 {
struct Keys {
    std::unordered_map<std::string,std::string> kv;
    bool has(const std::string& k) const { return kv.find(k)!=kv.end(); }
    const std::string& raw(const std::string& k) const {
        static const std::string empty;
        auto it = kv.find(k); return it==kv.end()? empty : it->second;
    }
    template<class T> T get(const std::string& k, T def) const {
        auto it = kv.find(k); if (it==kv.end()) return def;
        if constexpr (std::is_same_v<T,bool>) {
            if (it->second.empty()) return true;
        }
        std::istringstream ss(it->second); T out{}; ss>>out; return ss.fail()?def:out;
    }
    template<class T> bool equals(const std::string& k, const T& v) const {
        auto it = kv.find(k); if (it==kv.end()) return false;
        std::ostringstream os; os<<v; return it->second==os.str();
    }
};
}} // ns
