#pragma once

#include "../core/types.h"
#include "../core/logging.h"

#include <string>
#include <unordered_map>
#include <utility>

namespace fem {
namespace model {

// default set names
#define SET_NODE_ALL "NALL"
#define SET_ELEM_ALL "EALL"
#define SET_SURF_ALL "SFALL"
#define SET_LOAD_ALL "LALL"
#define SET_SUPP_ALL "SALL"

template<typename T = std::vector<ID>>
struct Sets {
    std::unordered_map<std::string, T> m_sets;
    std::string                        m_all;
    std::string                        m_active;

    template<typename... Args>
    Sets() : m_all(""), m_active("") {
    }

    template<typename... Args>
    explicit Sets(std::string def, Args&&... args)
        : m_all(std::move(def)) {
        activate(m_all, args...);
    }

    template<typename... Args>
    void activate(std::string name = "", Args&&... args) {
        if (name.empty()) {
            name = m_all;
        }
        auto it = m_sets.find(name);
        if (it != m_sets.end()) {
            m_active = name;
        } else {
            // Create a new entry in the sets map with forwarded arguments
            m_sets.emplace(std::piecewise_construct,
                           std::forward_as_tuple(name),
                           std::forward_as_tuple(std::forward<Args>(args)...));
            m_active = name;
        }
    }


    T& all() {
        logging::error(!m_all.empty(), "no all-set defined for given set");
        return m_sets[m_all];
    }
    T& current() {
        logging::error(!m_active.empty(), "active set is not defined for current _material");
        return m_sets[m_active];
    }
    bool has(const std::string& key){
        logging::error(!key.empty(), "empty keys are not valid for sets");
        return m_sets.find(key) != m_sets.end();
    }
    T& get(const std::string& key) {
        logging::error(has(key), key, "is not found within the set");
        return m_sets.at(key);
    }
    bool is_default_set() const{
        return m_all == m_active;
    }
    void add(ID id) {
        static_assert(std::is_same<T, std::vector<ID>>::value, "Cannot call add() for other than default type.");
        current().push_back(id);
        // dont add twice if active set is all
        if(m_all != m_active)
            all().push_back(id);
    }
};

}    // namespace model
}    // namespace fem
