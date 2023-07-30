#pragma once

#include "../core/types.h"

#include <string>
#include <unordered_map>
#include <utility>

namespace fem {
namespace model {

template<typename T = std::vector<ID>>
struct Sets {
    std::unordered_map<std::string, T> m_sets;
    std::string                        m_active;
    std::string                        m_all;

    explicit Sets(std::string def)
        : m_all(std::move(def)) {
        activate(m_all);
    }

    void activate(std::string name = "") {
        if (name.empty()) {
            name = m_all;
        }
        auto it = m_sets.find(name);
        if (it != m_sets.end()) {
            m_active = name;
        } else {
            // Create a new entry in the sets map
            m_sets.emplace(name, T());
            m_active = name;
        }
    }

    T& all() {
        return m_sets[m_all];
    }
    T& current() {
        return m_sets[m_active];
    }
    void add(ID id) {
        all().push_back(id);
        current().push_back(id);
    }
};

}    // namespace model
}    // namespace fem
