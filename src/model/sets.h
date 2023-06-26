#pragma once

#include "../core/types.h"

#include <string>
#include <unordered_map>
#include <utility>

namespace fem {
namespace model {
struct Sets {
    std::unordered_map<std::string, std::vector<ID>> m_sets;
    std::string                                      m_active;
    std::string                                      m_all;

    explicit Sets(std::string def)
        : m_all(std::move(def)){
        activate(m_all);
    }

    void activate(std::string name="") {
        if(name.empty()){
            name = m_all;
        }
        auto it = m_sets.find(name);
        if (it != m_sets.end()) {
            m_active = name;
        } else {
            // Create a new entry in the sets map
            m_sets.emplace(name, std::vector<ID>());
            m_active = name;
        }
    }

    std::vector<ID>& all() {
        return m_sets[m_all];
    }
    std::vector<ID>& ids() {
        return m_sets[m_active];
    }
    void add_id(ID id){
        all().push_back(id);
        ids().push_back(id);
    }
};
}    // namespace model
}    // namespace fem
