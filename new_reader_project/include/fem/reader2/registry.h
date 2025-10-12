
#pragma once
#include <functional>
#include <unordered_map>
#include "command.h"

namespace fem { namespace reader2 {
struct Registry {
    std::unordered_map<std::string, Command> map_;
    void command(const std::string& name, const std::function<void(Command&)>& fn){
        auto it = map_.find(name);
        if (it==map_.end()) it = map_.emplace(name, Command{name}).first;
        fn(it->second);
    }
    const Command* find(const std::string& name) const {
        auto it = map_.find(name); return it==map_.end()? nullptr : &it->second;
    }
};
}} // ns
