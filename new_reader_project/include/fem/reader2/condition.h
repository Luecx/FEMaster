
#pragma once
#include <memory>
#include <string>
#include "keys.h"

namespace fem { namespace reader2 {
struct ParentInfo {
    std::string command;
    Keys        keys;
};
struct Condition {
    virtual ~Condition() = default;
    virtual bool eval(const ParentInfo& parent, const Keys& self_keys) const = 0;
};
using CondPtr = std::shared_ptr<Condition>;
}} // ns
