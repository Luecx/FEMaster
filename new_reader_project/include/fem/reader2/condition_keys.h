
#pragma once
#include <set>
#include "condition.h"

namespace fem { namespace reader2 {
struct KeyPresent : Condition {
    std::string key; explicit KeyPresent(std::string k): key(std::move(k)) {}
    bool eval(const ParentInfo&, const Keys& k) const override { return k.has(key); }
};
struct KeyEquals : Condition {
    std::string key; std::set<std::string> allowed;
    KeyEquals(std::string k, std::set<std::string> v): key(std::move(k)), allowed(std::move(v)) {}
    bool eval(const ParentInfo&, const Keys& kk) const override { return allowed.count(kk.raw(key))>0; }
};
}} // ns
