
#pragma once
#include <set>
#include "condition.h"

namespace fem { namespace reader2 {
struct ParentIs : Condition {
    std::set<std::string> allowed;
    explicit ParentIs(std::set<std::string> s): allowed(std::move(s)) {}
    bool eval(const ParentInfo& p, const Keys&) const override { return allowed.count(p.command)>0; }
};
struct ParentKeyEquals : Condition {
    std::string key; std::set<std::string> allowed;
    ParentKeyEquals(std::string k, std::set<std::string> v): key(std::move(k)), allowed(std::move(v)) {}
    bool eval(const ParentInfo& p, const Keys&) const override {
        auto raw = p.keys.raw(key); if (raw.empty()) return false; return allowed.count(raw)>0;
    }
};
struct ParentHasKey : Condition {
    std::string key; explicit ParentHasKey(std::string k): key(std::move(k)) {}
    bool eval(const ParentInfo& p, const Keys&) const override { return p.keys.has(key); }
};
}} // ns
