
#pragma once
#include <vector>
#include "condition.h"

namespace fem { namespace reader2 {
struct And : Condition {
    std::vector<CondPtr> cs;
    explicit And(std::vector<CondPtr> v): cs(std::move(v)) {}
    bool eval(const ParentInfo& p, const Keys& k) const override {
        for (auto& c: cs) if (!c->eval(p,k)) return false; return true;
    }
};
struct Or : Condition {
    std::vector<CondPtr> cs;
    explicit Or(std::vector<CondPtr> v): cs(std::move(v)) {}
    bool eval(const ParentInfo& p, const Keys& k) const override {
        for (auto& c: cs) if (c->eval(p,k)) return true; return false;
    }
};
struct Not : Condition {
    CondPtr c; explicit Not(CondPtr v): c(std::move(v)) {}
    bool eval(const ParentInfo& p, const Keys& k) const override { return !c->eval(p,k); }
};
inline CondPtr AND(std::initializer_list<CondPtr> xs){ return std::make_shared<And>(std::vector<CondPtr>(xs)); }
inline CondPtr OR (std::initializer_list<CondPtr> xs){ return std::make_shared<Or >(std::vector<CondPtr>(xs)); }
inline CondPtr NOT(CondPtr x){ return std::make_shared<Not>(std::move(x)); }
}} // ns
