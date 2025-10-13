#include "scope_condition.h"

namespace fem::reader2::condition {

ScopeCondition ScopeCondition::At(std::initializer_list<fem::reader2::Scope> scopes) {
    ScopeCondition c; c._op = Op::Atom; c._atom.loc = Location::At;
    c._atom.scopes.insert(scopes.begin(), scopes.end());
    return c;
}

ScopeCondition ScopeCondition::Within(std::initializer_list<fem::reader2::Scope> scopes) {
    ScopeCondition c; c._op = Op::Atom; c._atom.loc = Location::Within;
    c._atom.scopes.insert(scopes.begin(), scopes.end());
    return c;
}

ScopeCondition& ScopeCondition::Where(KeyCondition expr) {
    if (_op != Op::Atom) { _children.clear(); _op = Op::Atom; _atom = {}; }
    _atom.expr = std::move(expr);
    return *this;
}

ScopeCondition ScopeCondition::And(std::initializer_list<ScopeCondition> xs) {
    ScopeCondition c; c._op = Op::And; c._children.assign(xs.begin(), xs.end()); return c;
}

ScopeCondition ScopeCondition::Or(std::initializer_list<ScopeCondition> xs) {
    ScopeCondition c; c._op = Op::Or; c._children.assign(xs.begin(), xs.end()); return c;
}

ScopeCondition ScopeCondition::Not(ScopeCondition x) {
    ScopeCondition c; c._op = Op::Not; c._children = {std::move(x)}; return c;
}

bool ScopeCondition::eval_atom(const fem::reader2::Context& ctx, bool nearest_first) const {
    const auto& frames = ctx.frames();

    auto matches_frame = [&](const fem::reader2::ScopeFrame& f) -> bool {
        if (_atom.scopes.find(f.scope) == _atom.scopes.end()) return false;
        if (!_atom.expr.has_value()) return true;
        return _atom.expr->eval(f.kv);
    };

    if (_atom.loc == Location::At) {
        if (frames.empty()) return false;
        return matches_frame(frames.back());
    }

    if (nearest_first) {
        for (auto it = frames.rbegin(); it != frames.rend(); ++it)
            if (matches_frame(*it)) return true;
    } else {
        for (const auto& f : frames)
            if (matches_frame(f)) return true;
    }
    return false;
}

bool ScopeCondition::satisfied_by(const fem::reader2::Context& ctx, bool nearest_first) const {
    switch (_op) {
        case Op::Atom:
            return eval_atom(ctx, nearest_first);
        case Op::And:
            for (const auto& c : _children) if (!c.satisfied_by(ctx, nearest_first)) return false;
            return true;
        case Op::Or:
            for (const auto& c : _children) if (c.satisfied_by(ctx, nearest_first)) return true;
            return false;
        case Op::Not:
            return _children.empty() ? true : !(_children.front().satisfied_by(ctx, nearest_first));
    }
    return false;
}

} // namespace fem::reader2::condition
