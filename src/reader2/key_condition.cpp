#include "key_condition.h"

namespace fem::reader2::condition {

// ---- KeyTest ----
KeyTest KeyTest::present(std::string k) { return {std::move(k), Type::Present, {}}; }
KeyTest KeyTest::missing(std::string k) { return {std::move(k), Type::Missing, {}}; }
KeyTest KeyTest::eq(std::string k, std::string v) {
    std::unordered_set<std::string> s{std::move(v)}; return {std::move(k), Type::Eq, std::move(s)};
}
KeyTest KeyTest::in(std::string k, std::initializer_list<std::string> vs) {
    return {std::move(k), Type::In, std::unordered_set<std::string>(vs.begin(), vs.end())};
}
KeyTest KeyTest::neq(std::string k, std::string v) {
    std::unordered_set<std::string> s{std::move(v)}; return {std::move(k), Type::Neq, std::move(s)};
}
KeyTest KeyTest::nin(std::string k, std::initializer_list<std::string> vs) {
    return {std::move(k), Type::Nin, std::unordered_set<std::string>(vs.begin(), vs.end())};
}

// ---- KeyCondition ----
KeyCondition KeyCondition::Leaf(KeyTest t) {
    KeyCondition e; e._op = Op::Leaf; e._test = std::move(t); return e;
}
KeyCondition KeyCondition::And(std::initializer_list<KeyCondition> xs) {
    KeyCondition e; e._op = Op::And; e._children.assign(xs.begin(), xs.end()); return e;
}
KeyCondition KeyCondition::Or(std::initializer_list<KeyCondition> xs) {
    KeyCondition e; e._op = Op::Or; e._children.assign(xs.begin(), xs.end()); return e;
}
KeyCondition KeyCondition::Not(KeyCondition x) {
    KeyCondition e; e._op = Op::Not; e._children = {std::move(x)}; return e;
}

bool KeyCondition::eval(const fem::reader2::Keyword::Map& kv) const {
    switch (_op) {
        case Op::Leaf: {
            const auto& t = _test;
            auto it = kv.find(t.key);
            switch (t.type) {
                case KeyTest::Type::Present: return it != kv.end();
                case KeyTest::Type::Missing: return it == kv.end();
                case KeyTest::Type::Eq:      return it != kv.end() && t.values.count(it->second);
                case KeyTest::Type::In:      return it != kv.end() && !t.values.empty() && t.values.count(it->second);
                case KeyTest::Type::Neq:     return it != kv.end() && !t.values.count(it->second);
                case KeyTest::Type::Nin:     return it != kv.end() && !t.values.count(it->second);
            }
            return false;
        }
        case Op::And:
            for (const auto& c : _children) if (!c.eval(kv)) return false;
            return true;
        case Op::Or:
            for (const auto& c : _children) if (c.eval(kv)) return true;
            return false;
        case Op::Not:
            return _children.empty() ? true : !(_children.front().eval(kv));
    }
    return false;
}

} // namespace fem::reader2::condition
