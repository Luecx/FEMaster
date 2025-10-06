// reader2/condition.h
#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <initializer_list>
#include <sstream>
#include "keyword.h"  // needed for kw.has/kw.get in inline lambdas
namespace fem::reader2 {

class Keyword;

/**
 * \brief Predicate wrapper used to guard command plans.
 */
struct Condition {
    using Fn = std::function<bool(const Keyword&)>;

    Fn         test;
    std::string text;   // human-readable, used in documentation

    bool operator()(const Keyword& kw) const
    {
        return test ? test(kw) : false;
    }

    // ---- Builders ----
    /// \brief Condition that always evaluates to true.
    static Condition always();
    /// \brief Condition requiring that a keyword contains the specified key.
    static Condition has(std::string key);
    /// \brief Condition requiring key/value equality.
    static Condition eq(std::string key, std::string val);
    /// \brief Condition requiring the value to be within a set of candidates.
    static Condition in(std::string key, std::initializer_list<std::string> vals);
    /// \brief Condition requiring inequality with the provided value.
    static Condition neq(std::string key, std::string val);
    /// \brief Negate another condition.
    static Condition not_(Condition c);
    /// \brief Logical AND across many conditions.
    static Condition all(std::initializer_list<Condition> xs);
    /// \brief Logical OR across many conditions.
    static Condition any(std::initializer_list<Condition> xs);
};

inline Condition Condition::always() {
    return Condition{
        [](const Keyword&){ return true; },
        "ALWAYS"
    };
}

inline Condition Condition::has(std::string key) {
    std::string txt = std::move(key);
    return Condition{
        [k = txt](const Keyword& kw){ return kw.has(k); },
        "has(" + txt + ")"
    };
}
inline Condition Condition::eq(std::string key, std::string val) {
    std::string k = std::move(key), v = std::move(val);
    return Condition{
        [k, v](const Keyword& kw){
            try { return kw.get<std::string>(k, "") == v; } catch(...) { return false; }
        },
        k + "=" + v
    };
}
inline Condition Condition::in(std::string key, std::initializer_list<std::string> vals) {
    std::string k = std::move(key);
    std::vector<std::string> vs(vals.begin(), vals.end());
    std::ostringstream oss; oss << k << " in {";
    for (size_t i=0;i<vs.size();++i){ if(i) oss<<", "; oss<<vs[i]; }
    oss << "}";
    return Condition{
        [k, vs](const Keyword& kw){
            try {
                const auto v = kw.get<std::string>(k, "");
                for (auto& x: vs) if (v==x) return true;
                return false;
            } catch(...) { return false; }
        },
        oss.str()
    };
}
inline Condition Condition::neq(std::string key, std::string val) {
    std::string k=std::move(key), v=std::move(val);
    return Condition{
        [k, v](const Keyword& kw){
            try { return kw.get<std::string>(k, "") != v; } catch(...) { return false; }
        },
        k + "!=" + v
    };
}
inline Condition Condition::not_(Condition c) {
    return Condition{
        [c](const Keyword& kw){ return !c(kw); },
        "!(" + c.text + ")"
    };
}
inline Condition Condition::all(std::initializer_list<Condition> xs) {
    std::vector<Condition> cs(xs.begin(), xs.end());
    std::ostringstream oss; oss << "(";
    for (size_t i=0;i<cs.size();++i){ if(i) oss<<" && "; oss<<cs[i].text; }
    oss << ")";
    return Condition{
        [cs](const Keyword& kw){
            for (auto& c: cs) if (!c(kw)) return false;
            return true;
        },
        oss.str()
    };
}
inline Condition Condition::any(std::initializer_list<Condition> xs) {
    std::vector<Condition> cs(xs.begin(), xs.end());
    std::ostringstream oss; oss << "(";
    for (size_t i=0;i<cs.size();++i){ if(i) oss<<" || "; oss<<cs[i].text; }
    oss << ")";
    return Condition{
        [cs](const Keyword& kw){
            for (auto& c: cs) if (c(kw)) return true;
            return false;
        },
        oss.str()
    };
}

} // namespace fem::reader2
