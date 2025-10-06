#pragma once
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace fem::reader2 {

struct KeyRules {
    struct Rule {
        bool required = false;
        std::optional<std::string> def;          // default (if not required)
        std::unordered_set<std::string> allowed; // allowed choices (already uppercased by parser)
        std::string description;                  // human doc
    };

    // name -> rule (names are already uppercased by line parser)
    std::unordered_map<std::string, Rule> rules;

    KeyRules& describe(std::string name, std::string text){
        rules[std::move(name)].description = std::move(text); return *this;
    }
    KeyRules& require(std::string name) {
        rules[std::move(name)].required = true; return *this;
    }
    KeyRules& optional(std::string name, std::optional<std::string> def = std::nullopt) {
        auto& r = rules[std::move(name)];
        r.required = false; r.def = std::move(def); return *this;
    }
    KeyRules& allowed(std::string name, std::initializer_list<std::string> xs) {
        auto& r = rules[name];
        r.allowed.insert(xs.begin(), xs.end()); return *this;
    }

    // Validate & materialize a final key map (string->string), applying defaults & enums.
    // Throws with a helpful message on violation.
    std::unordered_map<std::string, std::string>
    apply(const std::unordered_map<std::string,std::string>& input) const {
        std::unordered_map<std::string, std::string> out = input; // start with given keys
        for (const auto& [k, r] : rules) {
            auto it = out.find(k);
            if (it == out.end()) {
                if (r.required) {
                    throw std::runtime_error("Missing required key: " + k);
                }
                if (r.def) out[k] = *r.def;
                continue;
            }
            if (!r.allowed.empty() && !r.allowed.count(it->second)) {
                std::ostringstream oss;
                oss << "Invalid value '" << it->second << "' for key " << k << ". Allowed: {";
                bool first=true; for (auto& v : r.allowed){ if(!first) oss<<","; first=false; oss<<v; } oss<<"}";
                throw std::runtime_error(oss.str());
            }
        }
        return out;
    }

    // Typed getter helper
    template<class T>
    static T convert(const std::string& s){
        if constexpr (std::is_same_v<T,std::string>) return s;
        std::istringstream iss(s); T v{}; iss >> v;
        if (iss.fail()) throw std::runtime_error("Key conversion failed for value '"+s+"'");
        return v;
    }
};

} // namespace fem::reader2
