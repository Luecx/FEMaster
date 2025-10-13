/**
 * @file keyword.h
 * @brief Declarative specification for keyword-line arguments of a DSL command.
 *
 * `KeywordSpec` allows command registrations to declare expected keyword arguments,
 * including canonical names, alternative spellings, optional defaults, and limited
 * value domains. The DSL engine can use this metadata to normalize parsed keys,
 * enforce required fields, and provide better diagnostics before variant execution.
 */

#pragma once

#include <initializer_list>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <utility>

namespace fem {
namespace dsl {

struct KeywordSpec {
    struct Entry {
        std::string canonical;                 ///< Canonical name stored in the `Keys` map.
        bool required = false;                 ///< Whether the key must be present.
        bool is_flag = false;                  ///< True if the key is a flag (no value expected).
        bool has_default = false;              ///< Whether a default value should be injected.
        std::string default_value;             ///< Default value (if `has_default`).
        std::vector<std::string> allowed;      ///< Optional list of allowed values.
        std::string doc;                       ///< Optional documentation string.
        std::vector<std::string> alternatives; ///< Alternative spellings that map to `canonical`.
    };

    static KeywordSpec make() { return {}; }

    KeywordSpec& key(const std::string& canonical) {
        touch_entry(canonical);
        last_key_ = canonical;
        auto& e = entries_.at(canonical);
        e.canonical = canonical;
        e.is_flag = false;
        return *this;
    }

    KeywordSpec& flag(const std::string& canonical) {
        touch_entry(canonical);
        last_key_ = canonical;
        auto& e = entries_.at(canonical);
        e.canonical = canonical;
        e.is_flag = true;
        e.required = false;
        e.has_default = false;
        e.default_value.clear();
        return *this;
    }

    const std::unordered_map<std::string, Entry>& entries() const { return entries_; }
    const std::unordered_map<std::string, std::string>& alias_map() const { return alias_to_canonical_; }

    KeywordSpec& doc(std::string text) {
        last_entry("doc").doc = std::move(text);
        return *this;
    }

    KeywordSpec& required() {
        last_entry("required").required = true;
        return *this;
    }

    KeywordSpec& optional() {
        auto& e = last_entry("optional");
        e.required = false;
        e.has_default = false;
        e.default_value.clear();
        return *this;
    }

    KeywordSpec& optional(const std::string& def) {
        auto& e = last_entry("optional");
        e.required = false;
        e.has_default = true;
        e.default_value = def;
        return *this;
    }

    KeywordSpec& allowed(std::initializer_list<std::string> vals) {
        auto& e = last_entry("allowed");
        e.allowed.assign(vals.begin(), vals.end());
        return *this;
    }

    KeywordSpec& allowed(const std::vector<std::string>& vals) {
        auto& e = last_entry("allowed");
        e.allowed = vals;
        return *this;
    }

    KeywordSpec& alternative(const std::string& alias) {
        if (last_key_.empty()) {
            throw std::runtime_error("KeywordSpec::alternative(): no active key");
        }
        register_alias(last_key_, alias);
        return *this;
    }

private:
    void touch_entry(const std::string& canonical) {
        auto it = entries_.find(canonical);
        if (it == entries_.end()) {
            Entry e;
            e.canonical = canonical;
            entries_.emplace(canonical, std::move(e));
        }
    }

    void register_alias(const std::string& canonical, const std::string& alias) {
        if (alias.empty())
            throw std::runtime_error("Keyword alias cannot be empty");

        if (alias == canonical)
            return;

        if (entries_.find(alias) != entries_.end()) {
            throw std::runtime_error("Alias '" + alias + "' conflicts with an existing canonical key");
        }

        auto [it, inserted] = alias_to_canonical_.emplace(alias, canonical);
        if (!inserted && it->second != canonical) {
            throw std::runtime_error("Alias '" + alias + "' already bound to a different key");
        }

        auto& entry = entries_.at(canonical);
        if (std::find(entry.alternatives.begin(), entry.alternatives.end(), alias) == entry.alternatives.end()) {
            entry.alternatives.push_back(alias);
        }
    }

    std::unordered_map<std::string, Entry> entries_;
    std::unordered_map<std::string, std::string> alias_to_canonical_;
    std::string last_key_;

    Entry& last_entry(const char* api) {
        if (last_key_.empty()) {
            throw std::runtime_error(std::string("KeywordSpec::") + api + "(): no active key");
        }
        return entries_.at(last_key_);
    }
};

} // namespace dsl
} // namespace fem
