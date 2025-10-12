/**
 * @file keys.h
 * @brief Lightweight key/value map for keyword line attributes (e.g. `*ELASTIC,E=...,NU=...`).
 *
 * The `Keys` utility wraps a `std::unordered_map<std::string,std::string>` and provides:
 *  - presence checks (`has`)
 *  - raw access (`raw`)
 *  - type-safe parsing (`get<T>`) with sensible behavior for booleans
 *  - string equality convenience (`equals`)
 *  - factory for extracting keys directly from a normalized DSL `Line`
 *    (`static Keys from_keyword_line(const Line&)`)
 *
 * ## Boolean semantics
 * A boolean key is interpreted as:
 *  - **true** if the key exists **without a value** (flag) *or* its value is one of:
 *    `1, T, TRUE, Y, YES, ON`
 *  - **false** if the key exists and its value is one of:
 *    `0, F, FALSE, N, NO, OFF`
 *  - If the key does not exist, `get<bool>(...)` returns the provided default.
 *  - If the key exists but the value is unrecognized (not in the sets above),
 *    it falls back to **false** (explicit beats implicit).
 *
 * Normal (non-boolean) parsing uses stream extraction (`operator>>`) and returns
 * the supplied default if extraction fails.
 *
 * @author
 *   Finn Eggers (project)
 * @date
 *   12.10.2025
 */

#pragma once
#include <string>
#include <unordered_map>
#include <sstream>
#include <type_traits>
#include <utility>
#include <algorithm>

#include "line.h"

namespace fem {
namespace dsl {

/**
 * @struct Keys
 * @brief Wrapper around a case-sensitive key/value dictionary for keyword attributes.
 *
 * The container stores keys (e.g. `"TYPE"`) and optional values (empty value means a flag).
 * Values are not normalized here (normalization should happen in your line tokenizer).
 */
struct Keys {
    /** Underlying storage (empty string value denotes a present flag without explicit value). */
    std::unordered_map<std::string, std::string> _kv;

    /**
     * @brief Builds a `Keys` map from a normalized keyword `Line`.
     *
     * The `Line` must be classified as `KEYWORD_LINE` and its internal normalized buffer
     * is expected to be free of spaces, uppercased, and without the leading `*`. The buffer
     * is split by commas into the command name and subsequent key/value pairs:
     *
     * ```
     * COMMAND,KEY=VALUE,FLAG
     * ```
     *
     * - Pairs with an `=` are added as `key -> value`.
     * - Bare items (no `=`) are stored as flags with an **empty value** (`""`).
     *
     * @param l A normalized DSL `Line` representing a keyword line.
     * @return Extracted key/value pairs (flags are stored with empty values).
     */
    static Keys from_keyword_line(const Line& l) {
        Keys k;
        if (l.type() != KEYWORD_LINE) return k;

        // Normalized keyword buffer looks like: "COMMAND,KEY=VAL,FLAG"
        std::string s = l.line();
        std::stringstream ss(s);

        // Discard the command (first comma-separated token)
        std::string cmd;
        (void)std::getline(ss, cmd, ',');

        // Parse subsequent comma-separated items as key/value or flags
        std::string item;
        while (std::getline(ss, item, ',')) {
            const auto pos = item.find('=');
            if (pos == std::string::npos) {
                k._kv[item] = "";                // flag (no value)
            } else {
                k._kv[item.substr(0, pos)] = item.substr(pos + 1);
            }
        }
        return k;
    }

    /**
     * @brief Returns true if the key is present.
     */
    bool has(const std::string& k) const {
        return _kv.find(k) != _kv.end();
    }

    /**
     * @brief Returns the raw value for a key, or a reference to a static empty string.
     *
     * This does not perform any conversion or validation.
     */
    const std::string& raw(const std::string& k) const {
        static const std::string empty;
        auto it = _kv.find(k);
        return it == _kv.end() ? empty : it->second;
    }

    /**
     * @brief Parses the value for key `k` as type `T`, or returns `def` if the key is missing.
     *
     * For non-bool types, standard stream extraction (`operator>>`) is used.
     * For `bool`, special semantics apply (see file header). If the key is present but
     * the value is unrecognized for a boolean, the function returns `false`.
     *
     * @tparam T   Target type (supports arithmetic types and `std::string`; `bool` has custom rules).
     * @param k    Key to look up.
     * @param def  Default to return if the key is not present.
     */
    template<class T>
    T get(const std::string& k, T def) const {
        auto it = _kv.find(k);
        if (it == _kv.end()) return def;

        if constexpr (std::is_same_v<T, bool>) {
            // Flag without value -> true
            if (it->second.empty()) return true;

            // Uppercase copy for case-insensitive matching
            std::string v = it->second;
            std::transform(v.begin(), v.end(), v.begin(),
                           [](unsigned char c){ return static_cast<char>(std::toupper(c)); });

            // True/false token sets
            auto is_true_token = [&]{
                return    v == "1"
                       || v == "T"
                       || v == "TRUE"
                       || v == "Y"
                       || v == "YES"
                       || v == "ON";
            };
            auto is_false_token = [&]{
                return    v == "0"
                       || v == "F"
                       || v == "FALSE"
                       || v == "N"
                       || v == "NO"
                       || v == "OFF";
            };

            if (is_true_token())  return true;
            if (is_false_token()) return false;

            // Fallback: try integer parse (non-zero -> true)
            std::istringstream ss(v);
            long long x = 0;
            if (ss >> x) return x != 0;

            // Unknown explicit value -> treat as false (explicit beats implicit).
            return false;
        } else if constexpr (std::is_same_v<T, std::string>) {
            return it->second;
        } else {
            std::istringstream ss(it->second);
            T out{};
            ss >> out;
            return ss.fail() ? def : out;
        }
    }

    /**
     * @brief Compares the stored value for `k` to `v` after converting `v` to string.
     *
     * @tparam T  Comparable value type convertible via stream insertion (`operator<<`).
     * @return `true` if the key exists and the stringified value equals the stored raw string.
     *
     * @warning This does **not** apply boolean semantics; it compares the raw string.
     *          For boolean logic, use `get<bool>(k, def)`.
     */
    template<class T>
    bool equals(const std::string& k, const T& v) const {
        auto it = _kv.find(k);
        if (it == _kv.end()) return false;
        std::ostringstream os;
        os << v;
        return it->second == os.str();
    }
};

} // namespace dsl
} // namespace fem
