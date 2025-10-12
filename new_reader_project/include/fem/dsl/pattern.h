// ===== FILE: ./include/fem/dsl/pattern.h =====
/**
 * @file
 * @brief Defines the `Pattern` container describing token structures for segment parsing.
 *
 * A `Pattern` is an ordered collection of pattern elements (`Fixed<T,N>` etc.) that
 * collectively describe how a line (or set of lines) is to be parsed into typed values.
 *
 * Each element contributes a fixed number of tokens, and the engine enforces these
 * expectations when reading data lines. A pattern may optionally span multiple lines
 * when marked as multiline.
 *
 * Typical usage:
 * @code
 * Pattern::make()
 *   .fixed<int,4>().name("i").desc("indices")
 *   .fixed<double,2>().name("v").desc("values").on_missing(0.0)
 *   .allow_multiline();
 * @endcode
 *
 * The example defines a pattern that expects four integers followed by two doubles,
 * possibly split across multiple input lines. The fluent setters after `fixed(...)`
 * modify the last added element.
 */

#pragma once
#include <cstddef>
#include <vector>
#include <memory>
#include <stdexcept>
#include <typeinfo>
#include <string>
#include "pattern_element.h"

namespace fem {
namespace dsl {

/**
 * @class Pattern
 * @brief Aggregates one or more pattern elements defining a segment's data layout.
 *
 * The `Pattern` object acts as a builder for the expected token sequence within a segment.
 * Each call to `fixed<T,N>()` appends a new element describing `N` consecutive tokens of type `T`.
 * Immediately following calls to `name(...)`, `desc(...)`, `on_empty(...)`, and `on_missing(...)`
 * apply to the last added element (no temporaries needed).
 *
 * When `allow_multiline()` is set, the pattern may span across multiple input lines
 * until all required tokens have been read.
 */
struct Pattern {
    /**
     * @brief If true, the pattern may span multiple data lines.
     */
    bool _multiline = false;

    /**
     * @brief Ordered collection of pattern elements defining the token structure.
     */
    std::vector<std::shared_ptr<PatternElementBase>> _elems;

    /**
     * @brief Constructs a new empty pattern.
     *
     * @return A default-initialized `Pattern`.
     */
    static Pattern make() { return {}; }

    /**
     * @brief Enables multiline parsing for this pattern.
     *
     * @return Reference to `*this` for fluent chaining.
     */
    Pattern& allow_multiline() {
        _multiline = true;
        return *this;
    }

    /**
     * @brief Appends a fixed-size typed pattern element.
     *
     * @tparam T Value type (e.g., `int`, `double`, `std::string`).
     * @tparam N Number of expected tokens for this element.
     * @return Reference to `*this` for fluent chaining.
     */
    template<class T, std::size_t N>
    Pattern& fixed() {
        _elems.emplace_back(std::make_shared<Fixed<T, N>>());
        return *this;
    }

    /**
     * @brief Appends a single-value element equivalent to `Fixed<T,1>`.
     *
     * @tparam T Value type.
     * @return Reference to `*this` for fluent chaining.
     */
    template<class T>
    Pattern& one() {
        _elems.emplace_back(std::make_shared<Fixed<T, 1>>());
        return *this;
    }

    /**
     * @brief Sets a description for the last added element.
     *
     * @param d Description text to apply.
     * @return Reference to `*this` for fluent chaining.
     *
     * @throws std::runtime_error If the pattern has no elements or the element rejects the change.
     */
    Pattern& desc(const std::string& d) {
        auto* pe = last_or_throw("desc");
        if (!pe->set_desc(d)) {
            throw std::runtime_error("Pattern::desc(): last element does not accept desc()");
        }
        return *this;
    }

    /**
     * @brief Sets the base name for the last element's token(s).
     *
     * For array elements, printers may render `base1..baseN`.
     *
     * @param base Base identifier (no spaces).
     * @return Reference to `*this` for fluent chaining.
     *
     * @throws std::runtime_error If the pattern has no elements or the element rejects the change.
     */
    Pattern& name(const std::string& base) {
        auto* pe = last_or_throw("name");
        if (!pe->set_name(base)) {
            throw std::runtime_error("Pattern::name(): last element does not accept name()");
        }
        return *this;
    }

    /**
     * @brief Permits empty tokens and sets the default value for the last element.
     *
     * @tparam T Value type matching the last element's underlying type.
     * @param v Default value to substitute for empty tokens.
     * @return Reference to `*this` for fluent chaining.
     *
     * @throws std::runtime_error If the pattern has no elements or the type does not match.
     */
    template<class T>
    Pattern& on_empty(const T& v) {
        auto* pe = last_or_throw("on_empty");
        if (!pe->set_on_empty_any(&v, typeid(T))) {
            throw std::runtime_error("Pattern::on_empty(): type mismatch for last element");
        }
        return *this;
    }

    /**
     * @brief Permits missing tail tokens and sets the default value for the last element.
     *
     * @tparam T Value type matching the last element's underlying type.
     * @param v Default value to substitute for missing tokens.
     * @return Reference to `*this` for fluent chaining.
     *
     * @throws std::runtime_error If the pattern has no elements or the type does not match.
     */
    template<class T>
    Pattern& on_missing(const T& v) {
        auto* pe = last_or_throw("on_missing");
        if (!pe->set_on_missing_any(&v, typeid(T))) {
            throw std::runtime_error("Pattern::on_missing(): type mismatch for last element");
        }
        return *this;
    }

    /**
     * @brief Normalizes and completes a flat token vector against this pattern.
     *
     * For each element in order, this function:
     *  - replaces **empty tokens** with the element's default if `allow_empty()` is true,
     *    otherwise fails;
     *  - **appends the same default** for **missing** trailing tokens if `allow_missing()` is true,
     *    otherwise fails.
     *
     * If the input `toks` has **more tokens** than the pattern requires, the vector
     * is truncated to the required size.
     *
     * @param toks In/out token vector (already normalized by the tokenizer).
     * @param err_out On failure, receives a human-readable diagnostic.
     * @return `true` if normalization/completion succeeded, `false` otherwise.
     */
    bool normalize_and_complete_tokens(std::vector<std::string>& toks, std::string& err_out) const {
        std::size_t pos = 0;
        for (const auto& e_ptr : _elems) {
            const auto cnt = e_ptr->count();

            for (std::size_t i = 0; i < cnt; ++i) {
                const std::size_t idx = pos + i;

                if (idx < toks.size()) {
                    if (toks[idx].empty()) {
                        if (e_ptr->allow_empty()) {
                            toks[idx] = e_ptr->default_token_string();
                        } else {
                            err_out = "Empty token not allowed for element";
                            return false;
                        }
                    }
                } else {
                    if (e_ptr->allow_missing()) {
                        toks.push_back(e_ptr->default_token_string());
                    } else {
                        err_out = "Missing token(s) not allowed for element";
                        return false;
                    }
                }
            }

            pos += cnt;
        }

        if (toks.size() > pos) {
            toks.resize(pos);
        }
        return true;
    }

    /**
     * @brief Returns the total number of tokens this pattern expects.
     *
     * This equals the sum of `count()` for all contained elements.
     *
     * @return Total token count required by the pattern.
     */
    std::size_t required_tokens() const {
        std::size_t n = 0;
        for (const auto& e : _elems) {
            n += e->count();
        }
        return n;
    }

    /**
     * @brief Returns whether this pattern may span multiple data lines.
     *
     * @return True if multiline is enabled, false otherwise.
     */
    bool is_multiline() const {
        return _multiline;
    }

private:
    /**
     * @brief Returns a pointer to the last element or throws if the pattern is empty.
     *
     * @param api Name of the API for diagnostics.
     * @return Raw pointer to the last `PatternElementBase`.
     *
     * @throws std::runtime_error If no elements are present in the pattern.
     */
    PatternElementBase* last_or_throw(const char* api) {
        if (_elems.empty()) {
            throw std::runtime_error(std::string("Pattern::") + api + "(): no elements in pattern");
        }
        return _elems.back().get();
    }
};

} // namespace dsl
} // namespace fem
