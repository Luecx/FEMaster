// ===== FILE: ./include/fem/dsl/pattern_element.h =====
/**
 * @file
 * @brief Defines pattern element types used to describe segment data layouts.
 *
 * A pattern is composed of one or more elements that specify how many tokens
 * are expected and how they should be interpreted. This header provides:
 *
 *  - `PatternElementBase`: abstract base for pattern elements.
 *  - `Fixed<T,N>`: an element that reserves exactly `N` tokens of type `T`
 *                  (for `N==1` a scalar `T`, for `N>1` a `std::array<T,N>`).
 *
 * Defaults and tolerance:
 *  - `on_empty(v)`: permits empty tokens (`""`) and substitutes `v`.
 *  - `on_missing(v)`: permits missing tokens at the tail and substitutes `v`.
 *
 * Documentation hooks:
 *  - `desc(text)`: human-readable description of this element.
 *  - `name(base)`: base name for tokens (e.g., `"x"` → `x` or `x1,x2,...` for `N>1`).
 *  - `type_name()` returns a stable, human-readable type name (e.g., `"int"`, `"double"`).
 *
 * The engine consults these flags during parsing; help-printers use `desc(...)`,
 * `name(...)`, `count()`, and `type_name()` to render compact Data-Layouts.
 */

#pragma once
#include <array>
#include <string>
#include <memory>
#include <cstddef>
#include <utility>
#include <type_traits>
#include <typeinfo>
#include <sstream>

namespace fem {
namespace dsl {

/**
 * @brief Internal utilities to map C++ types to friendly names for docs.
 */
namespace detail {

/**
 * @tparam T
 * @brief Generic type-name mapper (default).
 *
 * @return A generic name for unrecognized types.
 */
template<class T>
struct TypeName {
    /**
     * @brief Returns a generic, human-readable type name.
     *
     * @return The string literal "value".
     */
    static const char* get() { return "value"; }
};

/**
 * @brief Specialization for int.
 */
template<>
struct TypeName<int> {
    /**
     * @brief Returns the human-readable type name.
     *
     * @return The string literal "int".
     */
    static const char* get() { return "int"; }
};

/**
 * @brief Specialization for long.
 */
template<>
struct TypeName<long> {
    /**
     * @brief Returns the human-readable type name.
     *
     * @return The string literal "long".
     */
    static const char* get() { return "long"; }
};

/**
 * @brief Specialization for unsigned.
 */
template<>
struct TypeName<unsigned> {
    /**
     * @brief Returns the human-readable type name.
     *
     * @return The string literal "unsigned".
     */
    static const char* get() { return "unsigned"; }
};

/**
 * @brief Specialization for float.
 */
template<>
struct TypeName<float> {
    /**
     * @brief Returns the human-readable type name.
     *
     * @return The string literal "float".
     */
    static const char* get() { return "float"; }
};

/**
 * @brief Specialization for double.
 */
template<>
struct TypeName<double> {
    /**
     * @brief Returns the human-readable type name.
     *
     * @return The string literal "double".
     */
    static const char* get() { return "double"; }
};

/**
 * @brief Specialization for std::string.
 */
template<>
struct TypeName<std::string> {
    /**
     * @brief Returns the human-readable type name.
     *
     * @return The string literal "string".
     */
    static const char* get() { return "string"; }
};

} // namespace detail

/**
 * @class PatternElementBase
 * @brief Abstract base class for pattern elements.
 *
 * All concrete pattern elements report how many tokens they consume via `count()`.
 * For documentation, elements can expose a type name, a human-readable description,
 * and an optional base-token name used by help printers.
 *
 * It also provides type-erased setters so `Pattern` can fluently modify the
 * last added element (desc, name, on_empty, on_missing) without exposing concrete types.
 *
 * Additionally, minimal query accessors (`allow_empty`, `allow_missing`, `default_token_string`)
 * are provided so the engine can apply defaults uniformly without knowing `T`.
 */
struct PatternElementBase {
    /**
     * @brief Virtual destructor.
     */
    virtual ~PatternElementBase() = default;

    /**
     * @brief Returns the number of tokens this element contributes to the pattern.
     *
     * @return Number of tokens reserved by this element.
     */
    virtual std::size_t count() const = 0;

    /**
     * @brief Returns a stable, human-readable type name.
     *
     * @return A string literal such as "int", "double", or "string".
     */
    virtual const char* type_name() const = 0;

    /**
     * @brief Returns the element's human-readable description.
     *
     * @return Description text or an empty string reference if unset.
     */
    virtual const std::string& description() const {
        static const std::string kEmpty;
        return kEmpty;
    }

    /**
     * @brief Returns the base name for tokens.
     *
     * If non-empty and `count() > 1`, printers may render `base1, base2, ..., baseN`.
     *
     * @return Base name or an empty string reference if unset.
     */
    virtual const std::string& name_base() const {
        static const std::string kEmpty;
        return kEmpty;
    }

    /**
     * @brief Type-erased setter for a human-readable description.
     *
     * @param d Description to set.
     * @return True on success, false if unsupported by the implementation.
     */
    virtual bool set_desc(const std::string& d) {
        (void)d;
        return false;
    }

    /**
     * @brief Type-erased setter for a base name used to label tokens.
     *
     * @param n Base name to set.
     * @return True on success, false if unsupported by the implementation.
     */
    virtual bool set_name(const std::string& n) {
        (void)n;
        return false;
    }

    /**
     * @brief Type-erased setter to permit empty tokens and set a default value.
     *
     * @param value Pointer to a value of the element's underlying type.
     * @param ti Type info describing the pointed-to value.
     * @return True if the type matches and the value was applied, otherwise false.
     */
    virtual bool set_on_empty_any(const void* value, const std::type_info& ti) {
        (void)value;
        (void)ti;
        return false;
    }

    /**
     * @brief Type-erased setter to permit missing tail tokens and set a default value.
     *
     * @param value Pointer to a value of the element's underlying type.
     * @param ti Type info describing the pointed-to value.
     * @return True if the type matches and the value was applied, otherwise false.
     */
    virtual bool set_on_missing_any(const void* value, const std::type_info& ti) {
        (void)value;
        (void)ti;
        return false;
    }

    /**
     * @brief Returns whether empty tokens are permitted for this element.
     *
     * @return True if `on_empty(...)` has been set, otherwise false.
     */
    virtual bool allow_empty() const { return false; }

    /**
     * @brief Returns whether missing tail tokens are permitted for this element.
     *
     * @return True if `on_missing(...)` has been set, otherwise false.
     */
    virtual bool allow_missing() const { return false; }

    /**
     * @brief Returns the stringified default token to use for both empty/missing positions.
     *
     * Implementations may store the default in typed form and stringify here.
     *
     * @return Default token string (may be empty if unset).
     */
    virtual std::string default_token_string() const { return {}; }
};

/**
 * @class Fixed
 * @brief Pattern element that reserves exactly `N` tokens of type `T`.
 *
 * Semantics:
 * - When `N == 1`, this element maps to a single value of type `T`.
 * - When `N > 1`, this element maps to `std::array<T,N>`.
 * - `desc(text)` attaches a human-readable description (for diagnostics/docs).
 * - `name(base)` sets an optional base name for the token(s): `base` or `base1..baseN`.
 * - `on_empty(v)` allows empty tokens (`""`) and replaces them with `v`.
 * - `on_missing(v)` allows missing tail tokens and fills them with `v`.
 *
 * The engine is responsible for enforcing and applying these flags during parsing.
 *
 * @tparam T Value type for tokens (e.g., `int`, `double`, `std::string`).
 * @tparam N Number of tokens reserved by this element (compile-time constant).
 */
template<class T, std::size_t N>
struct Fixed : PatternElementBase {
    /**
     * @brief Human-readable description used in diagnostics/help output (optional).
     */
    std::string _description;

    /**
     * @brief Optional base name used to label tokens (e.g., "x" → x / x1..xN).
     */
    std::string _name_base;

    /**
     * @brief If true, empty tokens are permitted and substituted with defaults.
     */
    bool _allow_empty = false;

    /**
     * @brief If true, missing tail tokens are permitted and substituted with defaults.
     */
    bool _allow_missing = false;

    /**
     * @brief Default used for single-value substitution.
     */
    T _default_single{};

    /**
     * @brief Default used to fill arrays (`N > 1`) or multiple empties.
     */
    std::array<T, N> _default_array{};

    /**
     * @brief Sets a short description used for docs and error messages.
     *
     * @param d Free-form text; kept verbatim.
     * @return Reference to `*this` for fluent chaining.
     */
    Fixed& desc(std::string d) {
        _description = std::move(d);
        return *this;
    }

    /**
     * @brief Sets a base name for this element's token(s).
     *
     * For `N==1`, printers will render `base`. For `N>1`, `base1, base2, ..., baseN`.
     *
     * @param base Base identifier (no spaces).
     * @return Reference to `*this` for fluent chaining.
     */
    Fixed& name(std::string base) {
        _name_base = std::move(base);
        return *this;
    }

    /**
     * @brief Permits empty tokens and substitutes the provided default.
     *
     * If a token is present but empty (`""`), the engine may replace it with `v`.
     * Also initializes the internal default array uniformly with `v`.
     *
     * @param v Default value to substitute for empty tokens.
     * @return Reference to `*this` for fluent chaining.
     */
    Fixed& on_empty(const T& v) {
        _allow_empty    = true;
        _default_single = v;
        _default_array.fill(v);
        return *this;
    }

    /**
     * @brief Permits missing tail tokens and substitutes the provided default.
     *
     * If fewer than `N` tokens are provided, the engine may fill the trailing
     * positions with `v`. Also initializes the internal default array uniformly with `v`.
     *
     * @param v Default value to substitute for missing tokens.
     * @return Reference to `*this` for fluent chaining.
     */
    Fixed& on_missing(const T& v) {
        _allow_missing  = true;
        _default_single = v;
        _default_array.fill(v);
        return *this;
    }

    /**
     * @brief Returns the number of tokens reserved by this element.
     *
     * @return The constant `N`.
     */
    std::size_t count() const override {
        return N;
    }

    /**
     * @brief Returns a friendly, stable type name for help output.
     *
     * @return A string literal describing the token type (e.g., "int").
     */
    const char* type_name() const override {
        return detail::TypeName<T>::get();
    }

    /**
     * @brief Returns the stored human-readable description.
     *
     * @return Description string reference (may be empty).
     */
    const std::string& description() const override {
        return _description;
    }

    /**
     * @brief Returns the base name used to label this element's token(s).
     *
     * @return Base name string reference (may be empty).
     */
    const std::string& name_base() const override {
        return _name_base;
    }

    /**
     * @brief Type-erased setter for the element description.
     *
     * @param d Description to set.
     * @return Always true for Fixed elements.
     */
    bool set_desc(const std::string& d) override {
        _description = d;
        return true;
    }

    /**
     * @brief Type-erased setter for the base name.
     *
     * @param n Base name to set.
     * @return Always true for Fixed elements.
     */
    bool set_name(const std::string& n) override {
        _name_base = n;
        return true;
    }

    /**
     * @brief Type-erased setter to permit empty tokens with default value.
     *
     * Applies only if `ti` matches the underlying type `T`.
     *
     * @param value Pointer to a value of type `T`.
     * @param ti Type info for the pointed-to value.
     * @return True if the type matches and the default was applied, false otherwise.
     */
    bool set_on_empty_any(const void* value, const std::type_info& ti) override {
        if (ti == typeid(T)) {
            on_empty(*reinterpret_cast<const T*>(value));
            return true;
        }
        return false;
    }

    /**
     * @brief Type-erased setter to permit missing tokens with default value.
     *
     * Applies only if `ti` matches the underlying type `T`.
     *
     * @param value Pointer to a value of type `T`.
     * @param ti Type info for the pointed-to value.
     * @return True if the type matches and the default was applied, false otherwise.
     */
    bool set_on_missing_any(const void* value, const std::type_info& ti) override {
        if (ti == typeid(T)) {
            on_missing(*reinterpret_cast<const T*>(value));
            return true;
        }
        return false;
    }

    /**
     * @brief Returns whether empty tokens are permitted for this element.
     *
     * @return True if `on_empty(...)` has been set, otherwise false.
     */
    bool allow_empty() const override {
        return _allow_empty;
    }

    /**
     * @brief Returns whether missing tail tokens are permitted for this element.
     *
     * @return True if `on_missing(...)` has been set, otherwise false.
     */
    bool allow_missing() const override {
        return _allow_missing;
    }

    /**
     * @brief Returns the stringified default token used for both empty/missing positions.
     *
     * @return Default token string created from `_default_single`.
     */
    std::string default_token_string() const override {
        std::ostringstream os;
        os << _default_single;
        return os.str();
    }
};

} // namespace dsl
} // namespace fem
