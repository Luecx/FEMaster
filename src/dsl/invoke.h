/**
 * @file detail/invoke.h
 * @brief Internal utilities for callable introspection and token-to-type conversion.
 *
 * This header is private to the DSL module. It provides:
 *  - `function_traits<F>` to extract a callable's parameter tuple (works with lambdas).
 *  - `is_std_array<T>` to detect `std::array<...,N>` arguments.
 *  - Scalar/array parsers that convert string tokens to typed values.
 *  - Token count helpers to validate lambda signatures against a `Pattern`.
 *  - `Invoker<Us...>` which applies a callable to tokens converted as `Us...`.
 *
 * Supported parameter types:
 *  - arithmetic scalars (`int`, `long`, `unsigned`, `float`, `double`, …)
 *  - `std::string` (returned as-is; adjust tokenizer if case sensitivity matters)
 *  - `std::array<T,N>` (consumes `N` tokens of base type `T`)
 *
 * @warning This header is implementation detail of the DSL. Do not include it directly
 *          from user code. Prefer the public DSL interfaces.
 *
 * @date 12.10.2025
 * @author
 *   Finn Eggers (project)
 */

#pragma once
#include <tuple>
#include <type_traits>
#include <array>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <functional>

namespace fem {
namespace dsl {
namespace detail {

/**
 * @brief Traits to extract the parameter tuple type from a callable (lambdas included).
 */
template<typename T>
struct function_traits;

// function pointer
template<typename R, typename... Args>
struct function_traits<R(*)(Args...)> {
    using args_tuple = std::tuple<Args...>;
};

// function reference
template<typename R, typename... Args>
struct function_traits<R(&)(Args...)> : function_traits<R(*)(Args...)> {};

// std::function
template<typename R, typename... Args>
struct function_traits<std::function<R(Args...)>> : function_traits<R(*)(Args...)> {};

// member function pointer (lambda operator())
template<typename C, typename R, typename... Args>
struct function_traits<R(C::*)(Args...) const> : function_traits<R(*)(Args...)> {};

// generic callable (lambda, functor): fall back to operator()
template<typename F>
struct function_traits : function_traits<decltype(&F::operator())> {};

/** @brief Trait for detecting `std::array`. */
template<class T> struct is_std_array : std::false_type {};
template<class U, std::size_t N> struct is_std_array<std::array<U,N>> : std::true_type {};
template<class T> inline constexpr bool is_std_array_v = is_std_array<T>::value;

/**
 * @brief Parses an arithmetic token into type `T`.
 * @tparam T Arithmetic type.
 * @throws std::runtime_error if parsing fails.
 */
template<class T>
inline T parse_arithmetic(const std::string& s) {
    static_assert(std::is_arithmetic_v<T>, "parse_arithmetic requires arithmetic type");
    std::istringstream ss(s);
    T v{};
    ss >> v;
    if (ss.fail()) throw std::runtime_error("parse_arithmetic failed for token: " + s);
    return v;
}

/**
 * @brief Parses a single token as a scalar `T` (`std::string` or arithmetic).
 */
template<class T>
inline T parse_scalar(const std::string& s) {
    if constexpr (std::is_same_v<T, std::string>) {
        // NOTE: If your tokenizer uppercases values, adjust there when case-sensitive strings are needed.
        return s;
    } else if constexpr (std::is_arithmetic_v<T>) {
        return parse_arithmetic<T>(s);
    } else {
        static_assert(!sizeof(T), "Unsupported scalar type in DSL invocation");
    }
}

/**
 * @brief Parses `N` tokens into `std::array<U,N>`, advancing the index.
 */
template<class U, std::size_t N>
inline std::array<U,N> parse_array(const std::vector<std::string>& toks, std::size_t& i) {
    std::array<U,N> a{};
    for (std::size_t k = 0; k < N; ++k) {
        if (i >= toks.size()) throw std::runtime_error("insufficient tokens for std::array");
        a[k] = parse_scalar<U>(toks[i++]);
    }
    return a;
}

/**
 * @brief Consumes tokens to produce a value of type `T`, advancing the index.
 */
template<class T>
inline auto take(const std::vector<std::string>& toks, std::size_t& i) {
    using Decayed = std::remove_cv_t<std::remove_reference_t<T>>;

    if constexpr (std::is_same_v<Decayed, std::string> || std::is_arithmetic_v<Decayed>) {
        if (i >= toks.size()) throw std::runtime_error("insufficient tokens for scalar/string");
        return parse_scalar<Decayed>(toks[i++]);
    } else if constexpr (is_std_array_v<Decayed>) {
        using U = typename Decayed::value_type;
        constexpr std::size_t N = std::tuple_size<Decayed>::value;
        return parse_array<U,N>(toks, i);
    } else {
        static_assert(!sizeof(T), "Unsupported parameter type in DSL invocation");
    }
}

/**
 * @brief Returns how many tokens are required to form type `T`.
 *        1 for scalars/strings, `N` for `std::array<*,N>`.
 */
template<class T>
constexpr std::size_t token_count_for() {
    using Decayed = std::remove_cv_t<std::remove_reference_t<T>>;

    if constexpr (is_std_array_v<Decayed>) {
        return std::tuple_size<Decayed>::value;
    } else {
        return 1;
    }
}

/** @brief Sum of token counts for a parameter pack. */
template<class... Us>
constexpr std::size_t token_count_for_pack() {
    return (token_count_for<Us>() + ... + 0);
}

/**
 * @brief Applies a callable `F` to tokens converted into the parameter types `Us...`.
 */
template<class... Us>
struct Invoker {
    template<class F>
    static void run(F&& f, const std::vector<std::string>& toks) {
        std::size_t i = 0;
        auto args = std::tuple<std::remove_cv_t<std::remove_reference_t<Us>>...>{ take<Us>(toks, i)... };
        std::apply(std::forward<F>(f), args);
    }
};

} // namespace detail
} // namespace dsl
} // namespace fem
