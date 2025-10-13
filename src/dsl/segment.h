/**
 * @file segment.h
 * @brief Declares `Segment`, a unit of data acquisition inside a variant.
 *
 * A `Segment` describes:
 *  - how many data lines may be consumed (`LineRange`),
 *  - how tokens are structured (`Pattern`),
 *  - and how parsed tokens are delivered to user code (`bind(...)`).
 *
 * The `.bind(F)` API deduces the lambda parameter types and converts tokens into
 * those types automatically. Supported parameter types are:
 *  - arithmetic scalars (`int`, `float`, `double`, `long`, `unsigned`, …)
 *  - `std::string`
 *  - `std::array<T,N>` (consumes `N` tokens of type `T`)
 *
 * Example:
 * @code
 * Segment::make()
 *   .range(LineRange{}.min(1).max(1))
 *   .pattern(Pattern::make().fixed<int,4>().fixed<double,1>())
 *   .bind([](std::array<int,4> ids, double k){
 *       // ...
 *   });
 * @endcode
 *
 * @see line_range.h
 * @see pattern.h
 * @see variant.h
 * @date 12.10.2025
 * @author
 *   Finn Eggers (project)
 */

#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <array>
#include <tuple>
#include <type_traits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <string>
#include <numeric>

#include "line_range.h"
#include "pattern.h"
#include "invoke.h"

namespace fem {
namespace dsl {

/**
 * @class Segment
 * @brief Unit of parsing: consumes a bounded number of lines, parses tokens via a pattern,
 *        and invokes a user callback with typed arguments.
 */
struct Segment {
    /** Line range constraint for this segment. */
    LineRange _range {};

    /** Token structure (sequence of pattern elements). */
    Pattern _pattern {};

    /** Type-erased invoker called with the flattened token list. */
    std::function<void(const std::vector<std::string>&)> _invoke;

    /**
     * @brief Factory for fluent construction.
     */
    static Segment make() { return {}; }

    /**
     * @brief Sets the line range for this segment.
     */
    Segment& range(LineRange r) {
        _range = r;
        return *this;
    }

    /**
     * @brief Sets the token pattern for this segment.
     */
    Segment& pattern(Pattern p) {
        _pattern = std::move(p);
        return *this;
    }

    /**
     * @brief Binds a user callable; parameter types are deduced from the callable signature.
     *
     * Supported parameter types:
     *  - arithmetic scalars (`int`, `double`, `float`, `long`, `unsigned`, …)
     *  - `std::string`
     *  - `std::array<T,N>` (consumes `N` tokens of type `T`)
     *
     * @tparam F Callable type (lambda, function object, function pointer, std::function)
     * @param fn Callable to be invoked with converted arguments.
     * @return Reference to `*this` for fluent chaining.
     */
    template<class F>
    Segment& bind(F&& fn) {
        using traits     = detail::function_traits<std::decay_t<F>>;
        using args_tuple = typename traits::args_tuple;
        constexpr std::size_t N = std::tuple_size<args_tuple>::value;
        return bind_impl(std::forward<F>(fn), std::make_index_sequence<N>{});
    }

private:
    template<class F, std::size_t... I>
    Segment& bind_impl(F&& fn, std::index_sequence<I...>) {
        using traits     = detail::function_traits<std::decay_t<F>>;
        using args_tuple = typename traits::args_tuple;
        using Inv        = detail::Invoker<std::tuple_element_t<I, args_tuple>...>;

        // Optional: sanity check against the pattern (can be relaxed if desired)
        const std::size_t need_lambda  = detail::token_count_for_pack<std::tuple_element_t<I, args_tuple>...>();
        const std::size_t need_pattern = _pattern.required_tokens();
        if (need_pattern != 0 && need_pattern != need_lambda) {
            std::ostringstream os;
            os << "Pattern token count (" << need_pattern
               << ") does not match lambda token count (" << need_lambda << ").";
            throw std::runtime_error(os.str());
        }

        _invoke = [fn_capture = std::forward<F>(fn)](const std::vector<std::string>& toks) mutable {
            Inv::run(fn_capture, toks);
        };
        return *this;
    }
};

} // namespace dsl
} // namespace fem
