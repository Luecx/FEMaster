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

namespace fem { namespace reader2 {

// ----------------------------- detail helpers (namespace scope) -----------------------------
namespace detail {

// ---------- function_traits: deduce lambda/function parameter types ----------
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

// generic callable (lambda, functor): use operator() if present
template<typename F>
struct function_traits : function_traits<decltype(&F::operator())> {};

// ---------- type checks ----------
template<class T> struct is_std_array : std::false_type {};
template<class U, std::size_t N> struct is_std_array<std::array<U,N>> : std::true_type {};
template<class T> inline constexpr bool is_std_array_v = is_std_array<T>::value;

// ---------- parsers ----------
template<class T>
inline T parse_arithmetic(const std::string& s){
    static_assert(std::is_arithmetic_v<T>, "parse_arithmetic requires arithmetic type");
    std::istringstream ss(s);
    T v{};
    ss >> v;
    if (ss.fail()) throw std::runtime_error(std::string("parse_arithmetic failed for token: ") + s);
    return v;
}

template<class T>
inline T parse_scalar(const std::string& s){
    if constexpr (std::is_same_v<T, std::string>) {
        // Achtung: Deine aktuelle Line-Implementierung uppercaset Values.
        // Falls Case-sensitiv nötig ist, dort anpassen.
        return s;
    } else if constexpr (std::is_arithmetic_v<T>) {
        return parse_arithmetic<T>(s);
    } else {
        static_assert(!sizeof(T), "Unsupported scalar type in Segment::bind");
    }
}

template<class U, std::size_t N>
inline std::array<U,N> parse_array(const std::vector<std::string>& toks, std::size_t& i){
    std::array<U,N> a{};
    for (std::size_t k=0; k<N; ++k){
        if (i >= toks.size()) throw std::runtime_error("insufficient tokens for std::array");
        a[k] = parse_scalar<U>(toks[i++]);
    }
    return a;
}

template<class T>
inline T take(const std::vector<std::string>& toks, std::size_t& i){
    if constexpr (std::is_same_v<T, std::string> || std::is_arithmetic_v<T>){
        if (i >= toks.size()) throw std::runtime_error("insufficient tokens for scalar/string");
        return parse_scalar<T>(toks[i++]);
    } else if constexpr (is_std_array_v<T>) {
        using U = typename T::value_type;
        constexpr std::size_t N = std::tuple_size<T>::value;
        return parse_array<U,N>(toks, i);
    } else {
        static_assert(!sizeof(T), "Unsupported parameter type in Segment::bind");
    }
}

// wie viele Tokens werden für Typ T benötigt?
template<class T>
constexpr std::size_t token_count_for(){
    if constexpr (is_std_array_v<T>) {
        return std::tuple_size<T>::value;
    } else {
        return 1;
    }
}

template<class... Us>
constexpr std::size_t token_count_for_pack(){
    return (token_count_for<Us>() + ... + 0);
}

template<class... Us>
struct Invoker {
    template<class F>
    static void run(F&& f, const std::vector<std::string>& toks){
        std::size_t i = 0;
        auto args = std::tuple<Us...>{ take<Us>(toks, i)... };
        std::apply(std::forward<F>(f), args);
    }
};

} // namespace detail
// -------------------------------------------------------------------------------------------

inline std::size_t pattern_required_tokens_inline(const Pattern& p){
    std::size_t n=0; for (auto& e : p.elems) n += e->count(); return n;
}

struct Segment {
    LineRange range_{};
    Pattern   pattern_{};

    // Tokens -> typed user-callback
    std::function<void(const std::vector<std::string>&)> invoke_;

    static Segment make(){ return {}; }

    Segment& range(LineRange r){ range_ = r; return *this; }
    Segment& pattern(Pattern p){ pattern_ = std::move(p); return *this; }

    // Accepts any callable F; argument types are deduced from F's signature.
    // Examples:
    //   .bind([](double E,int a,int b,int c){...})
    //   .bind([](std::array<int,4> ids,double k){...})
    //   .bind([](std::string name,float alpha){...})
    template<class F>
    Segment& bind(F&& fn){
        using traits     = detail::function_traits<std::decay_t<F>>;
        using args_tuple = typename traits::args_tuple;

        // Expand tuple to types...
        constexpr std::size_t N = std::tuple_size<args_tuple>::value;
        return bind_impl(std::forward<F>(fn), std::make_index_sequence<N>{});
    }

    private:
    template<class F, std::size_t... I>
    Segment& bind_impl(F&& fn, std::index_sequence<I...>){
        using traits     = detail::function_traits<std::decay_t<F>>;
        using args_tuple = typename traits::args_tuple;
        using Inv        = detail::Invoker<std::tuple_element_t<I, args_tuple>...>;

        // Optional: Konsistenzcheck Pattern vs. Lambda-Parameter
        const std::size_t need_lambda  = detail::token_count_for_pack<std::tuple_element_t<I, args_tuple>...>();
        const std::size_t need_pattern = pattern_required_tokens_inline(pattern_);
        if (need_pattern != 0 && need_pattern != need_lambda){
            // Hinweis: Manche Users wollen bewusst abweichen (z.B. separate Segmente).
            // Hier helfen klare Fehlermeldungen.
            std::ostringstream os;
            os << "Pattern token count (" << need_pattern
               << ") does not match lambda token count (" << need_lambda << ").";
            throw std::runtime_error(os.str());
        }

        // Store the erasured invoker
        invoke_ = [fn_capture = std::forward<F>(fn)](const std::vector<std::string>& toks) mutable {
            Inv::run(fn_capture, toks);
        };
        return *this;
    }
};

}} // namespace fem::reader2
