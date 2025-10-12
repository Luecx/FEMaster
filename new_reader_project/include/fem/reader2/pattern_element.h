
#pragma once
#include <array>
#include <string>
#include <memory>

namespace fem { namespace reader2 {
struct PatternElementBase {
    virtual ~PatternElementBase() = default;
    virtual std::size_t count() const = 0;
};
template<class T, std::size_t N>
struct Fixed : PatternElementBase {
    std::string description;
    bool allow_empty = false;
    bool allow_missing = false;
    T default_single{};
    std::array<T,N> default_array{};
    Fixed& desc(std::string d){ description = std::move(d); return *this; }
    Fixed& on_empty(const T& v){ allow_empty=true; default_single=v; default_array.fill(v); return *this; }
    Fixed& on_missing(const T& v){ allow_missing=true; default_single=v; default_array.fill(v); return *this; }
    std::size_t count() const override { return N; }
};
}} // ns
