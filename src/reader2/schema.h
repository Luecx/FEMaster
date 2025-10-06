// parser/schema.h
#pragma once
#include <functional>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "types.h"
#include "context.h"

namespace fem::reader2 {

struct Schema {
    // A handler tries to match a single DATA line (vector of strings)
    // Returns true if matched & executed
    using Handler = std::function<bool(Context&, const std::vector<std::string>&, const LineMeta&)>;

    struct AltDoc {
        std::string title;        // short label, e.g. "Node support (6 dof)"
        std::string description;  // free-form explanation
        std::vector<std::string> column_names;     // ["ID","UX","UY","UZ","RX","RY","RZ"]
        std::vector<std::string> column_meanings;  // same length; human help per column
        std::vector<std::string> column_units;     // optional units per column (may be empty strings)
    };

    std::vector<Handler> alternatives;
    std::vector<AltDoc>  alt_docs; // 1:1 with alternatives

    // --- Define directly in class (header-only) ---
    static Schema one_of(std::initializer_list<Schema> xs) {
        Schema s;
        for (auto& x : xs) {
            s.alternatives.insert(s.alternatives.end(), x.alternatives.begin(), x.alternatives.end());
            s.alt_docs.insert(s.alt_docs.end(), x.alt_docs.begin(), x.alt_docs.end());
        }
        return s;
    }

    template<class... Ts>
    struct Columns {
        bool na_nan = false;
        std::function<void(Context&, Ts..., const LineMeta&)> cb;

        // docs
        AltDoc doc;

        Columns& na_floats_as_nan(bool on=true){ na_nan = on; return *this; }
        Columns& bind(std::function<void(Context&, Ts..., const LineMeta&)> f){ cb = std::move(f); return *this; }

        // ---- Documentation helpers ----
        Columns& title(std::string t){ doc.title = std::move(t); return *this; }
        Columns& describe(std::string d){ doc.description = std::move(d); return *this; }
        Columns& cols_named(std::initializer_list<std::string> names){ doc.column_names.assign(names.begin(), names.end()); return *this; }
        Columns& cols_meaning(std::initializer_list<std::string> meanings){ doc.column_meanings.assign(meanings.begin(), meanings.end()); return *this; }
        Columns& cols_units(std::initializer_list<std::string> units){ doc.column_units.assign(units.begin(), units.end()); return *this; }
    };

    template<class... Ts>
    static Columns<Ts...> cols(){ return {}; }

    // Header-only implementation (no schema.cpp)
    template<class... Ts>
    static Schema make_from_columns(Columns<Ts...> c) {
        Schema s;

        // handler
        s.alternatives.push_back(
            [c](Context& cx, const std::vector<std::string>& vals, const LineMeta& meta) mutable -> bool {
                if (vals.size() != sizeof...(Ts)) return false;

                auto convert_cell = [&](auto tag, const std::string& str) {
                    using T = std::decay_t<decltype(tag)>;
                    if constexpr (std::is_same_v<T, std::string>) {
                        return str;
                    } else if constexpr (std::is_integral_v<T> && !std::is_same_v<T, bool>) {
                        if (str.empty()) throw std::runtime_error("Empty integer cell");
                        std::istringstream iss(str);
                        T v{}; iss >> v;
                        if (iss.fail()) throw std::runtime_error("Int conv fail: " + str);
                        return v;
                    } else if constexpr (std::is_floating_point_v<T>) {
                        if (str.empty()) {
                            return T(c.na_nan ? std::numeric_limits<T>::quiet_NaN() : T(0));
                        }
                        std::istringstream iss(str);
                        T v{}; iss >> v;
                        if (iss.fail()) throw std::runtime_error("Float conv fail: " + str);
                        return v;
                    } else {
                        static_assert(!sizeof(T*), "Unsupported type in Schema::Columns");
                    }
                };

                try {
                    size_t i = 0;
                    auto tup = std::tuple<Ts...>{ convert_cell(std::decay_t<Ts>{}, vals[i++])... };
                    if (c.cb) {
                        std::apply([&](auto&&... xs){
                            c.cb(cx, std::forward<decltype(xs)>(xs)..., meta);
                        }, tup);
                    }
                    return true;
                } catch(...) {
                    return false;
                }
            }
        );

        // docs
        s.alt_docs.push_back(std::move(c.doc));
        return s;
    }
};

// Utility for line ranges in CommandSpec
struct LineRange { size_t min_lines; size_t max_lines; };

} // namespace fem::reader2
