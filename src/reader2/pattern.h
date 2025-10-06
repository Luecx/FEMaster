#pragma once
/**
 * @file pattern.h
 * @brief Fluent builder capturing expected data layout for command segments.
 */

#include <array>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "context.h"
#include "pattern_element.h"
#include "keyword.h"   // for Keyword::Map
#include "types.h"

namespace fem::reader2 {

/** @brief Convert a string token into a strongly-typed value. */
template<class T>
inline T parse_cell(const std::string& token) {
    if constexpr (std::is_same_v<T, std::string>) {
        return token;
    } else {
        std::istringstream stream(token);
        T value{};
        stream >> value;
        if (stream.fail()) {
            throw std::runtime_error("Parse failed for '" + token + "'");
        }
        return value;
    }
}

template<class T>
struct is_std_array : std::false_type {};
template<class U, size_t N>
struct is_std_array<std::array<U, N>> : std::true_type {};

template<class T>
struct is_std_vector : std::false_type {};
template<class U>
struct is_std_vector<std::vector<U>> : std::true_type {};

/**
 * @brief Map C++ types to PatternElement::ValueKind.
 * Supports all integral types, all floating-point types, and std::string.
 */
template<class T>
struct TypeToValueKind {
    using U = std::remove_cv_t<std::remove_reference_t<T>>;
    static_assert(std::is_integral_v<U> || std::is_floating_point_v<U> || std::is_same_v<U, std::string>,
                  "Unsupported type in Pattern::fixed/var (expected integral, floating-point, or std::string).");

    static constexpr PatternElement::ValueKind value =
        std::is_same_v<U, std::string> ? PatternElement::ValueKind::Text
      : (std::is_floating_point_v<U>   ? PatternElement::ValueKind::Floating
                                       : PatternElement::ValueKind::Integer);
};

template<class T>
constexpr PatternElement::ValueKind type_to_kind_v =
    TypeToValueKind<std::remove_cv_t<std::remove_reference_t<T>>>::value;

/**
 * @brief Fluent builder that captures the expected data layout for command segments.
 */
class Pattern {
public:
    /** @brief Factory helper for fluent configuration. */
    static Pattern make();

    /** @brief Enable or disable multi-line consumption. */
    Pattern& allow_multiline(bool on = true);

    /** @brief Add a fixed-width field. */
    template<class T, size_t N = 1>
    Pattern& fixed(std::string name) {
        const auto type = N == 1 ? PatternElement::Type::Single : PatternElement::Type::Fixed;
        _elements.emplace_back(type,
                               type_to_kind_v<T>,
                               std::move(name),
                               N,
                               N,
                               std::string{} /* desc set via desc() */);
        return *this;
    }

    /** @brief Add a variable-length tail field. */
    template<class T>
    Pattern& var(std::string base, size_t min_count, size_t max_count) {
        if (min_count > max_count) {
            throw std::runtime_error("var(): min_count > max_count");
        }
        _elements.emplace_back(PatternElement::Type::Variable,
                               type_to_kind_v<T>,
                               std::move(base),
                               min_count,
                               max_count,
                               std::string{} /* desc set via desc() */);
        return *this;
    }

    /**
     * @brief Attach documentation to the most recently added element.
     * @throws std::runtime_error if called before any element was added.
     */
    Pattern& desc(std::string description);

    /** @brief Provide a summary for documentation output. */
    Pattern& summary(std::string text);
    /** @brief Provide additional notes for documentation output. */
    Pattern& notes(std::string text);

    /**
     * @brief Bind a callback that receives the current keyword key→value map and parsed values.
     * Signature: void(Context&, const Keyword::Map&, Ts..., const LineMeta&)
     */
    template<class... Ts, class F>
    Pattern& bind_kv(F callback) {
        _invoker_with_kv = [fun = std::move(callback)](Context& context, const Keyword::Map& kv, void* tuple_ptr, const LineMeta& meta) {
            auto& tuple = *static_cast<std::tuple<Ts...>*>(tuple_ptr);
            std::apply([&](auto&... xs) { fun(context, kv, xs..., meta); }, tuple);
        };
        _converter = [](const std::vector<std::string>& tokens, const std::vector<size_t>& counts) {
            return convert_to_tuple<Ts...>(tokens, counts);
        };
        _arity = sizeof...(Ts);
        return *this;
    }

    /** @brief Whether the pattern allows reading across multiple lines. */
    [[nodiscard]] bool multiline() const;

    /** @brief Validate structural invariants (only one variable element, at the end). */
    void validate() const;

    /** @brief Minimum tokens required to satisfy the pattern. */
    [[nodiscard]] size_t min_required_tokens() const;

    /** @brief Maximum tokens the pattern may consume. */
    [[nodiscard]] size_t max_allowed_tokens() const;

    /** @brief Retrieve column names for documentation up to the maximum extent. */
    [[nodiscard]] std::vector<std::string> column_names() const;

    /** @brief Summary text for documentation. */
    [[nodiscard]] const std::string& doc_summary() const;
    /** @brief Notes text for documentation. */
    [[nodiscard]] const std::string& doc_notes() const;

    /** @brief Access the configured elements. */
    [[nodiscard]] const std::vector<PatternElement>& elements() const;

    /** @brief Convert tokens into the tuple used by the bound callback. */
    std::shared_ptr<void> convert(const std::vector<std::string>& tokens,
                                  const std::vector<size_t>& counts) const;

    /**
     * @brief Invoke the bound callback with the converted tuple (always passes Keyword::Map).
     */
    void invoke(Context& context, const Keyword& kw, void* tuple_ptr, const LineMeta& meta) const;

    /** @brief Number of bound tuple entries. */
    [[nodiscard]] size_t arity() const;

private:
    template<class... Ts>
    static std::shared_ptr<void> convert_to_tuple(const std::vector<std::string>& tokens,
                                                  const std::vector<size_t>& counts) {
        auto tuple = std::make_shared<std::tuple<Ts...>>();
        size_t token_index = 0;
        size_t count_index = 0;

        auto consume_vector = [&](auto* out_vector, size_t quantity) {
            using ValueType = typename std::decay_t<decltype(*out_vector)>::value_type;
            out_vector->clear();
            out_vector->reserve(quantity);
            for (size_t j = 0; j < quantity; ++j) {
                out_vector->push_back(parse_cell<ValueType>(tokens.at(token_index++)));
            }
        };

        ([&] {
            using Type = Ts;
            if constexpr (is_std_array<Type>::value) {
                constexpr size_t N = std::tuple_size_v<Type>;
                for (size_t j = 0; j < N; ++j) {
                    std::get<Type>(*tuple)[j] = parse_cell<typename Type::value_type>(tokens.at(token_index++));
                }
                ++count_index;
            } else if constexpr (is_std_vector<Type>::value) {
                size_t quantity = counts.at(count_index++);
                consume_vector(&std::get<Type>(*tuple), quantity);
            } else {
                std::get<Type>(*tuple) = parse_cell<Type>(tokens.at(token_index++));
                ++count_index;
            }
        }(), ...);

        return tuple;
    }

    std::vector<PatternElement> _elements;
    bool                        _multiline = false;
    std::string                 _summary;
    std::string                 _notes;

    std::function<void(Context&, const Keyword::Map&, void*, const LineMeta&)> _invoker_with_kv;
    std::function<std::shared_ptr<void>(const std::vector<std::string>&,
                                        const std::vector<size_t>&)>           _converter;
    size_t                                                                        _arity = 0;
};

} // namespace fem::reader2
