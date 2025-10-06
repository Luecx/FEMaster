#pragma once

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
#include "types.h"

namespace fem::reader2 {

/**
 * \brief Convert a string token into a strongly-typed value.
 */
template<class T>
inline T parse_cell(const std::string& token)
{
    if constexpr (std::is_same_v<T, std::string>) {
        return token;
    } else {
        std::istringstream stream(token);
        T                  value{};
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

template<class T>
struct TypeToValueKind;

template<>
struct TypeToValueKind<int> {
    static constexpr PatternElement::ValueKind value = PatternElement::ValueKind::Integer;
};

template<>
struct TypeToValueKind<double> {
    static constexpr PatternElement::ValueKind value = PatternElement::ValueKind::Floating;
};

template<>
struct TypeToValueKind<std::string> {
    static constexpr PatternElement::ValueKind value = PatternElement::ValueKind::Text;
};

template<class T>
constexpr PatternElement::ValueKind type_to_kind_v =
    TypeToValueKind<std::remove_cv_t<std::remove_reference_t<T>>>::value;

/**
 * \brief Fluent builder that captures the expected data layout for command segments.
 */
class Pattern {
public:
    /**
     * \brief Factory helper for fluent configuration.
     */
    static Pattern make()
    {
        return Pattern{};
    }

    /**
     * \brief Enable or disable multi-line consumption.
     */
    Pattern& allow_multiline(bool on = true)
    {
        _multiline = on;
        return *this;
    }

    /**
     * \brief Add a fixed-width field.
     */
    template<class T, size_t N = 1>
    Pattern& fixed(std::string name)
    {
        const auto type = N == 1 ? PatternElement::Type::Single : PatternElement::Type::Fixed;
        _elements.emplace_back(type,
                               type_to_kind_v<T>,
                               std::move(name),
                               N,
                               N,
                               _pendingDescription,
                               _pendingUnit);
        clear_pending_docs();
        return *this;
    }

    /**
     * \brief Add a variable-length tail field.
     */
    template<class T>
    Pattern& var(std::string base, size_t min_count, size_t max_count)
    {
        if (min_count > max_count) {
            throw std::runtime_error("var(): min_count > max_count");
        }
        _elements.emplace_back(PatternElement::Type::Variable,
                               type_to_kind_v<T>,
                               std::move(base),
                               min_count,
                               max_count,
                               _pendingDescription,
                               _pendingUnit);
        clear_pending_docs();
        return *this;
    }

    /**
     * \brief Attach documentation that describes the next added element.
     */
    Pattern& desc(std::string description)
    {
        _pendingDescription = std::move(description);
        return *this;
    }

    /**
     * \brief Attach unit text for the next added element.
     */
    Pattern& unit(std::string unit_text)
    {
        _pendingUnit = std::move(unit_text);
        return *this;
    }

    /**
     * \brief Provide a summary for documentation output.
     */
    Pattern& summary(std::string text)
    {
        _summary = std::move(text);
        return *this;
    }

    /**
     * \brief Provide additional notes for documentation output.
     */
    Pattern& notes(std::string text)
    {
        _notes = std::move(text);
        return *this;
    }

    /**
     * \brief Bind a callback that will receive the parsed tuple of values.
     */
    template<class... Ts, class F>
    Pattern& bind(F callback)
    {
        _invoker = [fun = std::move(callback)](Context& context, void* tuple_ptr, const LineMeta& meta) {
            auto& tuple = *static_cast<std::tuple<Ts...>*>(tuple_ptr);
            std::apply([&](auto&... xs) { fun(context, xs..., meta); }, tuple);
        };
        _converter = [](const std::vector<std::string>& tokens, const std::vector<size_t>& counts) {
            return convert_to_tuple<Ts...>(tokens, counts);
        };
        _arity = sizeof...(Ts);
        return *this;
    }

    /**
     * \brief Whether the pattern allows reading across multiple lines.
     */
    [[nodiscard]] bool multiline() const
    {
        return _multiline;
    }

    /**
     * \brief Validate structural invariants (only one variable element, at the end).
     */
    void validate() const
    {
        size_t var_seen = 0;
        for (size_t index = 0; index < _elements.size(); ++index) {
            if (_elements[index].type() == PatternElement::Type::Variable) {
                ++var_seen;
                if (index + 1 != _elements.size()) {
                    throw std::runtime_error("Pattern: variable element must be last");
                }
            }
        }
        if (var_seen > 1) {
            throw std::runtime_error("Pattern: only one variable element supported");
        }
    }

    /**
     * \brief Minimum tokens required to satisfy the pattern.
     */
    [[nodiscard]] size_t min_required_tokens() const
    {
        size_t count = 0;
        for (const auto& element : _elements) {
            count += element.min_count();
        }
        return count;
    }

    /**
     * \brief Maximum tokens the pattern may consume.
     */
    [[nodiscard]] size_t max_allowed_tokens() const
    {
        size_t count = 0;
        for (const auto& element : _elements) {
            count += element.max_count();
        }
        return count;
    }

    /**
     * \brief Retrieve column names for documentation up to the maximum extent.
     */
    [[nodiscard]] std::vector<std::string> column_names() const
    {
        std::vector<std::string> names;
        auto push_sequence = [&](const std::string& base, size_t total) {
            if (total == 1) {
                names.push_back(base);
                return;
            }
            for (size_t i = 1; i <= total; ++i) {
                names.push_back(base + std::to_string(i));
            }
        };
        for (const auto& element : _elements) {
            push_sequence(element.base_name(), element.max_count());
        }
        return names;
    }

    /**
     * \brief Summary text for documentation.
     */
    [[nodiscard]] const std::string& doc_summary() const
    {
        return _summary;
    }

    /**
     * \brief Notes text for documentation.
     */
    [[nodiscard]] const std::string& doc_notes() const
    {
        return _notes;
    }

    /**
     * \brief Access the configured elements.
     */
    [[nodiscard]] const std::vector<PatternElement>& elements() const
    {
        return _elements;
    }

    /**
     * \brief Convert tokens into the tuple used by the bound callback.
     */
    std::shared_ptr<void> convert(const std::vector<std::string>& tokens,
                                  const std::vector<size_t>& counts) const
    {
        if (!_converter) {
            throw std::runtime_error("Pattern not bound");
        }
        return _converter(tokens, counts);
    }

    /**
     * \brief Invoke the bound callback with the converted tuple.
     */
    void invoke(Context& context, void* tuple_ptr, const LineMeta& meta) const
    {
        if (!_invoker) {
            throw std::runtime_error("Pattern not bound");
        }
        _invoker(context, tuple_ptr, meta);
    }

    /**
     * \brief Number of bound tuple entries.
     */
    [[nodiscard]] size_t arity() const
    {
        return _arity;
    }

private:
    template<class... Ts>
    static std::shared_ptr<void> convert_to_tuple(const std::vector<std::string>& tokens,
                                                  const std::vector<size_t>& counts)
    {
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

    void clear_pending_docs()
    {
        _pendingDescription.clear();
        _pendingUnit.clear();
    }

    std::vector<PatternElement> _elements;
    bool                        _multiline = false;
    std::string                 _summary;
    std::string                 _notes;
    std::string                 _pendingDescription;
    std::string                 _pendingUnit;
    std::function<std::shared_ptr<void>(const std::vector<std::string>&, const std::vector<size_t>&)> _converter;
    std::function<void(Context&, void*, const LineMeta&)>                                              _invoker;
    size_t                                                                                _arity = 0;
};

} // namespace fem::reader2
