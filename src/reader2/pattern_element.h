#pragma once
/**
 * @file pattern_element.h
 * @brief Describes the semantics of a single pattern element.
 */

#include <cstddef>
#include <string>
#include <utility>

namespace fem::reader2 {

/**
 * @brief Describes one slot in a pattern (single/fixed/variable, kind, counts, docs).
 */
class PatternElement {
public:
    enum class Type {
        Single,
        Fixed,
        Variable
    };

    enum class ValueKind {
        Integer,
        Floating,
        Text
    };

    /**
     * @brief Construct an element.
     * @param type        Element class (Single/Fixed/Variable)
     * @param kind        Value kind (Integer/Floating/Text)
     * @param base_name   Base name used in docs/binding
     * @param min_count   Minimum tokens consumed
     * @param max_count   Maximum tokens consumed
     * @param description Human-readable description
     */
    PatternElement(Type type,
                   ValueKind kind,
                   std::string base_name,
                   size_t min_count,
                   size_t max_count,
                   std::string description)
        : _type(type)
        , _kind(kind)
        , _baseName(std::move(base_name))
        , _minCount(min_count)
        , _maxCount(max_count)
        , _description(std::move(description)) {}

    /** @brief Identify the element class. */
    [[nodiscard]] Type type() const { return _type; }
    /** @brief Retrieve the value kind. */
    [[nodiscard]] ValueKind kind() const { return _kind; }
    /** @brief Fetch the base name. */
    [[nodiscard]] const std::string& base_name() const { return _baseName; }
    /** @brief Minimum number of tokens. */
    [[nodiscard]] size_t min_count() const { return _minCount; }
    /** @brief Maximum number of tokens. */
    [[nodiscard]] size_t max_count() const { return _maxCount; }
    /** @brief Description text for docs. */
    [[nodiscard]] const std::string& description() const { return _description; }

    /** @brief Set/override description text. */
    void set_description(std::string d) { _description = std::move(d); }

    /** @brief Effective consumed tokens given extra variable tail. */
    size_t effective_count(size_t extras) const {
        if (_type == Type::Variable) return _minCount + extras;
        return _maxCount;
    }

private:
    Type        _type;
    ValueKind   _kind;
    std::string _baseName;
    size_t      _minCount;
    size_t      _maxCount;
    std::string _description;
};

} // namespace fem::reader2
