#pragma once

#include <cstddef>
#include <string>
#include <utility>

namespace fem::reader2 {

/**
 * @brief Describes the semantics of a single pattern element.
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
     * @brief Construct an element describing a particular slot in the pattern.
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

    /** @brief Retrieve the value kind captured by the element. */
    [[nodiscard]] ValueKind kind() const { return _kind; }

    /** @brief Fetch the base name used for documentation and binding. */
    [[nodiscard]] const std::string& base_name() const { return _baseName; }

    /** @brief Minimum number of tokens the element expects. */
    [[nodiscard]] size_t min_count() const { return _minCount; }

    /** @brief Maximum number of tokens the element accepts. */
    [[nodiscard]] size_t max_count() const { return _maxCount; }

    /** @brief Human-readable description for documentation output. */
    [[nodiscard]] const std::string& description() const { return _description; }

    /** @brief Overwrite the human-readable description (used by Pattern::desc). */
    void set_description(std::string description) { _description = std::move(description); }

    /** @brief Compute how many tokens the element actually consumes with extras. */
    size_t effective_count(size_t extras) const {
        if (_type == Type::Variable) {
            return _minCount + extras;
        }
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
