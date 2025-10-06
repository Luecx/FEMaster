#pragma once
#include <iomanip>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "types.h"

namespace fem::reader2 {

/**
 * \brief Represents a fully parsed line from an input stream.
 */
struct Line {
public:
    /**
     * \brief Default constructor.
     */
    Line() = default;

    /**
     * \brief Parse the supplied raw string and populate the line metadata.
     * \param text Raw line content as read from file.
     * \return Reference to this line instance for chaining.
     */
    Line& assign(const std::string& text);

    /**
     * \brief Retrieve the trimmed raw text representation.
     */
    const std::string& raw() const {
        return _raw;
    }

    /**
     * \brief Retrieve the parsed data tokens for data lines.
     */
    const std::vector<std::string>& values() const {
        return _values;
    }

    /**
     * \brief Retrieve the parsed keyword key/value pairs.
     */
    const std::unordered_map<std::string, std::string>& keys() const {
        return _keys;
    }

    /**
     * \brief Command name for keyword lines.
     */
    const std::string& command() const {
        return _command;
    }

    /**
     * \brief Identify the parsed line classification.
     */
    LineType type() const {
        return _type;
    }

    /**
     * \brief Determine whether the line can be skipped during processing.
     */
    bool ignorable() const {
        return _type == LineType::COMMENT || _type == LineType::EMPTY;
    }

    /**
     * \brief Mark the line as representing end-of-file.
     */
    void mark_eof() {
        _type = LineType::EOF_MARK;
    }

private:
    std::string                                  _raw;                    ///< Trimmed raw string.
    std::vector<std::string>                     _values;                 ///< Data tokens for data lines.
    std::unordered_map<std::string, std::string> _keys;                   ///< Keyword key/value mapping.
    std::string                                  _command;                ///< Keyword command identifier.
    LineType                                     _type = LineType::EMPTY; ///< Current line type.
};

/**
 * \brief Stream the parsed line content for diagnostic purposes.
 */
std::ostream& operator<<(std::ostream& os, const Line& line);

} // namespace fem::reader2
