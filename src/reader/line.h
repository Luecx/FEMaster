#pragma once

#include <algorithm>
#include <cctype>
#include <ostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "../core/core.h"

namespace fem {
namespace reader {

/******************************************************************************
 * @brief Enum representing the type of line parsed from the input.
 ******************************************************************************/
enum LineType {
    COMMENT = -1,  ///< Line is a comment.
    EMPTY_LINE = 0, ///< Line is empty or whitespace.
    KEYWORD_LINE = 1, ///< Line contains a keyword with associated keys.
    DATA_LINE = 2, ///< Line contains data values.
    END_OF_FILE = 3, ///< End of file marker.
};

/******************************************************************************
 * @brief Represents a single line in the input file, handling comments,
 * keywords, and data.
 *
 * This struct can parse the line into relevant data, such as keyword commands
 * and associated keys, or raw data values. It supports access and querying
 * of the parsed information, and type-safe conversions of values.
 ******************************************************************************/
struct Line {
private:
    std::string m_line; ///< Raw line data.
    std::vector<std::string> m_values; ///< Parsed values for DATA_LINE.
    std::unordered_map<std::string, std::string> m_keys; ///< Parsed keys for KEYWORD_LINE.
    std::string m_command; ///< Command associated with a KEYWORD_LINE.
    LineType m_type = EMPTY_LINE; ///< The type of the line.

    //-------------------------------------------------------------------------
    // Private helper function
    //-------------------------------------------------------------------------
    /**
     * @brief Converts a string to a specified type T.
     *
     * @param s The string to convert.
     * @return The value of type T.
     */
    template <typename T>
    T convert_to(const std::string& s) const {
        T value;
        std::istringstream iss(s);
        iss >> value;
        logging::error(!iss.fail(), "Failed to convert value: " + s);
        return value;
    }

public:
    //-------------------------------------------------------------------------
    // Assignment operator to set the line content and parse its type.
    //-------------------------------------------------------------------------
    /**
     * @brief Assign a string to the Line object, parsing its type.
     *
     * @param line The string representing the line content.
     * @return Reference to the Line object after assignment.
     */
    Line& operator=(const std::string& line);

    //-------------------------------------------------------------------------
    // Getter functions
    //-------------------------------------------------------------------------
    /**
     * @brief Get the raw line content.
     *
     * @return Reference to the raw line string.
     */
    const std::string& line() const;

    /**
     * @brief Get the parsed values (for DATA_LINE type).
     *
     * @return Reference to the vector of values.
     */
    const std::vector<std::string>& values() const;

    /**
     * @brief Check if a specific key exists in the KEYWORD_LINE.
     *
     * @param key The key to check.
     * @return True if the key exists, false otherwise.
     */
    bool has_key(const std::string& key) const;

    /**
     * @brief Get the command name for KEYWORD_LINE.
     *
     * @return The command string.
     */
    const std::string& command() const;

    /**
     * @brief Get the type of the line (COMMENT, EMPTY_LINE, etc.).
     *
     * @return The LineType value.
     */
    LineType type() const;

    /**
     * @brief Check if the line is relevant (not a comment or empty line).
     *
     * @return True if the line is relevant, false otherwise.
     */
    bool ignorable() const {
        return m_type == COMMENT || m_type == EMPTY_LINE;
    }

    //-------------------------------------------------------------------------
    // Parsing functions for keys and values
    //-------------------------------------------------------------------------
    /**
     * @brief Get a required key value and throw an error if the key doesn't exist.
     *
     * @param key The key to lookup.
     * @return The value associated with the key.
     */
    template <typename T>
    T require(const std::string& key) const;

    /**
     * @brief Get a required value by checking multiple keys. Throws an error
     * if none of the keys exist.
     *
     * @param first_key The first key to lookup.
     * @param other_keys Additional keys to check.
     * @return The value associated with the first existing key.
     */
    template <typename T, typename... Args>
    T require(const std::string& first_key, const Args&... other_keys) const;

    /**
     * @brief Get a value from the keys or return a default value if the key doesn't exist.
     *
     * @param key The key to lookup.
     * @param default_value The value to return if the key is not found.
     * @return The value associated with the key, or the default value.
     */
    template <typename T>
    T parse(const std::string& key, const T& default_value) const;

    /**
     * @brief Get the number of values (for DATA_LINE type).
     *
     * @return The number of values.
     */
    size_t count_values() const;

    /**
     * @brief Get a specific value from the list of parsed values.
     *
     * @param index The index of the value in the list.
     * @param default_value The value to return if the index is out of bounds.
     * @return The value at the specified index, or the default value.
     */
    template <typename T>
    T get_value(size_t index, const T& default_value) const;

    //-------------------------------------------------------------------------
    // Setters
    //-------------------------------------------------------------------------
    /**
     * @brief Sets__OLD the end of file status to this line.
     */
    void eof();

    //-------------------------------------------------------------------------
    // Output operator
    //-------------------------------------------------------------------------
    /**
     * @brief Output the line information to an output stream.
     *
     * @param os The output stream.
     * @param line The Line object to output.
     * @return Reference to the output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Line& line) {
        os << "Line   : " << line.line() << "\n";
        os << "Type   : " << line.type() << "\n";

        if (line.type() == KEYWORD_LINE) {
            os << "Command: " << line.command() << "\n";
            os << "Keys   :\n";
            for (const auto& pair : line.m_keys) {
                os << std::right << std::setw(24) << pair.first << ": " << pair.second << "\n";
            }
        } else if (line.type() == DATA_LINE) {
            os << "Values :\n";
            for (size_t i = 0; i < line.count_values(); ++i) {
                os << std::right << std::setw(24) << i << ": " << line.get_value<std::string>(i, "") << "\n";
            }
        }

        return os;
    }
};

// Implementation of templated methods
template<typename T>
T Line::get_value(size_t index, const T& default_value) const {
    logging::error(m_type == DATA_LINE, "The 'get_value' function can only be used with DATA_LINE type.");
    if (index >= m_values.size()) {
        return default_value;
    }
    return convert_to<T>(m_values[index]);
}
template<typename T>
T Line::parse(const std::string& key, const T& default_value) const {
    logging::error(m_type == KEYWORD_LINE, "The 'parse' function can only be used with KEYWORD_LINE type.");
    auto it = m_keys.find(key);
    if (it == m_keys.end()) {
        return default_value;
    }
    return convert_to<T>(it->second);
}
template<typename T>
T Line::require(const std::string& key) const {
    logging::error(m_type == KEYWORD_LINE, "The 'require' function can only be used with KEYWORD_LINE type.");
    auto it = m_keys.find(key);
    logging::error(it != m_keys.end(), "Key does not exist: " + key);
    return convert_to<T>(it->second);
}

// Variadic require() function to accept multiple keys
template <typename T, typename... Args>
T Line::require(const std::string& first_key, const Args&... other_keys) const {
    logging::error(m_type == KEYWORD_LINE, "The 'require' function can only be used with KEYWORD_LINE type.");

    auto it = m_keys.find(first_key);

    if (it != m_keys.end()) {
        return convert_to<T>(it->second);  // Found the first key, return its value
    } else if constexpr (sizeof...(other_keys) > 0) {
        // Recursively check the next keys
        return require<T>(other_keys...);
    } else {
        // If none of the keys exist, throw an error
        logging::error(false, "None of the provided keys exist: " + first_key);
        return T();  // Return a default-constructed value (this line should never be reached)
    }
}

}
}    // namespace fem
