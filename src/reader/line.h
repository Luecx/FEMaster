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

enum LineType{
    COMMENT = -1,
    EMPTY_LINE = 0,
    KEYWORD_LINE = 1,
    DATA_LINE = 2,
    END_OF_FILE = 3,
};

bool relevant(LineType lt);

struct Line{
    private:
    // raw line data
    std::string m_line;

    // in case of a DATA_LINE, store the values
    std::vector<std::string> m_values;

    // in case of a KEYWORD_LINE; store the keys
    std::unordered_map<std::string, std::string> m_keys{};

    // store the command name in case of a keyword
    std::string m_command{};

    // store the type of the current line
    LineType m_type = EMPTY_LINE;

    public:

    // Private helper function to convert string to T
    template <typename T>
    T convert_to(const std::string& s) const {
        T value;
        std::istringstream iss(s);
        iss >> value;
        logging::error(!iss.fail(), "Failed to convert value: " + s);
        return value;
    }

    Line& operator=(const std::string& line);


    // Getter functions
    const std::string& line() const;

    const std::vector<std::string>& values() const;

    bool has_key(const std::string& key) const;

    const std::string& command() const;

    LineType type() const;

    // "require" function to get a value from the keys and throw an error if it does not exist.
    template <typename T>
    T require(const std::string& key) const;

    // Variadic require() function to accept multiple keys
    template <typename T, typename... Args>
    T require(const std::string& first_key, const Args&... other_keys) const;

    // "parse" function to get a value from the keys or return a default value if it does not exist.
    template <typename T>
    T parse(const std::string& key, const T& default_value) const;

    // "count_values" function to get the number of values
    size_t count_values() const;

    // "get_value" function to get a specific value from the list of values
    template <typename T>
    T get_value(size_t index, const T& default_value) const;

    // sets the end of file status to this line
    void eof();

    friend std::ostream& operator<<(std::ostream& os, const Line& line){
        os << "Line   : " << line.line() << "\n";
        os << "Type   : " << line.type() << "\n";

        if (line.type() == KEYWORD_LINE) {
            os << "Command: " << line.command() << "\n";
            os << "Keys   :\n";
            int i = 0;
            for (const auto& pair : line.m_keys) {
                os << std::right << std::setw(24) << pair.first << ": " << pair.second << "\n";
                ++i;
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
