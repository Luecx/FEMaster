#pragma once

#include <algorithm>
#include <cctype>
#include <ostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace fem {
namespace reader {

enum LineType{
    COMMENT = -1,
    EMPTY_LINE = 0,
    KEYWORD_LINE = 1,
    DATA_LINE = 2,
    END_OF_FILE = 3,
};

bool relevant(LineType lt){
    return lt == KEYWORD_LINE || lt == DATA_LINE;
}

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
        log_error(!iss.fail(), "Failed to convert value: " + s);
        return value;
    }

    Line& operator=(const std::string& line){
        // reset fields
        m_line = line;
        m_values.clear();
        m_keys.clear();
        m_command.clear();
        m_type = END_OF_FILE;

        // remove leading spaces
        m_line.erase(m_line.begin(), std::find_if(m_line.begin(), m_line.end(), [](int ch) {
                         return !std::isspace(ch);
                     }));

        // check if the m_line is a comment --> starts with //, #, **, !
        if (m_line.substr(0, 2) == "//" || m_line[0] == '#' || m_line.substr(0, 2) == "**" || m_line[0] == '!') {
            m_type = COMMENT;
        } else if (m_line.empty()) {
            m_type = EMPTY_LINE;
        } else {
            // check if it is a keyword line
            if (m_line[0] == '*') {
                m_type = KEYWORD_LINE;
                m_line.erase(std::remove(m_line.begin(), m_line.end(), ' '), m_line.end());
                m_line.erase(0, 1);  // remove the leading '*'
                std::transform(m_line.begin(), m_line.end(), m_line.begin(), ::toupper);  // convert to upper case
                std::stringstream ss(m_line);
                std::string item;
                if (std::getline(ss, m_command, ',')) {  // get the command name
                    while (std::getline(ss, item, ',')) {
                        size_t pos = item.find("=");
                        if (pos != std::string::npos) {
                            m_keys[item.substr(0, pos)] = item.substr(pos + 1);
                        }
                    }
                }
            } else {
                m_type = DATA_LINE;
                m_line = std::regex_replace(m_line, std::regex(" +"), " ");  // replace multiple spaces with a single space
                std::stringstream ss(m_line);
                std::string item;
                while (std::getline(ss, item, ' ')) {
                    m_values.push_back(item);
                }
            }
        }

        return *this;
    }

    // Getter functions
    const std::string& line() const {
        return m_line;
    }

    const std::vector<std::string>& values() const {
        return m_values;
    }

    bool has_key(const std::string& key) const {
        return m_keys.find(key) != m_keys.end();
    }

    const std::string& command() const {
        return m_command;
    }

    LineType type() const {
        return m_type;
    }

    // "require" function to get a value from the keys and throw an error if it does not exist.
    template <typename T>
    T require(const std::string& key) const {
        log_error(m_type == KEYWORD_LINE, "The 'require' function can only be used with KEYWORD_LINE type.");
        auto it = m_keys.find(key);
        log_error(it != m_keys.end(), "Key does not exist: " + key);
        return convert_to<T>(it->second);
    }

    // "parse" function to get a value from the keys or return a default value if it does not exist.
    template <typename T>
    T parse(const std::string& key, const T& default_value) const {
        log_error(m_type == KEYWORD_LINE, "The 'parse' function can only be used with KEYWORD_LINE type.");
        auto it = m_keys.find(key);
        if (it == m_keys.end()) {
            return default_value;
        }
        return convert_to<T>(it->second);
    }

    // "count_values" function to get the number of values
    size_t count_values() const {
        log_error(m_type == DATA_LINE, "The 'count_values' function can only be used with DATA_LINE type.");
        return m_values.size();
    }

    // "get_value" function to get a specific value from the list of values
    template <typename T>
    T get_value(size_t index, const T& default_value) const {
        log_error(m_type == DATA_LINE, "The 'get_value' function can only be used with DATA_LINE type.");
        if (index >= m_values.size()) {
            return default_value;
        }
        return convert_to<T>(m_values[index]);
    }

    // sets the end of file status to this line
    void eof(){
        this->m_type = END_OF_FILE;
    }

    friend std::ostream& operator<<(std::ostream& os, const Line& line) {
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

}
}    // namespace fem
