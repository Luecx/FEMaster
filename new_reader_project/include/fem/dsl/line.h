/**
 * @file line.h
 * @brief Tokenizes raw input lines into keyword/data/comment records for the DSL.
 *
 * Responsibilities:
 *  - Trim whitespace, detect line type (comment/empty/keyword/data)
 *  - Normalize keyword lines: remove spaces, uppercase, drop leading '*'
 *  - Normalize data lines: remove spaces, uppercase, split by ','
 *  - Expose parsed command (`*FOO` → `"FOO"`) and key/value pairs for keyword lines
 *  - Expose vector of values (tokens) for data lines
 *
 * @note Normalization uppercases tokens. If you need case-sensitive strings, adapt the tokenizer.
 *
 * @see file.h
 * @see keys.h
 * @date 12.10.2025
 */

#pragma once
#include <ostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace fem {
namespace dsl {

/**
 * @brief Line classification used by the tokenizer.
 */
enum LineType {
    COMMENT = -1,
    EMPTY_LINE = 0,
    KEYWORD_LINE = 1,
    DATA_LINE = 2,
    END_OF_FILE = 3,
};

/**
 * @class Line
 * @brief Holds the parsed view of a single input line.
 *
 * A `Line` is re-parsed upon assignment from a raw string. Keyword and data lines are
 * normalized by removing whitespace and converting to uppercase, ensuring deterministic
 * downstream parsing. Keyword lines additionally strip the leading `*` and split into a
 * command name plus comma-separated key/value pairs.
 *
 * Data lines are split into comma-separated tokens and exposed via `values()`.
 *
 * @warning Because normalization uppercases tokens, case-sensitive string payloads will
 *          need a different tokenizer or an opt-out to preserve case.
 */
struct Line {
private:
    std::string _line;                                  ///< Normalized raw line buffer.
    std::vector<std::string> _values;                   ///< Tokens for data lines.
    std::unordered_map<std::string, std::string> _keys; ///< Keys for keyword lines.
    std::string _command;                               ///< Keyword command (without '*').
    LineType _type = EMPTY_LINE;                        ///< Classification.

    template <typename T>
    T convert_to(const std::string& s) const {
        T value{};
        std::istringstream iss(s);
        iss >> value;
        if (iss.fail()) throw std::runtime_error("convert_to failed for: " + s);
        return value;
    }

public:
    /**
     * @brief Assigns a raw line and (re)parses it.
     *
     * This function trims leading/trailing whitespace, classifies the line, and, for
     * keyword/data lines, produces a normalized, uppercase representation suitable for
     * consistent subsequent processing across the DSL.
     */
    Line& operator=(const std::string& line) {
        _line = line;
        _values.clear();
        _keys.clear();
        _command.clear();
        _type = END_OF_FILE;

        // trim left
        _line.erase(_line.begin(), std::find_if(_line.begin(), _line.end(),
                                                [](int ch){ return !std::isspace(ch); }));
        // trim right
        _line.erase(std::find_if(_line.rbegin(), _line.rend(),
                                 [](int ch){ return !std::isspace(ch); }).base(), _line.end());

        // classify
        if (_line.rfind("//",0)==0 || (!_line.empty() && (_line[0]=='#' || _line[0]=='!')) || (_line.rfind("**",0)==0)) {
            _type = COMMENT;
        } else if (_line.empty()) {
            _type = EMPTY_LINE;
        } else {
            if (_line[0]=='*') {
                _type = KEYWORD_LINE;
                std::string s = _line;
                // remove spaces, uppercase, strip '*'
                s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
                for (auto& ch : s)
                    ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                if (!s.empty() && s[0]=='*')
                    s.erase(0,1);

                std::stringstream ss(s);
                std::string item;
                if (std::getline(ss, _command, ',')) {
                    while (std::getline(ss, item, ',')) {
                        size_t pos = item.find('=');
                        if (pos != std::string::npos)
                            _keys[item.substr(0,pos)] = item.substr(pos+1);
                        else
                            _keys[item] = "";
                    }
                }
                _line = s; // keep normalized keyword line
            } else {
                _type = DATA_LINE;
                std::string s = _line;
                // remove spaces, uppercase
                s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
                for (auto& ch : s)
                    ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                std::stringstream ss(s);
                std::string item;
                while (std::getline(ss, item, ','))
                    _values.push_back(item);
                _line = s;
            }
        }
        return *this;
    }

    /**
     * @brief Returns the normalized line buffer.
     */
    const std::string& line() const { return _line; }

    /**
     * @brief Returns the data tokens (for `DATA_LINE`).
     */
    const std::vector<std::string>& values() const { return _values; }

    /**
     * @brief Returns true if the keyword line has a specific key.
     */
    bool has_key(const std::string& key) const { return _keys.find(key)!=_keys.end(); }

    /**
     * @brief Returns the command (for `KEYWORD_LINE`).
     */
    const std::string& command() const { return _command; }

    /**
     * @brief Returns the line type.
     */
    LineType type() const { return _type; }

    /**
     * @brief Whether this line can be skipped (comment/empty).
     */
    bool ignorable() const { return _type==COMMENT || _type==EMPTY_LINE; }

    /**
     * @brief Returns the number of data tokens (only for `DATA_LINE`).
     */
    std::size_t count_values() const {
        if (_type != DATA_LINE)
            throw std::runtime_error("count_values only for DATA_LINE");
        return _values.size();
    }

    /**
     * @brief Returns data token `index` parsed as `T`, or `def` on failure/out-of-range.
     *
     * @tparam T Target type.
     * @param index Zero-based token index.
     * @param def Default returned on conversion failure or out-of-range.
     */
    template<typename T>
    T get_value(std::size_t index, const T& def) const {
        if (_type != DATA_LINE)
            throw std::runtime_error("get_value only for DATA_LINE");
        if (index >= _values.size()) return def;
        std::istringstream ss(_values[index]);
        T out{};
        ss >> out;
        return ss.fail() ? def : out;
    }

    /**
     * @brief Appends this line’s normalized tokens to an output vector.
     *
     * Used when aggregating multiple consecutive `DATA_LINE`s within a multiline
     * pattern or segment definition. All tokens are already normalized
     * (uppercase, whitespace removed) by the parser.
     *
     * @param out Destination vector to which tokens will be appended.
     *
     * @throws std::runtime_error If this line is not of type `DATA_LINE`.
     */
    void append_values(std::vector<std::string>& out) const {
        if (_type != DATA_LINE)
            throw std::runtime_error("append_values: expected DATA_LINE");
        out.insert(out.end(), _values.begin(), _values.end());
    }

    /**
     * @brief Marks this line as end-of-file sentinel.
     */
    void eof() { _type = END_OF_FILE; }

    /**
     * @brief Stream output for diagnostics (shows type, command/keys or values).
     */
    friend std::ostream& operator<<(std::ostream& os, const Line& line) {
        os << "Line   : " << line.line() << "\n";
        os << "Type   : " << line.type() << "\n";
        if (line.type() == KEYWORD_LINE) {
            os << "Command: " << line.command() << "\n";
            os << "Keys   :\n";
            for (const auto& p : line._keys)
                os << std::setw(10) << p.first << ": " << p.second << "\n";
        } else if (line.type() == DATA_LINE) {
            os << "Values :\n";
            for (std::size_t i=0; i<line._values.size(); ++i)
                os << std::setw(10) << i << ": " << line._values[i] << "\n";
        }
        return os;
    }
};

} // namespace dsl
} // namespace fem
