
#include "line.h"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace fem::reader2 {

static inline void ltrim(std::string& s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) { return !std::isspace(ch); }));
}

static inline void rtrim(std::string& s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) { return !std::isspace(ch); }).base(), s.end());
}

static inline void trim(std::string& s)
{
    ltrim(s);
    rtrim(s);
}

static inline void upcase(std::string& s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
}

Line& Line::assign(const std::string& s)
{
    _raw.clear();
    _values.clear();
    _keys.clear();
    _command.clear();
    _type = LineType::EMPTY;

    _raw = s;
    trim(_raw);

    if (_raw.empty()) {
        _type = LineType::EMPTY;
        return *this;
    }

    // comments: //, #, **, !  (case-insensitive)
    if (_raw.rfind("//", 0) == 0 || _raw.rfind("#", 0) == 0 || _raw.rfind("**", 0) == 0 || _raw.rfind("!", 0) == 0) {
        _type = LineType::COMMENT;
        return *this;
    }

    // keyword: starts with '*'
    if (!_raw.empty() && _raw[0] == '*') {
        std::string t = _raw.substr(1);
        trim(t);
        // split by commas into tokens, but keep empty tokens
        std::stringstream ss(t);
        std::string item;
        if (std::getline(ss, item, ',')) {
            // first token = command (uppercased, spaces trimmed)
            trim(item); upcase(item);
            _command = item;
            _type = LineType::KEYWORD;
            // rest = key=value OR bare keys
            while (std::getline(ss, item, ',')) {
                // do not uppercase values, only keys; but many workflows prefer all-caps
                // Here: uppercase whole item for consistency like your legacy reader
                trim(item);
                upcase(item);
                if (item.empty()) continue;
                auto pos = item.find('=');
                if (pos == std::string::npos) {
                    _keys[item] = "";
                } else {
                    std::string k = item.substr(0,pos);
                    std::string v = item.substr(pos+1);
                    trim(k); trim(v);
                    _keys[k] = v;
                }
            }
        } else {
            // only '*'
            _type = LineType::KEYWORD;
            _command.clear();
        }
        return *this;
    }

    // data line: uppercase for consistency, but DO NOT remove commas or collapse empties
    {
        std::string t = _raw;
        upcase(t); // matches your previous normalization
        std::stringstream ss(t);
        std::string item;
        while (std::getline(ss, item, ',')) {
            trim(item);
            // keep empty tokens as ""
            _values.push_back(item);
        }
        _type = LineType::DATA;
        return *this;
    }
}

std::ostream& operator<<(std::ostream& os, const Line& line)
{
    os << "Line   : " << line.raw() << "\n";
    os << "Type   : " << static_cast<int>(line.type()) << "\n";
    if (line.type() == LineType::KEYWORD) {
        os << "Command: " << line.command() << "\nKeys   :\n";
        for (const auto& kv : line.keys()) {
            os << std::setw(24) << kv.first << ": " << kv.second << "\n";
        }
    } else if (line.type() == LineType::DATA) {
        os << "Values :\n";
        for (size_t i = 0; i < line.values().size(); ++i) {
            os << std::setw(24) << i << ": " << line.values()[i] << "\n";
        }
    }
    return os;
}

} // namespace fem::reader2
