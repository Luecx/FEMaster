
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

namespace fem { namespace reader {

enum LineType {
    COMMENT = -1,
    EMPTY_LINE = 0,
    KEYWORD_LINE = 1,
    DATA_LINE = 2,
    END_OF_FILE = 3,
};

struct Line {
private:
    std::string m_line;
    std::vector<std::string> m_values;
    std::unordered_map<std::string, std::string> m_keys;
    std::string m_command;
    LineType m_type = EMPTY_LINE;

    template <typename T>
    T convert_to(const std::string& s) const {
        T value{};
        std::istringstream iss(s);
        iss >> value;
        if (iss.fail()) { throw std::runtime_error("convert_to failed for: " + s); }
        return value;
    }

public:
    Line& operator=(const std::string& line) {
        m_line = line; m_values.clear(); m_keys.clear(); m_command.clear(); m_type = END_OF_FILE;
        m_line.erase(m_line.begin(), std::find_if(m_line.begin(), m_line.end(), [](int ch){ return !std::isspace(ch);}));
        m_line.erase(std::find_if(m_line.rbegin(), m_line.rend(), [](int ch){ return !std::isspace(ch);} ).base(), m_line.end());
        if (m_line.rfind("//",0)==0 || (!m_line.empty() && (m_line[0]=='#' || m_line[0]=='!')) || (m_line.rfind("**",0)==0)) {
            m_type = COMMENT;
        } else if (m_line.empty()) {
            m_type = EMPTY_LINE;
        } else {
            if (m_line[0]=='*') {
                m_type = KEYWORD_LINE;
                std::string s = m_line;
                // remove spaces, uppercase
                s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
                for (auto& ch : s) ch = std::toupper(ch);
                if (!s.empty() && s[0]=='*') s.erase(0,1);
                std::stringstream ss(s);
                std::string item;
                if (std::getline(ss, m_command, ',')) {
                    while (std::getline(ss, item, ',')) {
                        size_t pos = item.find("=");
                        if (pos != std::string::npos) m_keys[item.substr(0,pos)] = item.substr(pos+1);
                        else m_keys[item] = "";
                    }
                }
                m_line = s; // keep normalized keyword string for re-parsing if needed
            } else {
                m_type = DATA_LINE;
                std::string s = m_line;
                // remove spaces, uppercase
                s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
                for (auto& ch : s) ch = std::toupper(ch);
                std::stringstream ss(s);
                std::string item;
                while (std::getline(ss, item, ',')) m_values.push_back(item);
                m_line = s;
            }
        }
        return *this;
    }

    const std::string& line() const { return m_line; }
    const std::vector<std::string>& values() const { return m_values; }
    bool has_key(const std::string& key) const { return m_keys.find(key)!=m_keys.end(); }
    const std::string& command() const { return m_command; }
    LineType type() const { return m_type; }
    bool ignorable() const { return m_type==COMMENT || m_type==EMPTY_LINE; }

    template<typename T>
    T require(const std::string& key) const {
        auto it = m_keys.find(key);
        if (it==m_keys.end()) throw std::runtime_error("Missing key: " + key);
        std::istringstream ss(it->second);
        T out{}; ss>>out; if (ss.fail()) throw std::runtime_error("Bad key value for: " + key);
        return out;
    }

    template<typename T>
    T parse(const std::string& key, const T& def) const {
        auto it = m_keys.find(key); if (it==m_keys.end()) return def;
        std::istringstream ss(it->second); T out{}; ss>>out; return ss.fail()?def:out;
    }

    std::size_t count_values() const {
        if (m_type != DATA_LINE) throw std::runtime_error("count_values only for DATA_LINE");
        return m_values.size();
    }

    template<typename T>
    T get_value(std::size_t index, const T& def) const {
        if (m_type != DATA_LINE) throw std::runtime_error("get_value only for DATA_LINE");
        if (index >= m_values.size()) return def;
        std::istringstream ss(m_values[index]); T out{}; ss>>out; return ss.fail()?def:out;
    }

    void eof() { m_type = END_OF_FILE; }

    friend std::ostream& operator<<(std::ostream& os, const Line& line) {
        os << "Line   : " << line.line() << "\n";
        os << "Type   : " << line.type() << "\n";
        if (line.type() == KEYWORD_LINE) {
            os << "Command: " << line.command() << "\n";
            os << "Keys   :\n";
            for (const auto& p : line.m_keys) os << std::setw(10) << p.first << ": " << p.second << "\n";
        } else if (line.type() == DATA_LINE) {
            os << "Values :\n";
            for (std::size_t i=0;i<line.m_values.size();++i) os << std::setw(10) << i << ": " << line.m_values[i] << "\n";
        }
        return os;
    }
};

}} // namespace fem::reader
