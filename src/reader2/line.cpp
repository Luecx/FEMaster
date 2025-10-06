#include "line.h"
#include <algorithm>
#include <cctype>
#include <regex>
#include <sstream>

namespace fem::reader2 {

static inline void ltrim(std::string& s){ s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch){return !std::isspace(ch);})); }
static inline void rtrim(std::string& s){ s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch){return !std::isspace(ch);}).base(), s.end()); }
static inline void trim(std::string& s){ ltrim(s); rtrim(s); }
static inline void upcase(std::string& s){ std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){return std::toupper(c);}); }

Line& Line::assign(const std::string& s){
    m_raw.clear(); m_vals.clear(); m_keys.clear(); m_cmd.clear(); m_type = LineType::EMPTY;

    m_raw = s;
    trim(m_raw);

    if (m_raw.empty()) { m_type = LineType::EMPTY; return *this; }

    // comments: //, #, **, !  (case-insensitive)
    if (m_raw.rfind("//",0)==0 || m_raw.rfind("#",0)==0 || m_raw.rfind("**",0)==0 || m_raw.rfind("!",0)==0){
        m_type = LineType::COMMENT; return *this;
    }

    // keyword: starts with '*'
    if (!m_raw.empty() && m_raw[0]=='*'){
        std::string t = m_raw.substr(1);
        trim(t);
        // split by commas into tokens, but keep empty tokens
        std::stringstream ss(t);
        std::string item;
        if (std::getline(ss, item, ',')){
            // first token = command (uppercased, spaces trimmed)
            trim(item); upcase(item);
            m_cmd = item;
            m_type = LineType::KEYWORD;
            // rest = key=value OR bare keys
            while (std::getline(ss, item, ',')){
                // do not uppercase values, only keys; but many workflows prefer all-caps
                // Here: uppercase whole item for consistency like your legacy reader
                trim(item); upcase(item);
                if (item.empty()) continue;
                auto pos = item.find('=');
                if (pos==std::string::npos) m_keys[item] = ""; else {
                    std::string k = item.substr(0,pos);
                    std::string v = item.substr(pos+1);
                    trim(k); trim(v);
                    m_keys[k] = v;
                }
            }
        } else {
            // only '*'
            m_type = LineType::KEYWORD; m_cmd.clear();
        }
        return *this;
    }

    // data line: uppercase for consistency, but DO NOT remove commas or collapse empties
    {
        std::string t = m_raw;
        upcase(t); // matches your previous normalization
        std::stringstream ss(t);
        std::string item;
        while (std::getline(ss, item, ',')){
            trim(item);
            // keep empty tokens as ""
            m_vals.push_back(item);
        }
        m_type = LineType::DATA;
        return *this;
    }
}

std::ostream& operator<<(std::ostream& os, const Line& L){
    os << "Line   : " << L.raw() << "\n";
    os << "Type   : " << static_cast<int>(L.type()) << "\n";
    if (L.type()==LineType::KEYWORD){
        os << "Command: " << L.command() << "\nKeys   :\n";
        for (auto& kv : L.keys()) os << std::setw(24) << kv.first << ": " << kv.second << "\n";
    } else if (L.type()==LineType::DATA){
        os << "Values :\n";
        for (size_t i=0;i<L.values().size();++i) os << std::setw(24) << i << ": " << L.values()[i] << "\n";
    }
    return os;
}

} // namespace fem::reader2
