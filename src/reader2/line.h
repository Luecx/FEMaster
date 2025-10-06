#pragma once
#include <string>
#include <unordered_map>
#include <vector>
#include <ostream>
#include <iomanip>
#include "types.h"

namespace fem::reader2 {

struct Line {
private:
    std::string m_raw;                // raw, trimmed, no surrounding spaces
    std::vector<std::string> m_vals;  // for DATA
    std::unordered_map<std::string,std::string> m_keys; // for KEYWORD
    std::string m_cmd;                // for KEYWORD
    LineType m_type = LineType::EMPTY;

public:
    Line() = default;

    Line& assign(const std::string& s);

    const std::string& raw() const { return m_raw; }
    const std::vector<std::string>& values() const { return m_vals; }
    const std::unordered_map<std::string,std::string>& keys() const { return m_keys; }
    const std::string& command() const { return m_cmd; }
    LineType type() const { return m_type; }

    bool ignorable() const { return m_type==LineType::COMMENT || m_type==LineType::EMPTY; }
    void mark_eof() { m_type = LineType::EOF_MARK; }
};

std::ostream& operator<<(std::ostream& os, const Line& L);

} // namespace fem::reader2
