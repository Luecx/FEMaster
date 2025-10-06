#pragma once
#include <string>
#include <vector>
#include <limits>

namespace fem::reader2 {

inline constexpr size_t kInf = std::numeric_limits<size_t>::max();

/// Basic line kinds detected by the parser
enum class LineType { COMMENT=-1, EMPTY=0, KEYWORD=1, DATA=2, EOF_MARK=3 };

/// Metadata attached to a parsed line
struct LineMeta {
    std::string file;
    int         line_no = 0;
};

/// Alias for symbolic scope identifiers (e.g. "ROOT", "MATERIAL", "LOADCASE").
using Scope = std::string;

} // namespace fem::reader2
