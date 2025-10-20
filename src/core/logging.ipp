/**
 * @file logging.ipp
 * @brief Implements templated logging helpers with automatic line wrapping.
 *
 * @see src/core/logging.h
 */

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fem {
namespace logging {

inline int indentation_level = 0;
inline std::size_t g_console_width = 102; // default; configurable

inline void set_console_width(std::size_t cols) {
    g_console_width = std::max<std::size_t>(cols, 20);
}

inline std::size_t get_console_width() {
    return g_console_width;
}

inline void up() {
    ++indentation_level;
}

inline void down() {
    if (indentation_level > 0) {
        --indentation_level;
    }
}

inline std::string get_indentation() {
    std::string result;
    for (int i = 0; i < indentation_level; ++i) {
        result += "  "; // 2 spaces per level
    }
    return result;
}

namespace detail {

// --- Small helpers ----------------------------------------------------------

inline bool is_space_or_tab(char c) {
    return c == ' ' || c == '\t';
}

inline std::string leading_ws_of(const std::string& s) {
    std::size_t i = 0;
    while (i < s.size() && is_space_or_tab(s[i])) ++i;
    return s.substr(0, i);
}

// Split by '\n' but keep empty segments (so trailing/leading newlines are honored)
inline std::vector<std::string> split_lines_hard(const std::string& s) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\n') {
            parts.emplace_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
    parts.emplace_back(s.substr(start));
    return parts;
}

// Emit a single physical line to std::cout
inline void emit_line(const std::string& prefix, const std::string& text) {
    std::cout << prefix << text << std::endl;
}

/**
 * @brief Wrap a single logical segment (no '\n' inside) across multiple
 *        physical lines. The first line uses `base_prefix`, continuation lines
 *        use `base_prefix + seg_leading_ws`.
 *
 * @param base_prefix    The fixed prefix: e.g. "[INFO]" + indentation + " "
 * @param segment        The text segment (no '\n')
 * @param max_cols       Console width in characters
 */
inline void wrap_and_emit_segment(const std::string& base_prefix,
                                  const std::string& segment,
                                  std::size_t max_cols)
{
    // leading WS of this segment (spaces/tabs only) — replicated on wrapped lines
    const std::string seg_leading_ws = leading_ws_of(segment);

    // For the *first* physical line of this segment:
    std::string cur_prefix = base_prefix;
    std::size_t cur_prefix_len = cur_prefix.size();

    // For *continuation* lines:
    const std::string cont_prefix = base_prefix + seg_leading_ws;
    const std::size_t cont_prefix_len = cont_prefix.size();

    // Work on a view into segment
    std::size_t pos = 0;
    const std::size_t N = segment.size();

    // Helper to compute available width for current prefix
    auto avail_for_text = [&](std::size_t prefix_len) -> std::size_t {
        if (max_cols <= prefix_len) return 0;
        return max_cols - prefix_len;
    };

    bool first_line = true;

    while (pos < N) {
        // Decide prefix for this physical line
        const std::string& prefix = first_line ? cur_prefix : cont_prefix;
        const std::size_t prefix_len = first_line ? cur_prefix_len : cont_prefix_len;

        std::size_t avail = avail_for_text(prefix_len);
        if (avail == 0) {
            // Degenerate case: no room; emit empty line to avoid infinite loop.
            emit_line(prefix, "");
            first_line = false;
            continue;
        }

        // Remaining text of this segment for this line
        const char* data = segment.data() + pos;
        std::size_t remaining = N - pos;

        if (remaining <= avail) {
            // Fits entirely
            emit_line(prefix, std::string(data, remaining));
            break;
        }

        // Try to break at last whitespace within the window
        std::size_t break_pos = avail;
        // Scan backward to find a natural break on space/tab
        std::size_t i = avail;
        while (i > 0 && !is_space_or_tab(data[i - 1])) --i;
        if (i > 0) {
            break_pos = i;
        } else {
            // No whitespace: hard-break at avail
            break_pos = avail;
        }

        // Emit the slice
        emit_line(prefix, std::string(data, break_pos));

        // Advance position; if we broke on space(s), skip consecutive spaces/tabs
        pos += break_pos;
        while (pos < N && is_space_or_tab(segment[pos])) ++pos;

        first_line = false;
    }

    // Special case: empty segment -> still emit a blank line with prefix
    if (N == 0) {
        emit_line(base_prefix, "");
    }
}

/**
 * @brief Compose and emit a wrapped multi-line message.
 *
 * We first build:
 *   base_prefix = prefix + get_indentation() + " "
 * Then we stream the content (args...) into a string, split at hard '\n', and
 * wrap each segment separately. For each segment, its own leading whitespace
 * (spaces/tabs) is detected and reused for continuation lines.
 *
 * Empty content leads to a single empty line with prefix.
 */
template<typename... Args>
void log_wrapped_with_prefix(const std::string& prefix, Args&&... args) {
    // Build the *content* in a separate stream to keep formatting manipulators applied only there.
    std::ostringstream content_stream;
    (content_stream << ... << std::forward<Args>(args));
    const std::string content = content_stream.str();

    const std::string base_prefix = prefix + logging::get_indentation() + " ";
    const std::size_t max_cols = get_console_width();

    if (content.empty()) {
        // Always emit at least one blank line with prefix
        std::cout << base_prefix << std::endl;
        return;
    }

    // Split on hard newlines and wrap each segment
    const auto segments = split_lines_hard(content);
    for (std::size_t s = 0; s < segments.size(); ++s) {
        wrap_and_emit_segment(base_prefix, segments[s], max_cols);
    }
}

// --- original templated streaming helpers (kept for error() path) -----------
template<typename T>
inline std::string process(T value, std::ostringstream& stream) {
    stream << value;
    return stream.str();
}

template<>
inline std::string process<int>(int value, std::ostringstream& stream) {
    stream << value;
    return stream.str();
}

template<>
inline std::string process<float>(float value, std::ostringstream& stream) {
    stream << value;
    return stream.str();
}

template<>
inline std::string process<double>(double value, std::ostringstream& stream) {
    stream << value;
    return stream.str();
}

template<>
inline std::string process<std::ios_base& (*)(std::ios_base&)>(std::ios_base& (*manip)(std::ios_base&), std::ostringstream& stream) {
    stream << manip;
    return stream.str();
}

} // namespace detail

// --- Public APIs ------------------------------------------------------------

template<typename... Args>
void warning(bool condition, Args... args) {
    if (!condition) {
        detail::log_wrapped_with_prefix("[WARNING]", std::forward<Args>(args)...);
    }
}

template<typename... Args>
void info(bool condition, Args... args) {
    if (condition) {
        detail::log_wrapped_with_prefix("[INFO]", std::forward<Args>(args)...);
    }
}

template<typename... Args>
void error(bool condition, Args... args) {
    if (!condition) {
        // Build content first to include in [ERROR] header
        std::ostringstream s;
        (detail::process(args, s), ...);
        detail::log_wrapped_with_prefix("[ERROR]", s.str());
        std::exit(-1);
    }
}

} // namespace logging
} // namespace fem
