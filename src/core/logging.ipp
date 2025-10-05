/******************************************************************************
 * @file logging.ipp
 * @brief Implements templated logging helpers.
 *
 * @see src/core/logging.h
 ******************************************************************************/

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace fem {
namespace logging {

int indentation_level = 0;

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
        result += "  ";
    }
    return result;
}

namespace detail {

template<typename T>
std::string process(T value, std::ostringstream& stream) {
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

template<typename T, typename... Args>
void log_impl(const std::string& prefix, T value, Args... args) {
    std::ostringstream stream;
    stream << prefix << logging::get_indentation() << " ";
    process(value, stream);
    ((process(args, stream)), ...);
    std::cout << stream.str() << std::endl;
}

} // namespace detail

template<typename... Args>
void warning(bool condition, Args... args) {
    if (!condition) {
        detail::log_impl("[WARNING]", args...);
    }
}

template<typename... Args>
void info(bool condition, Args... args) {
    if (condition) {
        detail::log_impl("[INFO]", args...);
    }
}

template<typename... Args>
void error(bool condition, Args... args) {
    if (!condition) {
        std::ostringstream stream;
        ((detail::process(args, stream)), ...);
        std::cout << "[ERROR] " << stream.str() << std::endl;
        std::exit(-1);
    }
}

} // namespace logging
} // namespace fem
