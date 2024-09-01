#include <iostream>
#include <string>
#include <sstream>
#include <type_traits>
#include <iomanip>

namespace logging{


template<typename T>
std::string inline process(T t, std::ostringstream& oss) {
    oss << t;
    return oss.str();
}

template<> inline std::string process<int>(int i, std::ostringstream& oss) {
    oss << std::to_string(i);
    return oss.str();
}

template<> inline std::string process<float>(float f, std::ostringstream& oss) {
    oss << f;
    return oss.str();
}

template<> inline std::string process<double>(double d, std::ostringstream& oss) {
    oss << d;
    return oss.str();
}

template<> inline std::string process<std::ios_base& (*)(std::ios_base&)>(std::ios_base& (*f)(std::ios_base&), std::ostringstream& oss) {
    oss << f;
    return "";
}

template<typename T, typename... Args>
void inline log_impl(const std::string& log_type, T t, Args... args) {
    std::ostringstream oss;
    oss << log_type << logging::get_indentation() << " ";
    process(t, oss);
    ((process(args, oss)), ...);
    std::cout << oss.str() << std::endl;
}

template<typename... Args>
void inline warning(bool condition, Args... args) {
    if (!condition) {
        log_impl("[WARNING]", args...);
    }
}

template<typename... Args>
void inline info(bool condition, Args... args) {
    if (condition) {
        log_impl("[INFO]", args...);
    }
}

template<typename... Args>
void inline error(bool condition, Args... args) {
    if (!condition) {
        std::ostringstream oss;
        ((process(args, oss)), ...);
        std::string message = oss.str();
        std::cout << "[ERROR] " << message << std::endl;
        std::exit(-1);
    }
}

}