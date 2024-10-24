#pragma once

#include <iostream>
#include <stdexcept>
#include <string>

namespace fem::logging{
inline int indentation_level = 0;

inline void up() {
    indentation_level++;
}

inline void down() {
    if (indentation_level > 0) {
        indentation_level--;
    }
}

inline std::string get_indentation() {
    std::string res{};
    for(int i = 0; i < indentation_level; i++){
        res += "  ";
    }
    return res;
}

// Forward declarations of the functions
template<typename... Args> inline void warning(bool condition, Args... args);
template<typename... Args> inline void info(bool condition, Args... args);
template<typename... Args> inline void error(bool condition, Args... args);

} // namespace fem::logging

#include "logging.ipp"

