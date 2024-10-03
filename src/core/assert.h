#pragma once

#ifdef NDEBUG
#define runtime_assert(condition, message)                                                                             \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cout << "runtime assertion triggered on: " << #condition << "\n";                                     \
            std::cout << "message: " << message << "\n";                                                               \
            std::cout << "line   : " << __LINE__ << "\n";                                                              \
            std::cout << "file   : " << __FILE__ << "\n";                                                              \
            throw std::runtime_error(message);                                                                         \
        }                                                                                                              \
    } while (false)
#else
#define runtime_assert(condition, message)
#endif

#define runtime_check(condition, message)                                                                              \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cout << "runtime check triggered on: " << #condition << "\n";                                         \
            std::cout << "message: " << message << "\n";                                                               \
            std::cout << "line   : " << __LINE__ << "\n";                                                              \
            std::cout << "file   : " << __FILE__ << "\n";                                                              \
            throw std::runtime_error(message);                                                                         \
        }                                                                                                              \
    } while (false)
