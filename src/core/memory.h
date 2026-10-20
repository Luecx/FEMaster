/**
 * @file memory.h
 * @brief Cross-platform memory query utilities.
 *
 * This header provides a function to estimate the available physical RAM
 * on the system in bytes. It supports Windows, macOS, and Linux platforms.
 */
#pragma once

#include <cstdint>
#include "logging.h"

#if defined(_WIN32)
    #include <windows.h>
#elif defined(__APPLE__)
    #include <mach/mach.h>
#elif defined(__linux__)
    #include <fstream>
    #include <sstream>
    #include <string>
#endif


/**
 * @brief Returns an estimate of available physical RAM in bytes.
 * @return Available memory in bytes, or 0 if the query fails.
 */
inline std::uint64_t available_memory_bytes() {
#if defined(_WIN32)
    MEMORYSTATUSEX status{};
    status.dwLength = sizeof(status);

    if (!GlobalMemoryStatusEx(&status)) {
        return 0;
    }

    return status.ullAvailPhys;

#elif defined(__APPLE__)

    const mach_port_t host = mach_host_self();

    vm_statistics64_data_t stats{};
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;

    if (host_statistics64(
            host,
            HOST_VM_INFO64,
            reinterpret_cast<host_info64_t>(&stats),
            &count) != KERN_SUCCESS) {
        return 0;
            }

    vm_size_t page_size = 0;
    if (host_page_size(host, &page_size) != KERN_SUCCESS) {
        return 0;
    }

    // Free + reclaimable inactive pages (approximation).
    return (std::uint64_t(stats.free_count) +
            std::uint64_t(stats.inactive_count)) * page_size;

#elif defined(__linux__)

    std::ifstream file("/proc/meminfo");
    std::string line;

    while (std::getline(file, line)) {
        if (line.rfind("MemAvailable:", 0) == 0) {
            std::istringstream stream(line.substr(13));
            std::uint64_t kib = 0;

            if (stream >> kib) {
                return kib * 1024;
            }
        }
    }

    return 0;

#else
    logging::error(false, "Unsupported platform, cannot determine available memory");
    return 0;
#endif
}