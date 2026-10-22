/**
 * @file memory.h
 * @brief Cross-platform memory query utilities.
 *
 * This header provides a function to estimate the available physical RAM
 * on the system in bytes through fem::memory. It supports Windows, macOS
 * and Linux platforms; callers use the estimate to budget temporary storage.
 * Linux also observes visible cgroup v1/v2 memory limits and ancestor usage.
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

namespace fem::memory {

/**
 * @brief Returns an estimate of available physical RAM in bytes.
 *
 * Linux caps host MemAvailable by the remaining hard-limit budget of the
 * process memory cgroup and its visible ancestors. Limits and usage are queried
 * together on each call; the result is a snapshot, not a memory reservation.
 *
 * @return Available memory in bytes, or 0 if the query fails or a limit is exhausted.
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

    // Query host memory before applying process-specific cgroup limits
    std::ifstream file("/proc/meminfo");
    std::string   line;
    std::uint64_t available = 0;

    while (std::getline(file, line)) {
        if (line.rfind("MemAvailable:", 0) == 0) {
            std::istringstream stream(line.substr(13));
            std::uint64_t      kib = 0;
            if (stream >> kib) {
                available = kib * 1024;
            }
            break;
        }
    }
    if (available == 0) {
        return 0;
    }

    // Mountinfo encodes whitespace and backslashes in paths as octal escapes
    auto decode_path = [](const std::string& encoded) {
        std::string path;
        for (std::size_t i = 0; i < encoded.size(); ++i) {
            if (encoded[i] == '\\' && i + 3 < encoded.size() &&
                encoded[i + 1] >= '0' && encoded[i + 1] <= '7' &&
                encoded[i + 2] >= '0' && encoded[i + 2] <= '7' &&
                encoded[i + 3] >= '0' && encoded[i + 3] <= '7') {
                path += static_cast<char>((encoded[i + 1] - '0') * 64 +
                                          (encoded[i + 2] - '0') * 8  +
                                          (encoded[i + 3] - '0'));
                i += 3;
            } else {
                path += encoded[i];
            }
        }
        return path;
    };

    // Locate the actual controller mounts, including combined v1 controllers
    std::ifstream groups("/proc/self/cgroup");
    std::ifstream mounts("/proc/self/mountinfo");
    if (!groups || !mounts) {
        return 0;
    }
    while (std::getline(groups, line)) {
        const auto first  = line.find(':');
        const auto second = first == std::string::npos ? first : line.find(':', first + 1);
        if (second == std::string::npos) {
            continue;
        }
        const std::string controllers = "," + line.substr(first + 1, second - first - 1) + ",";
        const std::string group       = line.substr(second + 1);
        const bool        unified     = controllers == ",,";
        if (!unified && controllers.find(",memory,") == std::string::npos) {
            continue;
        }

        // A membership outside the namespace root cannot be resolved safely
        if ((group + "/").find("/../") != std::string::npos) {
            return 0;
        }

        mounts.clear();
        mounts.seekg(0);
        std::string mount_line;
        while (std::getline(mounts, mount_line)) {
            const auto separator = mount_line.find(" - ");
            if (separator == std::string::npos) {
                continue;
            }
            std::istringstream metadata(mount_line.substr(separator + 3));
            std::string        type, source, options;
            metadata >> type >> source >> options;
            if (unified ? type != "cgroup2" :
                (type != "cgroup" || ("," + options + ",").find(",memory,") == std::string::npos)) {
                continue;
            }

            // Membership and mount roots share the process cgroup namespace.
            // Only mounts containing this membership can impose its budget.
            std::istringstream fields(mount_line.substr(0, separator));
            std::string        mount_id, parent_id, device, root, mount;
            if (!(fields >> mount_id >> parent_id >> device >> root >> mount)) {
                continue;
            }
            root  = decode_path(root);
            mount = decode_path(mount);
            // Inherited mounts can refer outside the current cgroup namespace
            if ((root + "/").find("/../") != std::string::npos) {
                return 0;
            }
            if (root.empty() || mount.empty() || group.empty() || group[0] != '/') {
                continue;
            }
            if (root != "/" && group != root && group.rfind(root + "/", 0) != 0) {
                continue;
            }
            const std::string relative = root == "/" ? group : group.substr(root.size());
            const std::string leaf     = relative.empty() || relative == "/" ? mount :
                                         (mount == "/" ? relative : mount + relative);
            std::string       path     = leaf;

            // A child may be unlimited while an ancestor enforces a hard limit.
            // Stop at the mounted root; ancestors outside this mount are hidden.
            while (true) {
                bool applies = true;
                if (!unified && path != leaf) {
                    std::ifstream hierarchy_file(path + "/memory.use_hierarchy");
                    int hierarchical = 0;
                    applies = (hierarchy_file >> hierarchical) && hierarchical != 0;
                }
                std::ifstream limit_file(path + (unified ? "/memory.max" : "/memory.limit_in_bytes"));
                std::string   limit_text;
                if (applies && (limit_file >> limit_text) && limit_text != "max") {
                    std::istringstream limit_stream(limit_text);
                    std::uint64_t      limit = 0;
                    std::uint64_t      usage = 0;
                    std::ifstream usage_file(path + (unified ? "/memory.current" : "/memory.usage_in_bytes"));
                    if (limit_text[0] == '-' || !(limit_stream >> limit) || !(usage_file >> usage)) {
                        return 0;
                    }
                    const std::uint64_t remaining = usage < limit ? limit - usage : 0;
                    if (remaining < available) {
                        available = remaining;
                    }
                }
                if (path == mount || available == 0) {
                    break;
                }
                const auto parent = path.rfind('/');
                path.resize(parent == 0 ? 1 : parent);
            }
        }
    }
    return available;

#else
    logging::error(false, "Unsupported platform, cannot determine available memory");
    return 0;
#endif
}

} // namespace fem::memory
