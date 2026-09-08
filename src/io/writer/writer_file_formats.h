#pragma once

#include <string>

namespace fem {
namespace io {
namespace writer {

struct WriterFileFormats {
    bool res = true;
    bool frd = true;
    bool femr = false;
    std::string result_compression = "none";
};

} // namespace writer
} // namespace io
} // namespace fem
