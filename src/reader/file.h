#pragma once

#include "../core/core.h"
#include "line.h"

#include <fstream>

namespace fem {
namespace reader {

class File;
using FilePtr = std::unique_ptr<File>;

class File {

    Line          m_line;
    FilePtr       sub_file;
    std::ifstream stream;

    private:
    void open_sub_file(const std::string& name);

    public:
    File(const std::string& file);

    Line& next();
    Line& next_line();

    bool is_eof();
};
}    // namespace reader
}    // namespace fem
