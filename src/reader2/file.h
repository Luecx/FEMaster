#pragma once
#include <fstream>
#include <memory>
#include <string>
#include "line.h"

namespace fem::reader2 {

class File {
public:
    explicit File(const std::string& path);
    ~File();

    // Read next *raw* line and parse it into Line (including INCLUDE handling at keyword level)
    Line& next();

    // Read next non-ignorable line (skips EMPTY/COMMENT); sets EOF_MARK at end
    Line& next_effective();

    bool eof() const { return m_eof; }

    const std::string& path() const { return m_path; }
    int line_number() const { return m_line_no; }

private:
    Line m_cur;
    std::ifstream m_stream;
    std::unique_ptr<File> m_included; // nested include
    std::string m_path;
    int m_line_no = 0;
    bool m_eof = false;

    void open_include(const std::string& path);
};

} // namespace fem::reader2
