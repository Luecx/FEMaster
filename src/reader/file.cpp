//
// Created by Luecx on 04.09.2023.
//
#include "file.h"

namespace fem::reader {

void fem::reader::File::open_sub_file(const std::string& name) {
    sub_file = std::make_unique<File>(name);
}
fem::reader::File::File(const std::string& file) {
    stream.open(file);
    logging::error(stream.is_open(), "cannot open file: ", file);
}
Line& fem::reader::File::next() {
    // try reading from the sub file first
    if (sub_file && !sub_file->is_eof()) {
        Line& l = sub_file->next();
        if (l.type() != END_OF_FILE) {
            return l;
        }
    }
    std::string str;
    if (std::getline(stream, str)) {
        m_line = str;
    } else {
        m_line = "";
        m_line.eof();
    }
    // already parse opening a new file
    if (m_line.type() == KEYWORD_LINE && m_line.command() == "INCLUDE") {
        open_sub_file(m_line.require<std::string>("SRC"));
        return next();
    }

    return m_line;
}

Line& fem::reader::File::next_line() {
    do {
        m_line = next();
    } while (m_line.ignorable());
    return m_line;
}

bool fem::reader::File::is_eof() {
    return m_line.type() == END_OF_FILE;
}

}    // namespace fem::reader