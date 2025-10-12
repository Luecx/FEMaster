
#pragma once
#include <fstream>
#include <memory>
#include <string>
#include "line.h"

namespace fem { namespace reader {

class File;
using FilePtr = std::unique_ptr<File>;

class File {
    Line m_line;
    FilePtr sub_file;
    std::ifstream stream;

    void open_sub_file(const std::string& name) {
        sub_file = std::make_unique<File>(name);
    }

public:
    explicit File(const std::string& file) {
        stream.open(file);
        if (!stream.is_open()) throw std::runtime_error("cannot open file: " + file);
    }

    Line& next() {
        if (sub_file && !sub_file->is_eof()) {
            Line& l = sub_file->next();
            if (l.type() != END_OF_FILE) return l;
        }
        std::string str;
        if (std::getline(stream, str)) m_line = str;
        else { m_line = ""; m_line.eof(); }
        if (m_line.type()==KEYWORD_LINE && m_line.command()=="INCLUDE") {
            // crude: require key SRC
            // We re-parse again from m_line.line() to find SRC=...
            // But for demo keep it simple: assume "*INCLUDE,SRC=path"
            std::string s = m_line.line();
            auto pos = s.find("SRC=");
            if (pos != std::string::npos) open_sub_file(s.substr(pos+4));
            return next();
        }
        return m_line;
    }

    Line& next_line() {
        do { m_line = next(); } while (m_line.ignorable());
        return m_line;
    }

    bool is_eof() { return m_line.type()==END_OF_FILE; }
};

}} // namespace fem::reader
