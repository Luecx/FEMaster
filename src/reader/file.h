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
    void open_sub_file(const std::string& name){
        sub_file = std::make_unique<File>(name);
    }

    public:
    File(const std::string& file) {
        stream.open(file);
        log_error(stream.is_open(), "cannot open file: ", file);
    }

    Line& next() {
        // try reading from the sub file first
        if(sub_file && !sub_file->is_eof()){
            Line& l = sub_file->next();
            if(l.type() != END_OF_FILE){
                return l;
            }
        }
        std::string str;
        if (std::getline(stream, str)){
            m_line = str;
        }else{
            m_line = "";
            m_line.eof();
        }
        // already parse opening a new file
        if(m_line.type() == KEYWORD_LINE && m_line.command() == "INCLUDE"){
            open_sub_file(m_line.require<std::string>("SRC"));
            return next();
        }

        return m_line;
    }

    bool is_eof(){
        return m_line.type() == END_OF_FILE;
    }
};
}    // namespace reader
}    // namespace fem
