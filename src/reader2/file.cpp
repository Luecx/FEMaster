#include "file.h"
#include <stdexcept>

namespace fem::reader2 {

File::File(const std::string& path): m_path(path) {
    m_stream.open(path);
    if (!m_stream.is_open()) throw std::runtime_error("Cannot open file: "+path);
}

File::~File() = default;

void File::open_include(const std::string& path){
    m_included = std::make_unique<File>(path);
}

Line& File::next(){
    // try include first
    if (m_included){
        Line& inc = m_included->next();
        if (inc.type()!=LineType::EOF_MARK) return inc;
        m_included.reset();
    }

    std::string s;
    if (std::getline(m_stream, s)){
        ++m_line_no;
        m_cur.assign(s);
        // auto-include if keyword *INCLUDE,SRC=...
        if (m_cur.type()==LineType::KEYWORD && m_cur.command()=="INCLUDE"){
            auto it = m_cur.keys().find("SRC");
            if (it!=m_cur.keys().end()) {
                open_include(it->second);
                return next();
            }
        }
        return m_cur;
    }
    m_eof = true;
    m_cur.assign("");
    m_cur.mark_eof();
    return m_cur;
}

Line& File::next_effective(){
    do { next(); } while (!eof() && m_cur.ignorable());
    return m_cur;
}

} // namespace fem::reader2
