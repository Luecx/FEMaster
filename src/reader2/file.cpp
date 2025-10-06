#include "file.h"
#include <stdexcept>

namespace fem::reader2 {

File::File(const std::string& path) : _path(path)
{
    _stream.open(path);
    if (!_stream.is_open()) {
        throw std::runtime_error("Cannot open file: " + path);
    }
}

File::~File() = default;

void File::open_include(const std::string& path)
{
    _included = std::make_unique<File>(path);
}

Line& File::next()
{
    // try include first
    if (_included) {
        Line& inc = _included->next();
        if (inc.type() != LineType::EOF_MARK) {
            return inc;
        }
        _included.reset();
    }

    std::string s;
    if (std::getline(_stream, s)) {
        ++_lineNumber;
        _current.assign(s);
        // auto-include if keyword *INCLUDE,SRC=...
        if (_current.type() == LineType::KEYWORD && _current.command() == "INCLUDE") {
            auto it = _current.keys().find("SRC");
            if (it != _current.keys().end()) {
                open_include(it->second);
                return next();
            }
        }
        return _current;
    }
    _eof = true;
    _current.assign("");
    _current.mark_eof();
    return _current;
}

Line& File::next_effective()
{
    do {
        next();
    } while (!eof() && _current.ignorable());
    return _current;
}

} // namespace fem::reader2
