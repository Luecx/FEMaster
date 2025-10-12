/**
 * @file file.h
 * @brief Lightweight file reader that yields normalized `Line` objects.
 *
 * Features:
 *  - Streams lines from a file, normalizes, and classifies them (`LineType`)
 *  - Skips ignorable lines (`COMMENT`, `EMPTY_LINE`) via `next_line()`
 *  - Supports nested includes using the keyword `*INCLUDE,SRC=path`
 *
 * @see line.h
 * @date 12.10.2025
 */

#pragma once
#include <fstream>
#include <memory>
#include <string>
#include <stdexcept>
#include "line.h"

namespace fem {
namespace dsl {

class File;
using FilePtr = std::unique_ptr<File>;

/**
 * @class File
 * @brief Streaming facade for reading and normalizing input lines.
 *
 * The `File` class provides sequential access to normalized `Line` records.
 * It transparently resolves nested includes when encountering
 * `*INCLUDE,SRC=...` by opening the referenced file as a sub-stream.
 */
class File {
private:
    Line _line;                 ///< Last parsed/returned line (normalized).
    FilePtr _sub_file;          ///< Active sub-file for `*INCLUDE`.
    std::ifstream _stream;      ///< Underlying file stream.

    /**
     * @brief Opens a sub-file for nested include processing.
     *
     * @param name Path to the file specified by `SRC=...`.
     */
    void open_sub_file(const std::string& name) {
        _sub_file = std::make_unique<File>(name);
    }

public:
    /**
     * @brief Constructs a reader bound to a file path.
     *
     * @param file Path to the input deck file to open.
     *
     * @throws std::runtime_error If the file cannot be opened.
     */
    explicit File(const std::string& file) {
        _stream.open(file);
        if (!_stream.is_open())
            throw std::runtime_error("cannot open file: " + file);
    }

    /**
     * @brief Returns the next raw line (including comments/empty) and follows includes.
     *
     * This function:
     *  - Yields pending lines from an active sub-file first.
     *  - Reads one raw line from the current stream and normalizes it via `Line::operator=`.
     *  - If the line is a keyword `INCLUDE`, it attempts to open `SRC=...` as a sub-file and
     *    immediately continues reading from there.
     *
     * @return Reference to the internally stored `Line`.
     */
    Line& next() {
        // Prefer sub-file if active
        if (_sub_file && !_sub_file->is_eof()) {
            Line& l = _sub_file->next();
            if (l.type() != END_OF_FILE)
                return l;
        }

        std::string str;
        if (std::getline(_stream, str)) _line = str;
        else { _line = ""; _line.eof(); }

        // Handle *INCLUDE,SRC=...
        if (_line.type() == KEYWORD_LINE && _line.command() == "INCLUDE") {
            // Crude extraction: parse the normalized keyword buffer for "SRC="
            const std::string s = _line.line();
            auto pos = s.find("SRC=");
            if (pos != std::string::npos)
                open_sub_file(s.substr(pos + 4));
            return next();
        }

        return _line;
    }

    /**
     * @brief Returns the next non-ignorable line, skipping comments and empty lines.
     *
     * @return Reference to the internally stored `Line`.
     */
    Line& next_line() {
        do { _line = next(); } while (_line.ignorable());
        return _line;
    }

    /**
     * @brief Indicates whether the most recently returned line is `END_OF_FILE`.
     *
     * @return `true` if the last returned line was `END_OF_FILE`, otherwise `false`.
     */
    bool is_eof() {
        return _line.type() == END_OF_FILE;
    }
};

} // namespace dsl
} // namespace fem
