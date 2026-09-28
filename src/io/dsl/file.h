/**
 * @file file.h
 * @brief Lightweight file reader that yields normalized `Line` objects.
 *
 * Features:
 *  - Streams lines from a file, normalizes, and classifies them (`LineType`)
 *  - Skips ignorable lines (`COMMENT`, `EMPTY_LINE`) via `next_line()`
 *  - Supports nested includes using `*INCLUDE,SRC=path` and Abaqus `*INCLUDE,INPUT=path`
 *
 * @see line.h
 * @date 12.10.2025
 */

#pragma once
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <stdexcept>
#include "line.h"

namespace fem {
namespace io {
namespace dsl {

class File;
using FilePtr = std::unique_ptr<File>;

/**
 * @class File
 * @brief Streaming facade for reading and normalizing input lines.
 *
 * The `File` class provides sequential access to normalized `Line` records.
 * It transparently resolves nested includes when encountering
 * `*INCLUDE,SRC=...` or Abaqus-style `*INCLUDE,INPUT=...` by opening the
 * referenced file as a sub-stream.
 */
class File {
private:
    Line                               _line;          ///< Last parsed/returned line (normalized).
    FilePtr                            _sub_file;      ///< Active sub-file for `*INCLUDE`.

    // Path of the current input file
    std::filesystem::path _path;

    std::ifstream                      _stream;        ///< Underlying file stream.
    std::shared_ptr<const std::string> _source;        ///< Shared source path for returned lines.
    std::size_t                        _line_number{}; ///< Current one-based source line number.

    /**
     * @brief Opens a sub-file for nested include processing.
     *
     * Relative paths are resolved against the directory of the including file.
     * If no such file exists, the previous working-directory-relative behavior
     * is retained as a fallback for backwards compatibility.
     */
    void open_sub_file(const std::string& name) {
        std::filesystem::path path(name);
        if (path.is_relative()) {
            const std::filesystem::path local = _path.parent_path() / path;
            if (std::filesystem::exists(local))
                path = local;
        }
        _sub_file = std::make_unique<File>(path.lexically_normal().string());
    }

public:
    /**
     * @brief Constructs a reader bound to a file path.
     *
     * @param file Path to the input deck file to open.
     *
     * @throws std::runtime_error If the file cannot be opened.
     */
    explicit File(const std::string& file)
        : _path(file),
          _stream(file),
          _source(std::make_shared<const std::string>(file)) {
        if (!_stream.is_open())
            throw std::runtime_error("cannot open file: " + file);
    }

    /**
     * @brief Returns the next raw line (including comments/empty) and follows includes.
     *
     * This function:
     *  - Yields pending lines from an active sub-file first.
     *  - Reads one raw line from the current stream and normalizes it via `Line::operator=`.
     *  - If the line is a keyword `INCLUDE`, it opens `INPUT=...` or `SRC=...` as a sub-file
     *    and immediately continues reading from there.
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
        if (std::getline(_stream, str)) {
            _line = str;
            _line.set_location(_source, ++_line_number);
        } else {
            _line = "";
            _line.eof();
            _line.set_location(_source, _line_number + 1);
        }

        // Parse INCLUDE from the original source line so the file-system path keeps its case
        if (_line.type() == KEYWORD_LINE && _line.command() == "INCLUDE") {
            auto trim = [](std::string value) {
                const auto not_space = [](unsigned char ch) { return !std::isspace(ch); };
                value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
                value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
                return value;
            };

            std::string path;
            for (std::size_t pos = str.find(','); pos != std::string::npos;) {
                const std::size_t next = str.find(',', pos + 1);
                std::string token = trim(str.substr(pos + 1, next == std::string::npos ? next : next - pos - 1));

                const std::size_t eq = token.find('=');
                if (eq != std::string::npos) {
                    std::string key = trim(token.substr(0, eq));
                    std::transform(key.begin(), key.end(), key.begin(),
                                   [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });

                    if (key == "INPUT" || key == "SRC") {
                        path = trim(token.substr(eq + 1));
                        if (path.size() >= 2 &&
                            ((path.front() == '"' && path.back() == '"') ||
                             (path.front() == '\'' && path.back() == '\'')))
                            path = path.substr(1, path.size() - 2);
                        break;
                    }
                }

                pos = next;
            }

            if (path.empty())
                throw std::runtime_error("INCLUDE requires INPUT=... or SRC=... at " + _line.location().str());

            open_sub_file(path);
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
} // namespace io
} // namespace fem
