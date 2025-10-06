#pragma once

#include <fstream>
#include <memory>
#include <string>

#include "line.h"

namespace fem::reader2 {

/**
 * \brief Lightweight wrapper around an input file that expands INCLUDE statements.
 */
class File {
public:
    /**
     * \brief Construct and immediately open the provided filesystem path.
     * \param path Absolute or relative path to the input file.
     */
    explicit File(const std::string& path);

    /**
     * \brief Default virtual destructor for safe inheritance.
     */
    ~File();

    /**
     * \brief Read the next raw line, honouring nested INCLUDEs.
     * \return Reference to the reusable Line buffer that now holds the parsed information.
     */
    Line& next();

    /**
     * \brief Fetch the next non-ignorable line, skipping blanks and comments.
     * \return Reference to the reusable Line buffer populated with the next effective line.
     */
    Line& next_effective();

    /**
     * \brief Query whether the logical end of the file (including includes) has been reached.
     */
    bool eof() const
    {
        return _eof;
    }

    /**
     * \brief Access the filesystem path that backs this File instance.
     */
    const std::string& path() const
    {
        return _path;
    }

    /**
     * \brief Retrieve the current line number within the active physical file.
     */
    int line_number() const
    {
        return _lineNumber;
    }

private:
    Line                     _current;      ///< Reusable line buffer.
    std::ifstream            _stream;       ///< Stream for the current physical file.
    std::unique_ptr<File>    _included;     ///< Active include file, if any.
    std::string              _path;         ///< Path to the current file.
    int                      _lineNumber = 0; ///< Current 1-based line number.
    bool                     _eof = false;  ///< True once EOF has been produced.

    /**
     * \brief Open a nested include file specified by the SRC key.
     */
    void open_include(const std::string& path);
};

} // namespace fem::reader2
