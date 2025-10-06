#pragma once
#include <optional>
#include <string>
#include <vector>
#include <functional>

#include "file.h"
#include "registry.h"
#include "context.h"
#include "command_spec.h"
#include "pattern.h"

namespace fem::reader2 {

class Reader {
public:
    /**
     * \brief Execute the reader pipeline for the provided file.
     */
    void run(const std::string& filepath);

private:
    /// \brief Tracks an open scope with its exit callback.
    struct ActiveBlock {
        std::string                 scope_name;
        std::function<void(Context&)> on_exit_cb;
    };

    Context                 _context;
    std::optional<Line>     _look;            ///< Current buffered line.
    size_t                  _tokenIndex = 0; ///< Index within the buffered data line.
    std::vector<ActiveBlock> _blocks;        ///< Stack of open blocks.

    /// \brief Peek the next non-ignorable line without consuming it fully.
    Line& peek_effective(File& file);

    // token helpers for DATA lines
    /// \brief Ensure the buffered line exists and represents data.
    bool ensure_data_line(File& file);
    /// \brief Count remaining tokens on the buffered data line.
    size_t tokens_available() const;
    /// \brief Transfer up to count tokens into the output vector.
    size_t take_tokens(size_t count, std::vector<std::string>& out, const std::string& file, int& last_line_no);
    /// \brief Release the buffered line once consumed.
    void maybe_pop_line();

    /// \brief Bubble scopes until one accepts the specified command.
    size_t ensure_scope_accepts_or_bubble(const std::string& name);
};

} // namespace fem::reader2
