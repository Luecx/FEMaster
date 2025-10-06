// parser/reader.h
#pragma once
/**
 * @file reader.h
 * @brief Reader that streams a File, manages scopes, and dispatches to CommandSpec handlers.
 */

#include <optional>
#include <string>
#include <vector>
#include <functional>

#include "file.h"
#include "registry.h"
#include "context.h"
#include "command_spec.h"

namespace fem::reader2 {

    class Reader {
    public:
        /// Run the reader on a given file path.
        void run(const std::string& filepath);

    private:
        struct ActiveBlock {
            std::string scope_name;                   // e.g. "MATERIAL", "LOADCASE", etc.
            std::function<void(Context&)> on_exit_cb; // optional on-exit for the opening command
        };

        Context cx_;
        std::optional<Line> look_;
        Line scratch_; // to return a stable ref when consuming look_
        std::vector<ActiveBlock> blocks_; // stack of open scopes (top = current)

        Line& next_effective(File& f);
        Line& peek_effective(File& f);

        // Bubble up until we find a scope that allows `name`; call on_exit while popping.
        // Returns the scope index (0..top) where the command belongs. Throws if not found at ROOT.
        size_t ensure_scope_accepts_or_bubble(const std::string& name);
    };

} // namespace fem::reader2
