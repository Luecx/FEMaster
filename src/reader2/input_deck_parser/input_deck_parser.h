#pragma once
/**
 * @file InputDeckParser.h
 * @brief Two-stage input-deck front-end (analyse + process) built on fem::reader2.
 *
 * Keeps the new library untouched. No helper register_* functions, no adapter callbacks.
 * All parsing registrations live inline in process().
 */

#include <memory>
#include <string>
#include <vector>

#include "../model/model.h"          // fem::model::ModelPtr
#include "writer.h"                  // your existing Writer

// --- new parsing library (do not modify) ---
#include "../command_api.h"   // Registry, CommandRegistrar, Pattern, Condition, …
#include "../reader.h"        // fem::reader2::Reader
#include "../file.h"          // fem::reader2::File
#include "../line.h"          // fem::reader2::Line

namespace fem::reader2 {

/**
 * @brief High-level input-deck parser that mirrors the legacy two-phase API.
 *
 * Phase 1: analyse()  — fast scan to determine maxima and basic info.
 * Phase 2: process()  — register commands/patterns and parse via fem::reader2::Reader.
 */
class InputDeckParser {
    public:
    /**
     * @brief Construct from input and output paths.
     */
    InputDeckParser(const std::string& file_path, const std::string& output_path);

    /**
     * @brief Execute the full pipeline: analyse → build model → process → close writer.
     */
    void read();

    private:
    // -------- legacy-like state we still need --------
    struct ProblemData {
        int highest_node_id     = -1;
        int highest_element_id  = -1;
        int highest_surface_id  = -1;
        int current_loadcase_num = 1;
    };

    std::string          m_file_path;
    std::string          m_output_path;
    Writer               m_writer;
    ProblemData          m_data{};
    fem::model::ModelPtr m_model;

    // -------- phases --------
    void analyse();  ///< tolerant scan using fem::reader2::File/Line
    void process();  ///< inline registration + fem::reader2::Reader
};

} // namespace fem::reader2
