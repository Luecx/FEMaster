/**
 * @file parser_abq.h
 * @brief Declares the Abaqus input-deck parser.
 *
 * The Abaqus reader is intentionally empty at this stage. It provides the same
 * high-level run entry point as the native FEMaster parser so the executable can
 * select the input syntax through --format. Abaqus keywords are added
 * incrementally in dedicated follow-up changes.
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include "../writer/writers.h"

#include <string>

namespace fem::io::reader {

/**
 * @class ParserAbq
 * @brief Entry point for Abaqus input-deck parsing.
 *
 * No Abaqus keywords are registered yet. Therefore, an empty input deck is
 * accepted while any keyword is reported as unknown by the DSL engine.
 */
class ParserAbq {
public:
    ParserAbq() = default;
    ~ParserAbq() = default;

    /**
     * Parse one Abaqus input deck.
     *
     * @param input_path Input-deck path.
     * @param output_path Output base path reserved for later result handling.
     * @param writer_formats Requested result-writer formats.
     */
    void run(const std::string&                   input_path,
             const std::string&                   output_path,
             const io::writer::WriterFileFormats& writer_formats = io::writer::WriterFileFormats());
};

} // namespace fem::io::reader
