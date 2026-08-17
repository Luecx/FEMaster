/**
 * @file parser_abq.cpp
 * @brief Implements the empty Abaqus input-deck parser.
 *
 * The reader currently registers no commands. This establishes the format
 * selection and parser entry point without introducing any Abaqus keyword
 * semantics yet.
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#include "parser_abq.h"

#include "../dsl/engine.h"
#include "../dsl/file.h"
#include "../dsl/registry.h"

namespace fem::io::reader {

void ParserAbq::run(const std::string&                   input_path,
                    const std::string&                   output_path,
                    const io::writer::WriterFileFormats& writer_formats) {
    (void) output_path;
    (void) writer_formats;

    io::dsl::Registry registry;
    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
}

} // namespace fem::io::reader
