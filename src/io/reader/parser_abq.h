/**
 * @file parser_abq.h
 * @brief Declares the Abaqus input-deck parser.
 *
 * The Abaqus reader reuses the native four-stage parser pipeline and replaces
 * only the command registration for each stage.
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include "parser.h"

namespace fem::io::reader {

/**
 * @class ParserAbq
 * @brief Abaqus syntax specialization of the staged input parser.
 */
class ParserAbq : public Parser {
public:
    ParserAbq() = default;
    ~ParserAbq() override = default;

protected:
    void configure_count_stage(io::dsl::Registry& registry, CountData& count) override;
    void configure_topology_stage(io::dsl::Registry& registry) override;
    void configure_field_stage(io::dsl::Registry& registry) override;
    void configure_data_stage(io::dsl::Registry& registry) override;
};

} // namespace fem::io::reader
