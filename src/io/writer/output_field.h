/**
 * @file output_field.h
 * @brief Declares semantic result quantities understood by FEMaster output requests.
 *
 * Output fields describe what a user asks FEMaster to report. They deliberately
 * do not describe where the data comes from: an analysis step may provide a
 * field directly, while the output request handler may recover another field
 * recursively from already available dependencies.
 *
 * Input-deck aliases such as Abaqus/CalculiX U, RF, S and E are normalized to
 * this enum before analysis output is written. Result-file names remain the
 * existing FEMaster names.
 *
 * @see output_request_handler.h
 */

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace fem::io::writer {

/**
 * @brief Identifies one reportable FEMaster result quantity.
 *
 * The first group contains primary solution/state fields commonly supplied by
 * analysis steps. The remaining entries are either analysis-specific direct
 * quantities or values that can be recovered from other fields.
 */
enum class OutputField : std::uint8_t {
    DISPLACEMENT,
    VELOCITY,
    ACCELERATION,
    TEMPERATURE,

    EXTERNAL_FORCES,
    INTERNAL_FORCES,
    REACTION_FORCES,

    STRESS,
    STRAIN,
    PEEQ,
    STRESS_TOP,
    STRESS_BOT,
    SHELL_RESULTANTS,
    LOCAL_SECTION_FORCES,
    SHEAR_FLOW,
    COMPLIANCE,
    VOLUME,
    HEAT_FLUX,

    MODE_SHAPE,
    PARTICIPATION,
    EIGENVALUES,
    EIGENFREQUENCIES,
    FREQUENCIES,

    BUCKLING_MODE,
    BUCKLING_FACTORS,

    DISPLACEMENT_REAL,
    DISPLACEMENT_IMAG,
    STRESS_REAL,
    STRESS_IMAG,
    STRAIN_REAL,
    STRAIN_IMAG,

    LAMBDA,

    DENS_GRAD,
    DENSITY,
    ORIENTATION_GRAD,
    ORIENTATION,

    COUNT
};

// Stable FEMaster output naming and input-deck request translation
std::string_view              output_field_name(OutputField field);
std::optional<OutputField>    output_field_from_request(std::string token);

// Step/frame lifetime classification
bool output_field_is_step_field(OutputField field);

} // namespace fem::io::writer
