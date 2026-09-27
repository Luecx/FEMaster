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
    THERMAL_FREE_STRAIN,

    STRESS,
    STRAIN,
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

/**
 * @brief Returns the unchanged FEMaster result-file base name for a field.
 */
std::string_view output_field_name(OutputField field);

/**
 * @brief Maps an Abaqus/CalculiX or FEMaster request token to an output field.
 *
 * Both compact solver-style names (for example U, RF, S, E) and the existing
 * FEMaster field names are accepted. The mapping does not choose a writer
 * format and therefore treats NODE FILE and NODE OUTPUT identically.
 */
std::optional<OutputField> output_field_from_request(std::string token);

/**
 * @brief Returns true for quantities written once for the complete analysis step.
 */
bool output_field_is_step_field(OutputField field);

} // namespace fem::io::writer
