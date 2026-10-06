/**
 * @file output_field.cpp
 * @brief Implements output-field names and input-request aliases.
 *
 * The translation in this file is intentionally explicit. Output keywords are
 * a small, stable user-facing vocabulary, and a direct switch/table is easier
 * to audit than a generic registration layer.
 *
 * @see output_field.h
 */

#include "output_field.h"

#include <algorithm>
#include <cctype>
#include <unordered_map>

namespace fem::io::writer {

/**
 * Returns the established FEMaster result-file base name for one semantic
 * output quantity. Frame indices are appended by OutputRequestHandler and are
 * intentionally not part of this mapping.
 */
std::string_view output_field_name(OutputField field) {
    switch (field) {
        case OutputField::DISPLACEMENT:         return "DISPLACEMENT";
        case OutputField::VELOCITY:             return "VELOCITY";
        case OutputField::ACCELERATION:         return "ACCELERATION";
        case OutputField::TEMPERATURE:          return "TEMPERATURE";
        case OutputField::EXTERNAL_FORCES:      return "EXTERNAL_FORCES";
        case OutputField::INTERNAL_FORCES:      return "INTERNAL_FORCES";
        case OutputField::REACTION_FORCES:      return "REACTION_FORCES";
        case OutputField::STRESS:               return "STRESS";
        case OutputField::STRAIN:               return "STRAIN";
        case OutputField::PEEQ:                 return "PEEQ";
        case OutputField::STRESS_TOP:           return "STRESS_TOP";
        case OutputField::STRESS_BOT:           return "STRESS_BOT";
        case OutputField::SHELL_RESULTANTS:     return "SHELL_RESULTANTS";
        case OutputField::LOCAL_SECTION_FORCES: return "LOCAL_SECTION_FORCES";
        case OutputField::SHEAR_FLOW:           return "SHEAR_FLOW";
        case OutputField::COMPLIANCE:           return "COMPLIANCE";
        case OutputField::VOLUME:               return "VOLUME";
        case OutputField::HEAT_FLUX:            return "HEAT_FLUX";
        case OutputField::MODE_SHAPE:           return "MODE_SHAPE";
        case OutputField::PARTICIPATION:        return "PARTICIPATION";
        case OutputField::EIGENVALUES:          return "EIGENVALUES";
        case OutputField::EIGENFREQUENCIES:     return "EIGENFREQUENCIES";
        case OutputField::FREQUENCIES:          return "FREQUENCIES";
        case OutputField::BUCKLING_MODE:        return "BUCKLING_MODE";
        case OutputField::BUCKLING_FACTORS:     return "BUCKLING_FACTORS";
        case OutputField::DISPLACEMENT_REAL:    return "DISPLACEMENT_REAL";
        case OutputField::DISPLACEMENT_IMAG:    return "DISPLACEMENT_IMAG";
        case OutputField::STRESS_REAL:          return "STRESS_REAL";
        case OutputField::STRESS_IMAG:          return "STRESS_IMAG";
        case OutputField::STRAIN_REAL:          return "STRAIN_REAL";
        case OutputField::STRAIN_IMAG:          return "STRAIN_IMAG";
        case OutputField::LAMBDA:               return "LAMBDA";
        case OutputField::DENS_GRAD:            return "DENS_GRAD";
        case OutputField::DENSITY:              return "DENSITY";
        case OutputField::ORIENTATION_GRAD:     return "ORIENTATION_GRAD";
        case OutputField::ORIENTATION:          return "ORIENTATION";
        case OutputField::COUNT:                break;
    }
    return "UNKNOWN";
}

/**
 * Normalizes one input-deck output variable.
 *
 * Abaqus and CalculiX share the important compact field identifiers used here.
 * FEMaster's long result names are accepted as aliases as well so native decks
 * do not need a second output vocabulary.
 */
std::optional<OutputField> output_field_from_request(std::string token) {
    std::transform(token.begin(), token.end(), token.begin(),
        [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    static const std::unordered_map<std::string, OutputField> fields = {
        {"U",                    OutputField::DISPLACEMENT},
        {"V",                    OutputField::VELOCITY},
        {"A",                    OutputField::ACCELERATION},
        {"NT",                   OutputField::TEMPERATURE},
        {"NT11",                 OutputField::TEMPERATURE},
        {"TEMP",                 OutputField::TEMPERATURE},
        {"CF",                   OutputField::EXTERNAL_FORCES},
        {"NFORC",                OutputField::INTERNAL_FORCES},
        {"RF",                   OutputField::REACTION_FORCES},
        {"S",                    OutputField::STRESS},
        {"E",                    OutputField::STRAIN},
        {"PEEQ",                 OutputField::PEEQ},
        {"SF",                   OutputField::LOCAL_SECTION_FORCES},
        {"HFL",                  OutputField::HEAT_FLUX},

        {"DISPLACEMENT",         OutputField::DISPLACEMENT},
        {"VELOCITY",             OutputField::VELOCITY},
        {"ACCELERATION",         OutputField::ACCELERATION},
        {"TEMPERATURE",          OutputField::TEMPERATURE},
        {"EXTERNAL_FORCES",      OutputField::EXTERNAL_FORCES},
        {"INTERNAL_FORCES",      OutputField::INTERNAL_FORCES},
        {"REACTION_FORCES",      OutputField::REACTION_FORCES},
        {"STRESS",               OutputField::STRESS},
        {"STRAIN",               OutputField::STRAIN},
        {"STRESS_TOP",           OutputField::STRESS_TOP},
        {"STRESS_BOT",           OutputField::STRESS_BOT},
        {"SHELL_RESULTANTS",     OutputField::SHELL_RESULTANTS},
        {"LOCAL_SECTION_FORCES", OutputField::LOCAL_SECTION_FORCES},
        {"SHEAR_FLOW",           OutputField::SHEAR_FLOW},
        {"COMPLIANCE",           OutputField::COMPLIANCE},
        {"VOLUME",               OutputField::VOLUME},
        {"HEAT_FLUX",            OutputField::HEAT_FLUX},
        {"MODE_SHAPE",           OutputField::MODE_SHAPE},
        {"PARTICIPATION",        OutputField::PARTICIPATION},
        {"EIGENVALUES",          OutputField::EIGENVALUES},
        {"EIGENFREQUENCIES",     OutputField::EIGENFREQUENCIES},
        {"FREQUENCIES",          OutputField::FREQUENCIES},
        {"BUCKLING_MODE",        OutputField::BUCKLING_MODE},
        {"BUCKLING_FACTORS",     OutputField::BUCKLING_FACTORS},
        {"DISPLACEMENT_REAL",    OutputField::DISPLACEMENT_REAL},
        {"DISPLACEMENT_IMAG",    OutputField::DISPLACEMENT_IMAG},
        {"STRESS_REAL",          OutputField::STRESS_REAL},
        {"STRESS_IMAG",          OutputField::STRESS_IMAG},
        {"STRAIN_REAL",          OutputField::STRAIN_REAL},
        {"STRAIN_IMAG",          OutputField::STRAIN_IMAG},
        {"LAMBDA",               OutputField::LAMBDA},
        {"DENS_GRAD",            OutputField::DENS_GRAD},
        {"DENSITY",              OutputField::DENSITY},
        {"ORIENTATION_GRAD",     OutputField::ORIENTATION_GRAD},
        {"ORIENTATION",          OutputField::ORIENTATION}
    };

    const auto it = fields.find(token);
    if (it == fields.end()) return std::nullopt;
    return it->second;
}

/**
 * Step fields have one value for the complete analysis rather than one value
 * per result frame. Everything else follows frame lifetime.
 */
bool output_field_is_step_field(OutputField field) {
    switch (field) {
        case OutputField::PARTICIPATION:
        case OutputField::EIGENVALUES:
        case OutputField::EIGENFREQUENCIES:
        case OutputField::FREQUENCIES:
        case OutputField::BUCKLING_FACTORS:
            return true;
        default:
            return false;
    }
}

} // namespace fem::io::writer
