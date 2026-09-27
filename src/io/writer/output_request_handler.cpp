/**
 * @file output_request_handler.cpp
 * @brief Implements recursive field dependency resolution and requested output.
 *
 * The dependency graph is intentionally expressed as direct switch statements.
 * Result quantities are few, domain-specific and stable; keeping requirements
 * and recovery visible in one place makes numerical dependencies easy to audit.
 *
 * A requested field is resolved in the following order:
 *  1. a frame field supplied directly by the active analysis,
 *  2. a step field supplied directly by the active analysis,
 *  3. a field already recovered during the current frame,
 *  4. recursive recovery from its declared requirements.
 *
 * @see output_request_handler.h
 */

#include "output_request_handler.h"

#include "writers.h"

#include "../../core/logging.h"
#include "../../model/model.h"

#include <utility>

namespace fem::io::writer {

/**
 * Binds the model used for all derived-field recovery.
 */
void OutputRequestHandler::bind(model::Model* model) {
    model_ = model;
}

/**
 * Replaces the step's predefined output selection.
 *
 * Defaults are copied to the active request set only while the input deck has
 * not switched to explicit output requests.
 */
void OutputRequestHandler::set_defaults(std::initializer_list<OutputField> fields) {
    defaults_.fill(false);
    for (const OutputField field : fields) defaults_[index(field)] = true;

    if (!explicit_requests_) requests_ = defaults_;
}

/**
 * Restores the owning step's predefined output selection.
 */
void OutputRequestHandler::use_defaults() {
    requests_          = defaults_;
    explicit_requests_ = false;
}

/**
 * Selects nonlinear strain kinematics for stress/strain recovery.
 */
void OutputRequestHandler::set_nonlinear(bool nonlinear) {
    nonlinear_ = nonlinear;
}

/**
 * Installs an optional preparation hook for derived model recovery.
 *
 * Linear procedures leave this empty. Nonlinear procedures use it to reset
 * constitutive trial history before independent post-processing evaluations,
 * preserving the same state semantics as the former direct recovery calls.
 */
void OutputRequestHandler::set_before_compute(std::function<void()> callback) {
    before_compute_ = std::move(callback);
}

/**
 * Starts an explicit request set exactly once.
 *
 * Abaqus uses *OUTPUT, FIELD before NODE/ELEMENT OUTPUT, while CalculiX may
 * start directly with NODE FILE/EL FILE. Calling this method from either path
 * gives both syntaxes the same replacement semantics.
 */
void OutputRequestHandler::begin_explicit_requests() {
    if (explicit_requests_) return;

    requests_.fill(false);
    explicit_requests_ = true;
}

/**
 * Adds one requested quantity. The first explicit request replaces defaults.
 */
void OutputRequestHandler::request(OutputField field) {
    begin_explicit_requests();
    requests_[index(field)] = true;
}

/**
 * Returns whether a field belongs to the active request set.
 */
bool OutputRequestHandler::requested(OutputField field) const {
    return requests_[index(field)];
}

/**
 * Supplies a value whose lifetime covers the complete analysis step.
 */
void OutputRequestHandler::provide_step(OutputField field, model::Field& value) {
    step_fields_[index(field)] = &value;
}

/**
 * Starts a new result frame.
 *
 * Step fields and request state deliberately survive. Only frame-owned direct
 * references and lazily recovered fields are invalidated.
 */
void OutputRequestHandler::begin_frame(Precision frame_value, std::string suffix) {
    frame_fields_.fill(nullptr);
    computed_.fill(std::nullopt);
    resolving_.fill(false);

    frame_value_  = frame_value;
    frame_suffix_ = std::move(suffix);
}

/**
 * Supplies a direct field for the current frame.
 */
void OutputRequestHandler::provide(OutputField field, model::Field& value) {
    frame_fields_[index(field)] = &value;
}

/**
 * Declares the immediate dependencies of one recoverable output field.
 *
 * Requirements may themselves be derived. resolve() therefore walks this list
 * recursively until it reaches fields supplied directly by the analysis step.
 */
std::vector<OutputField> OutputRequestHandler::requirements(OutputField field) const {
    switch (field) {
        case OutputField::MODE_SHAPE:
        case OutputField::BUCKLING_MODE:
            return {OutputField::DISPLACEMENT};

        case OutputField::STRESS:
        case OutputField::STRAIN:
        case OutputField::STRESS_TOP:
        case OutputField::STRESS_BOT:
        case OutputField::SHELL_RESULTANTS:
            return {OutputField::DISPLACEMENT, OutputField::THERMAL_FREE_STRAIN};

        case OutputField::LOCAL_SECTION_FORCES:
        case OutputField::SHEAR_FLOW:
        case OutputField::COMPLIANCE:
        case OutputField::ORIENTATION_GRAD:
            return {OutputField::DISPLACEMENT};

        case OutputField::HEAT_FLUX:
            return {OutputField::TEMPERATURE};

        case OutputField::STRESS_REAL:
        case OutputField::STRAIN_REAL:
            return {OutputField::DISPLACEMENT_REAL, OutputField::THERMAL_FREE_STRAIN};

        case OutputField::STRESS_IMAG:
        case OutputField::STRAIN_IMAG:
            return {OutputField::DISPLACEMENT_IMAG, OutputField::THERMAL_FREE_STRAIN};

        default:
            return {};
    }
}

/**
 * Resolves one field recursively.
 *
 * Aliases such as MODE_SHAPE and BUCKLING_MODE intentionally return the
 * displacement field directly. No duplicate node-sized matrix is created just
 * to attach a different output name; naming is handled only while writing.
 */
model::Field& OutputRequestHandler::resolve(OutputField field) {
    const std::size_t i = index(field);

    if (frame_fields_[i] != nullptr) return *frame_fields_[i];
    if (step_fields_[i]  != nullptr) return *step_fields_[i];
    if (computed_[i].has_value())    return *computed_[i];

    logging::error(model_ != nullptr,
        "OutputRequestHandler: no model is bound");
    logging::error(!resolving_[i],
        "OutputRequestHandler: cyclic output dependency while resolving ",
        output_field_name(field));

    resolving_[i] = true;

    for (const OutputField requirement : requirements(field)) {
        resolve(requirement);
    }

    // Modal and buckling vectors are semantically distinct outputs but are
    // physically the displacement-like eigenvector supplied by the step.
    if (field == OutputField::MODE_SHAPE || field == OutputField::BUCKLING_MODE) {
        resolving_[i] = false;
        return resolve(OutputField::DISPLACEMENT);
    }

    compute(field);
    resolving_[i] = false;

    logging::error(computed_[i].has_value(),
        "OutputRequestHandler: field ", output_field_name(field),
        " is not available for the active analysis step");

    return *computed_[i];
}

/**
 * Writes every requested frame-lifetime quantity.
 */
void OutputRequestHandler::write_frame(ResultWriters& writer, const model::ModelData* model_data) {
    for (std::size_t i = 0; i < field_count; ++i) {
        const auto field = static_cast<OutputField>(i);
        if (!requests_[i] || output_field_is_step_field(field)) continue;

        write(field, writer, model_data, frame_suffix_, frame_value_);
    }
}

/**
 * Writes every requested step-lifetime quantity exactly once.
 */
void OutputRequestHandler::write_step(ResultWriters& writer, const model::ModelData* model_data) {
    for (std::size_t i = 0; i < field_count; ++i) {
        const auto field = static_cast<OutputField>(i);
        if (!requests_[i] || !output_field_is_step_field(field)) continue;

        write(field, writer, model_data, {}, std::numeric_limits<Precision>::quiet_NaN());
    }
}

/**
 * Converts the dense enum value to its array index.
 */
std::size_t OutputRequestHandler::index(OutputField field) {
    return static_cast<std::size_t>(field);
}

/**
 * Computes one missing output field from already resolved dependencies.
 *
 * Recovery routines that naturally produce a pair populate both cache entries
 * at once. A later request for the companion field therefore reuses the same
 * recovery work instead of evaluating the model twice.
 */
void OutputRequestHandler::compute(OutputField field) {
    // Each derived recovery must observe the same accepted physical state.
    // This is a no-op for linear analyses and resets nonlinear trial history
    // when the owning step supplied a preparation callback.
    if (before_compute_) before_compute_();

    switch (field) {
        case OutputField::STRESS:
        case OutputField::STRAIN: {
            auto [stress, strain] = model_->compute_stress_nodal(
                resolve(OutputField::DISPLACEMENT),
                nonlinear_,
                thermal_free_strain());

            computed_[index(OutputField::STRESS)] = std::move(stress);
            computed_[index(OutputField::STRAIN)] = std::move(strain);
            return;
        }

        case OutputField::STRESS_TOP:
        case OutputField::STRESS_BOT: {
            auto [top, bot] = model_->compute_stress_top_bot(
                resolve(OutputField::DISPLACEMENT),
                nonlinear_,
                thermal_free_strain());

            computed_[index(OutputField::STRESS_TOP)] = std::move(top);
            computed_[index(OutputField::STRESS_BOT)] = std::move(bot);
            return;
        }

        case OutputField::SHELL_RESULTANTS:
            computed_[index(field)] = model_->compute_shell_resultants(
                resolve(OutputField::DISPLACEMENT),
                thermal_free_strain());
            return;

        case OutputField::LOCAL_SECTION_FORCES:
            computed_[index(field)] =
                model_->compute_section_forces(resolve(OutputField::DISPLACEMENT));
            return;

        case OutputField::SHEAR_FLOW:
            computed_[index(field)] =
                model_->compute_shear_flow(resolve(OutputField::DISPLACEMENT));
            return;

        case OutputField::COMPLIANCE:
            computed_[index(field)] =
                model_->compute_compliance(resolve(OutputField::DISPLACEMENT));
            return;

        case OutputField::ORIENTATION_GRAD:
            computed_[index(field)] =
                model_->compute_compliance_angle_derivative(resolve(OutputField::DISPLACEMENT));
            return;

        case OutputField::VOLUME:
            computed_[index(field)] = model_->compute_volumes();
            return;

        case OutputField::HEAT_FLUX:
            computed_[index(field)] =
                model_->compute_heat_flux(resolve(OutputField::TEMPERATURE));
            return;

        case OutputField::STRESS_REAL:
        case OutputField::STRAIN_REAL: {
            auto [stress, strain] = model_->compute_stress_nodal(
                resolve(OutputField::DISPLACEMENT_REAL),
                false,
                thermal_free_strain());

            computed_[index(OutputField::STRESS_REAL)] = std::move(stress);
            computed_[index(OutputField::STRAIN_REAL)] = std::move(strain);
            return;
        }

        case OutputField::STRESS_IMAG:
        case OutputField::STRAIN_IMAG: {
            auto [stress, strain] = model_->compute_stress_nodal(
                resolve(OutputField::DISPLACEMENT_IMAG),
                false,
                thermal_free_strain());

            computed_[index(OutputField::STRESS_IMAG)] = std::move(stress);
            computed_[index(OutputField::STRAIN_IMAG)] = std::move(strain);
            return;
        }

        default:
            return;
    }
}

/**
 * Returns the thermal eigendeformation used by structural recovery.
 *
 * Every structural step supplies THERMAL_FREE_STRAIN. Unsupported procedures
 * supply an empty field, which explicitly maps to the existing nullptr meaning
 * of the model recovery API.
 */
const model::Field* OutputRequestHandler::thermal_free_strain() {
    model::Field& thermal = resolve(OutputField::THERMAL_FREE_STRAIN);
    return thermal.rows > 0 ? &thermal : nullptr;
}

/**
 * Resolves and writes one field while preserving the existing FEMaster name.
 *
 * Empty fields are valid dependency placeholders (most notably empty thermal
 * free strain) and are not emitted as result datasets.
 */
void OutputRequestHandler::write(OutputField field,
                                 ResultWriters& writer,
                                 const model::ModelData* model_data,
                                 const std::string& suffix,
                                 Precision frame_value) {
    model::Field& value = resolve(field);
    if (value.rows == 0) return;

    const std::string name = std::string(output_field_name(field)) + suffix;
    const model::ModelData* data =
        value.domain == model::FieldDomain::UNKNOWN ? nullptr : model_data;

    writer.write_field(value, name, data, frame_value);
}

} // namespace fem::io::writer
