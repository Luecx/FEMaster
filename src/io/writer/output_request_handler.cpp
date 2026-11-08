/**
 * @file output_request_handler.cpp
 * @brief Implements domain-aware lazy output recovery and frame-local caching.
 *
 * The handler first uses an exact (quantity, domain) match from its cache.
 * Selected nodal requests may recover the element-nodal representation and
 * project it with Field::to_nodal(). Otherwise a direct switch evaluates the
 * requested physical quantity in its natural domain, recursively obtaining
 * any required primary fields. Recoveries that return multiple fields store
 * them together so subsequent requests reuse the same model evaluation.
 *
 * @see output_request_handler.h
 */

#include "output_request_handler.h"

#include "writers.h"

#include "../../core/logging.h"
#include "../../model/model.h"

#include <algorithm>
#include <utility>

namespace fem::io::writer {

void OutputRequestHandler::bind(model::Model* model) {
    model_ = model;
}

/**
 * Adds one explicit field/domain request, retaining the order of first use.
 *
 * Duplicate requests do not cause duplicate output blocks.
 */
void OutputRequestHandler::request(OutputField field, model::FieldDomain domain) {
    const auto already_requested = std::find_if(
        requests_.begin(), requests_.end(),
        [&](const OutputKey& key) { return key.field == field && key.domain == domain; });

    if (already_requested == requests_.end()) {
        requests_.push_back({field, domain});
    }
}

void OutputRequestHandler::request(std::initializer_list<OutputKey> fields) {
    for (const OutputKey& key : fields) {
        request(key.field, key.domain);
    }
}

bool OutputRequestHandler::empty() const {
    return requests_.empty();
}

void OutputRequestHandler::set_nonlinear(bool nonlinear) {
    nonlinear_ = nonlinear;
}

void OutputRequestHandler::set_perturbation_base(const model::Field* displacement) {
    perturbation_base_ = displacement;
}

void OutputRequestHandler::set_linearization_base(const model::Field* displacement,
                                                  const model::Field* temperature) {
    linearization_base_ = displacement;
    linearization_temperature_ = temperature;
}

void OutputRequestHandler::set_before_compute(std::function<void()> callback) {
    before_compute_ = std::move(callback);
}

/**
 * Starts a new frame with an empty cache and fresh primary-field references.
 *
 * Requests persist between frames. Numeric recovery settings are scoped to
 * the owning analysis and persist for its complete run.
 */
void OutputRequestHandler::begin_frame(Precision frame_value, std::string suffix) {
    fields_.clear();
    frame_value_  = frame_value;
    frame_suffix_ = std::move(suffix);
}

/**
 * Supplies one analysis-owned field without copying its numerical values.
 *
 * The field's own domain selects its cache key. Provided data always wins
 * over an existing computed representation of the same quantity and domain.
 */
void OutputRequestHandler::provide(OutputField field, model::Field& value) {
    fields_[{field, value.domain}].provided = &value;
}

/**
 * Stores a result in its native domain, preserving directly supplied data.
 *
 * Pair recoveries such as stress/strain can populate both cache entries
 * without overwriting a primary result already provided by the analysis.
 */
void OutputRequestHandler::store(OutputField field, model::Field value) {
    auto& entry = fields_[{field, value.domain}];
    if (!entry.provided) {
        entry.computed = std::move(value);
    }
}

/**
 * Returns an existing result or computes the requested physical representation.
 *
 * Exact cache hits are returned immediately. For eligible NODE requests, the
 * ELEMENT_NODAL representation is computed first and projected only once.
 * All other requests are evaluated in their native domain by the explicit
 * numerical switch below. The switch itself expresses every dependency
 * through recursive compute() calls; no separate dependency graph is needed.
 *
 * Constitutive results use the analysis-supplied linearization/perturbation
 * context. A supplied before_compute hook restores nonlinear trial history
 * before each independent material evaluation.
 */
model::Field& OutputRequestHandler::compute(OutputField field, model::FieldDomain domain) {
    const OutputKey key{field, domain};

    // Reuse directly supplied and already computed values
    if (auto it = fields_.find(key); it != fields_.end()) {
        if (it->second.provided || it->second.computed.has_value()) {
            return it->second.value();
        }
    }

    logging::error(model_ != nullptr,
        "OutputRequestHandler: no model is bound");

    // Only quantities with meaningful element-nodal recovery can be averaged
    if (domain == model::FieldDomain::NODE) {
        bool project_element_nodal = false;
        switch (field) {
            case OutputField::STRESS:
            case OutputField::STRAIN:
            case OutputField::PEEQ:
            case OutputField::STRESS_TOP:
            case OutputField::STRESS_BOT:
            case OutputField::SHELL_RESULTANTS:
            case OutputField::LOCAL_SECTION_FORCES:
            case OutputField::SHEAR_FLOW:
            case OutputField::HEAT_FLUX:
            case OutputField::STRESS_REAL:
            case OutputField::STRAIN_REAL:
            case OutputField::STRESS_IMAG:
            case OutputField::STRAIN_IMAG:
                project_element_nodal = true;
                break;
            default:
                break;
        }

        if (project_element_nodal) {
            auto& source = compute(field, model::FieldDomain::ELEMENT_NODAL);
            store(field, source.to_nodal(*model_->_data, std::string(output_field_name(field))));
            return fields_.at(key).value();
        }
    }

    // Mode-shape labels are aliases of the displacement field of their frame
    if (field == OutputField::MODE_SHAPE || field == OutputField::BUCKLING_MODE) {
        logging::error(domain == model::FieldDomain::NODE,
            "OutputRequestHandler: mode shapes require NODE output");
        return compute(OutputField::DISPLACEMENT, model::FieldDomain::NODE);
    }

    // Enforce the actual native domain instead of silently substituting one
    switch (field) {
        case OutputField::STRESS:
        case OutputField::STRAIN:
        case OutputField::PEEQ:
        case OutputField::STRESS_TOP:
        case OutputField::STRESS_BOT:
        case OutputField::SHELL_RESULTANTS:
        case OutputField::LOCAL_SECTION_FORCES:
        case OutputField::SHEAR_FLOW:
        case OutputField::HEAT_FLUX:
        case OutputField::STRESS_REAL:
        case OutputField::STRAIN_REAL:
        case OutputField::STRESS_IMAG:
        case OutputField::STRAIN_IMAG:
            logging::error(domain == model::FieldDomain::ELEMENT_NODAL,
                "OutputRequestHandler: field ", output_field_name(field),
                " requires ELEMENT_NODAL recovery");
            break;

        case OutputField::COMPLIANCE:
        case OutputField::VOLUME:
        case OutputField::ORIENTATION_GRAD:
            logging::error(domain == model::FieldDomain::ELEMENT,
                "OutputRequestHandler: field ", output_field_name(field),
                " requires ELEMENT output");
            break;

        default:
            logging::error(false,
                "OutputRequestHandler: field ", output_field_name(field),
                " must be provided by the analysis in the requested domain");
    }

    // Evaluate from the same accepted material state in nonlinear analyses
    if (before_compute_) before_compute_();

    const model::Field* target_temperature = model_->_data->temperature_target
        ? model_->_data->temperature_target.get()
        : model_->_data->temperature.get();
    const model::Field* base_temperature = linearization_base_
        ? linearization_temperature_ : model_->_data->temperature.get();

    // Linear perturbation fields are increments relative to the committed base
    if (perturbation_base_) {
        auto& increment = compute(OutputField::DISPLACEMENT, model::FieldDomain::NODE);
        model::Field base  = *perturbation_base_;
        model::Field total = base;
        total += increment;

        switch (field) {
            case OutputField::STRESS:
            case OutputField::STRAIN: {
                auto [stress, strain] = model_->compute_stress(
                    total, target_temperature, &base, base_temperature);
                if (before_compute_) before_compute_();
                auto [base_stress, base_strain] = model_->compute_stress(
                    base, base_temperature, &base, base_temperature);

                stress -= base_stress;
                strain -= base_strain;
                store(OutputField::STRESS, std::move(stress));
                store(OutputField::STRAIN, std::move(strain));
                return fields_.at(key).value();
            }

            case OutputField::STRESS_TOP:
            case OutputField::STRESS_BOT: {
                auto [top, bot] = model_->compute_stress_top_bot(
                    total, target_temperature, &base, base_temperature);
                if (before_compute_) before_compute_();
                auto [base_top, base_bot] = model_->compute_stress_top_bot(
                    base, base_temperature, &base, base_temperature);

                top -= base_top;
                bot -= base_bot;
                store(OutputField::STRESS_TOP, std::move(top));
                store(OutputField::STRESS_BOT, std::move(bot));
                return fields_.at(key).value();
            }

            case OutputField::SHELL_RESULTANTS: {
                auto resultants = model_->compute_shell_resultants(
                    total, target_temperature, &base, base_temperature);
                if (before_compute_) before_compute_();
                resultants -= model_->compute_shell_resultants(
                    base, base_temperature, &base, base_temperature);
                store(field, std::move(resultants));
                return fields_.at(key).value();
            }

            case OutputField::LOCAL_SECTION_FORCES: {
                auto forces = model_->compute_section_forces(
                    total, target_temperature, &base, base_temperature);
                if (before_compute_) before_compute_();
                forces -= model_->compute_section_forces(
                    base, base_temperature, &base, base_temperature);
                store(field, std::move(forces));
                return fields_.at(key).value();
            }

            case OutputField::SHEAR_FLOW: {
                auto flow = model_->compute_shear_flow(total);
                flow -= model_->compute_shear_flow(base);
                store(field, std::move(flow));
                return fields_.at(key).value();
            }

            case OutputField::PEEQ: {
                auto peeq = model_->compute_peeq_element_nodal();
                peeq.set_zero();
                store(field, std::move(peeq));
                return fields_.at(key).value();
            }

            case OutputField::COMPLIANCE:
            case OutputField::ORIENTATION_GRAD:
                logging::error(false,
                    "Linear perturbation does not support compliance or orientation gradients");
                break;

            default:
                break;
        }
    }

    // Recover the requested quantity in the native model-defined domain
    switch (field) {
        case OutputField::STRESS:
        case OutputField::STRAIN: {
            auto& displacement = compute(OutputField::DISPLACEMENT, model::FieldDomain::NODE);
            auto [stress, strain] = model_->compute_stress(
                displacement,
                target_temperature,
                linearization_base_ ? linearization_base_ : (nonlinear_ ? &displacement : nullptr),
                linearization_base_ ? base_temperature : (nonlinear_ ? target_temperature : nullptr));

            store(OutputField::STRESS, std::move(stress));
            store(OutputField::STRAIN, std::move(strain));
            break;
        }

        case OutputField::STRESS_TOP:
        case OutputField::STRESS_BOT: {
            auto& displacement = compute(OutputField::DISPLACEMENT, model::FieldDomain::NODE);
            auto [top, bot] = model_->compute_stress_top_bot(
                displacement,
                target_temperature,
                linearization_base_ ? linearization_base_ : (nonlinear_ ? &displacement : nullptr),
                linearization_base_ ? base_temperature : (nonlinear_ ? target_temperature : nullptr));

            store(OutputField::STRESS_TOP, std::move(top));
            store(OutputField::STRESS_BOT, std::move(bot));
            break;
        }

        case OutputField::SHELL_RESULTANTS: {
            auto& displacement = compute(OutputField::DISPLACEMENT, model::FieldDomain::NODE);
            store(field, model_->compute_shell_resultants(
                displacement,
                target_temperature,
                linearization_base_ ? linearization_base_ : (nonlinear_ ? &displacement : nullptr),
                linearization_base_ ? base_temperature : (nonlinear_ ? target_temperature : nullptr)));
            break;
        }

        case OutputField::LOCAL_SECTION_FORCES: {
            auto& displacement = compute(OutputField::DISPLACEMENT, model::FieldDomain::NODE);
            store(field, model_->compute_section_forces(
                displacement,
                target_temperature,
                linearization_base_ ? linearization_base_ : (nonlinear_ ? &displacement : nullptr),
                linearization_base_ ? base_temperature : (nonlinear_ ? target_temperature : nullptr)));
            break;
        }

        case OutputField::SHEAR_FLOW:
            store(field, model_->compute_shear_flow(
                compute(OutputField::DISPLACEMENT, model::FieldDomain::NODE)));
            break;

        case OutputField::COMPLIANCE:
            store(field, model_->compute_compliance(
                compute(OutputField::DISPLACEMENT, model::FieldDomain::NODE)));
            break;

        case OutputField::ORIENTATION_GRAD:
            store(field, model_->compute_compliance_angle_derivative(
                compute(OutputField::DISPLACEMENT, model::FieldDomain::NODE)));
            break;

        case OutputField::PEEQ:
            store(field, model_->compute_peeq_element_nodal());
            break;

        case OutputField::VOLUME:
            store(field, model_->compute_volumes());
            break;

        case OutputField::HEAT_FLUX:
            store(field, model_->compute_heat_flux(
                compute(OutputField::TEMPERATURE, model::FieldDomain::NODE)));
            break;

        case OutputField::STRESS_REAL:
        case OutputField::STRAIN_REAL: {
            auto [stress, strain] = model_->compute_stress(
                compute(OutputField::DISPLACEMENT_REAL, model::FieldDomain::NODE),
                target_temperature, nullptr, nullptr);
            store(OutputField::STRESS_REAL, std::move(stress));
            store(OutputField::STRAIN_REAL, std::move(strain));
            break;
        }

        case OutputField::STRESS_IMAG:
        case OutputField::STRAIN_IMAG: {
            auto [stress, strain] = model_->compute_stress(
                compute(OutputField::DISPLACEMENT_IMAG, model::FieldDomain::NODE),
                target_temperature, nullptr, nullptr);
            store(OutputField::STRESS_IMAG, std::move(stress));
            store(OutputField::STRAIN_IMAG, std::move(strain));
            break;
        }

        default:
            break;
    }

    // Native recoveries must produce exactly the representation they declare
    const auto it = fields_.find(key);
    logging::error(it != fields_.end() && (it->second.provided || it->second.computed.has_value()),
        "OutputRequestHandler: recovery produced no ", output_field_name(field), " field");
    return it->second.value();
}

/**
 * Writes the requested frame fields in their explicitly selected domains.
 *
 * The ELEMENT_NODAL name prefix preserves existing RES/FRD conventions.
 * A field in any other domain keeps its semantic base name.
 */
void OutputRequestHandler::write_frame(ResultWriters& writer, const model::ModelData* model_data) {
    for (const OutputKey& key : requests_) {
        if (output_field_is_step_field(key.field)) continue;

        model::Field& value = compute(key.field, key.domain);
        if (value.rows == 0) continue;

        const std::string prefix = key.domain == model::FieldDomain::ELEMENT_NODAL
            ? "ELEMENT_NODAL_" : "";
        const std::string name = prefix + std::string(output_field_name(key.field)) + frame_suffix_;
        const model::ModelData* data = value.domain == model::FieldDomain::UNKNOWN
            ? nullptr : model_data;

        writer.write_field(value, name, data, frame_value_);
    }
}

} // namespace fem::io::writer
