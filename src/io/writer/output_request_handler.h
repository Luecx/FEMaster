/**
 * @file output_request_handler.h
 * @brief Declares lazy recovery and caching of finite-element output fields.
 *
 * Each requested result is identified by its physical quantity and FieldDomain.
 * Analysis steps provide primary fields; the handler recovers derived fields only
 * when an output requires them. Recoveries retain their natural domain and a
 * NODE request projects an ELEMENT_NODAL result only where this is meaningful.
 *
 * The handler owns computed fields for the current frame and borrows provided
 * primary fields. No default requests or dependency registry are maintained:
 * every load case selects its own defaults when the input deck selected none.
 *
 * @see OutputField
 * @see model::Field
 */

#pragma once

#include "output_field.h"

#include "../../core/types_num.h"
#include "../../data/field.h"

#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace fem::model {
struct Model;
struct ModelData;
}

namespace fem::io::writer {

struct ResultWriters;

/**
 * @brief Identifies a result quantity together with its spatial domain.
 *
 * The domain distinguishes nodal values, independent element-nodal values,
 * element totals and element integration/material-point data. It is part of
 * the cache key, so different representations can coexist within one frame.
 */
struct OutputKey {
    OutputField        field  = OutputField::DISPLACEMENT;
    model::FieldDomain domain = model::FieldDomain::NODE;

    bool operator<(const OutputKey& other) const {
        if (field != other.field) return field < other.field;
        return domain < other.domain;
    }
};

/**
 * @brief Stores one borrowed or computed result field for the current frame.
 *
 * A directly provided field remains owned by the analysis. A lazily recovered
 * or projected field is owned by the handler and remains valid until the next
 * begin_frame() call. Provided data takes precedence when both exist.
 */
struct CachedOutputField {
    model::Field*               provided = nullptr;
    std::optional<model::Field> computed;

    model::Field& value() {
        return provided ? *provided : *computed;
    }
};

/**
 * @brief Resolves requested model results and dispatches them to result writers.
 *
 * Each analysis step owns one instance. The parser appends explicit requests
 * while reading the step; the load case adds defaults only if no requests exist.
 * Primary results are provided at each frame, and compute() performs missing
 * numerical recoveries in their native domain. All computed representations
 * are cached until the next frame. Independent spectral/step summary fields
 * are written directly by their owning analysis procedure.
 *
 * Linearized, perturbation and nonlinear recovery use the existing model-state
 * context; these settings do not affect request selection or field ownership.
 */
class OutputRequestHandler {
    // Bound model and requested field/domain combinations
    model::Model*          model_ = nullptr;
    std::vector<OutputKey> requests_;

    // Frame-local fields: provided values are borrowed and computed values owned
    std::map<OutputKey, CachedOutputField> fields_;

    // Current frame metadata
    Precision   frame_value_ = std::numeric_limits<Precision>::quiet_NaN();
    std::string frame_suffix_;

    // Numerical recovery context provided by the owning load case
    bool nonlinear_ = false;
    const model::Field* perturbation_base_ = nullptr;
    const model::Field* linearization_base_ = nullptr;
    const model::Field* linearization_temperature_ = nullptr;
    std::function<void()> before_compute_;

public:
    // Bind the compiled model and configure requests
    void bind(model::Model* model);
    void request(OutputField field, model::FieldDomain domain = model::FieldDomain::NODE);
    void request(std::initializer_list<OutputKey> fields);
    bool empty() const;

    // Select the constitutive recovery configuration
    void set_nonlinear(bool nonlinear);
    void set_perturbation_base(const model::Field* displacement);
    void set_linearization_base(const model::Field* displacement,
                                const model::Field* temperature);
    void set_before_compute(std::function<void()> callback);

    // Supply current-frame primary results; provided fields are not copied
    void begin_frame(Precision frame_value = std::numeric_limits<Precision>::quiet_NaN(),
                     std::string suffix = {});
    void provide(OutputField field, model::Field& value);

    // Return a cached field, recover one natively or project element-nodal data
    model::Field& compute(OutputField field, model::FieldDomain domain);

    // Write all explicitly requested frame fields, using the returned domain
    void write_frame(ResultWriters& writer, const model::ModelData* model_data);

private:
    // Save a recovered field in its actual domain without overriding provided data
    void store(OutputField field, model::Field value);
};

} // namespace fem::io::writer
