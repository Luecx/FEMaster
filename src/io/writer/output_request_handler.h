/**
 * @file output_request_handler.h
 * @brief Declares lazy dependency-based result recovery for one analysis step.
 *
 * Every load case owns one OutputRequestHandler. The load case provides only
 * fields that exist directly in its solved state, while the handler recursively
 * resolves requested derived quantities through the model. Step-lifetime and
 * frame-lifetime storage are separated so modal/buckling summaries can coexist
 * with repeatedly changing frame fields.
 *
 * @see output_field.h
 * @see model::Field
 */

#pragma once

#include "output_field.h"

#include "../../core/types_num.h"
#include "../../data/field.h"

#include <array>
#include <initializer_list>
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
 * @brief Resolves and writes requested result quantities for one analysis step.
 *
 * The handler owns request state and lazily recovered fields. Fields supplied
 * through provide()/provide_step() remain owned by the analysis step and are
 * referenced only until the next frame reset or the end of the step.
 */
class OutputRequestHandler {
    static constexpr std::size_t field_count =
        static_cast<std::size_t>(OutputField::COUNT);

    // Analysis dependency used only when a requested field must be recovered
    model::Model* model_ = nullptr;

    // Request state: defaults remain available so PRESELECT can restore them
    std::array<bool, field_count> defaults_{};
    std::array<bool, field_count> requests_{};
    bool explicit_requests_ = false;

    // Non-owning primary/direct fields and owned lazily recovered frame fields
    std::array<model::Field*, field_count> step_fields_{};
    std::array<model::Field*, field_count> frame_fields_{};
    std::array<std::optional<model::Field>, field_count> computed_{};

    // Recursion guard for accidental cyclic requirement definitions
    std::array<bool, field_count> resolving_{};

    // Current frame metadata used to preserve existing result names/values
    Precision   frame_value_ = std::numeric_limits<Precision>::quiet_NaN();
    std::string frame_suffix_;

    // Kinematic recovery mode selected by the owning analysis step
    bool nonlinear_ = false;

public:
    // Analysis binding and step-defined default requests
    void bind(model::Model* model);
    void set_defaults(std::initializer_list<OutputField> fields);
    void use_defaults();
    void set_nonlinear(bool nonlinear);

    // Explicit requests supplied by the input deck
    void begin_explicit_requests();
    void request(OutputField field);
    bool requested(OutputField field) const;

    // Direct quantities supplied by the active analysis
    void provide_step(OutputField field, model::Field& value);
    void begin_frame(Precision frame_value = std::numeric_limits<Precision>::quiet_NaN(),
                     std::string suffix = {});
    void provide(OutputField field, model::Field& value);

    // Recursive dependency resolution
    std::vector<OutputField> requirements(OutputField field) const;
    model::Field& resolve(OutputField field);

    // Write all requested quantities belonging to the respective lifetime
    void write_frame(ResultWriters& writer, const model::ModelData* model_data);
    void write_step(ResultWriters& writer, const model::ModelData* model_data);

private:
    static std::size_t index(OutputField field);

    // Compute one missing field after all declared requirements are available
    void compute(OutputField field);

    // Empty thermal strain is the explicit "thermal recovery not active" value
    const model::Field* thermal_free_strain();

    // Shared writer path preserving established result names
    void write(OutputField field,
               ResultWriters& writer,
               const model::ModelData* model_data,
               const std::string& suffix,
               Precision frame_value);
};

} // namespace fem::io::writer
