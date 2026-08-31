#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text()


def write(path, text):
    (ROOT / path).write_text(text)


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one occurrence, found {count}")
    return text.replace(old, new, 1)


def replace_count(text, old, new, expected, label):
    count = text.count(old)
    if count != expected:
        raise RuntimeError(f"{label}: expected {expected} occurrences, found {count}")
    return text.replace(old, new)


def block(text, start_marker, end_marker):
    start = text.index(start_marker)
    end = text.index(end_marker, start)
    return start, end, text[start:end]


# -----------------------------------------------------------------------------
# J2: keep frozen-history evaluate() separate from source->target integrate().
# -----------------------------------------------------------------------------
path = "src/material/isotropic_j2_state.cpp"
text = read(path)

start = text.index("void copy_source_to_target(")
end = text.index("\n} // namespace\n", start)
text = text[:start] + text[end:]

marker = "// -----------------------------------------------------------------------------\n// Source -> target state integration\n// -----------------------------------------------------------------------------\n"
start = text.index(marker)
end = text.rindex("\n} // namespace fem::material")
text = text[:start] + text[end:]
write(path, text)

path = "src/material/isotropic_j2_elasticity.cpp"
text = read(path)

store_fn = """void store_state(Precision* state, const State& values) {
    for (Index i = 0; i < state_count; ++i) {
        state[i] = values[static_cast<std::size_t>(i)];
    }
}
"""
validate_fn = store_fn + """
void validate_state_target(const Precision* state, Precision* target_state) {
    logging::error(state != nullptr, "J2: source material state is null");
    logging::error(target_state != nullptr, "J2: target material state is null");
    logging::error(state != target_state,
                   "J2: source and target material state must not alias");
}
"""
text = replace_once(text, store_fn, validate_fn, "J2 target validation")

signatures = [
("""void IsotropicJ2Elasticity::evaluate(const VolumeStrainLinearized& strain,
                                     Precision*                    state,
                                     VolumeStressCauchy&           stress,
                                     Mat6&                         tangent) const {""",
 """void IsotropicJ2Elasticity::integrate(const VolumeStrainLinearized& strain,
                                      const Precision*              state,
                                      Precision*                    target_state,
                                      VolumeStressCauchy&           stress,
                                      Mat6&                         tangent) const {"""),
("""void IsotropicJ2Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     Precision*                       state,
                                     VolumeStressPK2&                 stress,
                                     Mat6&                            tangent) const {""",
 """void IsotropicJ2Elasticity::integrate(const VolumeStrainGreenLagrange& strain,
                                      const Precision*                   state,
                                      Precision*                       target_state,
                                      VolumeStressPK2&                 stress,
                                      Mat6&                            tangent) const {"""),
("""void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                     Precision*                           state,
                                     ShellMaterialStressCauchy&            stress,
                                     Mat5&                                 tangent) const {""",
 """void IsotropicJ2Elasticity::integrate(const ShellMaterialStrainLinearized& strain,
                                      const Precision*                     state,
                                      Precision*                           target_state,
                                      ShellMaterialStressCauchy&            stress,
                                      Mat5&                                 tangent) const {"""),
("""void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                     Precision*                              state,
                                     ShellMaterialStressPK2&                 stress,
                                     Mat5&                                   tangent) const {""",
 """void IsotropicJ2Elasticity::integrate(const ShellMaterialStrainGreenLagrange& strain,
                                      const Precision*                          state,
                                      Precision*                              target_state,
                                      ShellMaterialStressPK2&                 stress,
                                      Mat5&                                   tangent) const {"""),
("""void IsotropicJ2Elasticity::evaluate(const AxialStrainLinearized& strain,
                                     Precision*                   state,
                                     AxialStressCauchy&           stress,
                                     Precision&                   tangent) const {""",
 """void IsotropicJ2Elasticity::integrate(const AxialStrainLinearized& strain,
                                      const Precision*             state,
                                      Precision*                   target_state,
                                      AxialStressCauchy&           stress,
                                      Precision&                   tangent) const {"""),
("""void IsotropicJ2Elasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                     Precision*                      state,
                                     AxialStressPK2&                 stress,
                                     Precision&                      tangent) const {""",
 """void IsotropicJ2Elasticity::integrate(const AxialStrainGreenLagrange& strain,
                                      const Precision*                  state,
                                      Precision*                      target_state,
                                      AxialStressPK2&                 stress,
                                      Precision&                      tangent) const {"""),
]
for i, (old, new) in enumerate(signatures):
    text = replace_once(text, old, new, f"J2 integrate signature {i}")

text = replace_count(
    text,
    "    const State committed = load_state(state);",
    "    validate_state_target(state, target_state);\n    const State committed = load_state(state);",
    6,
    "J2 source target validation calls"
)
text = replace_count(
    text,
    "    store_state(state, update.state);",
    "    store_state(target_state, update.state);",
    6,
    "J2 target writes"
)
write(path, text)


# -----------------------------------------------------------------------------
# Solid elements: central source/target row access and explicit integration.
# -----------------------------------------------------------------------------
path = "src/model/solid/element_solid.h"
text = read(path)
text = replace_once(
    text,
    "    Mat6 material_tangent_reference(Precision r, Precision s, Precision t, Precision* state);",
    "    const Precision* material_state_source(Index ip);\n"
    "    Precision*       material_state_target(Index ip);\n\n"
    "    Mat6 material_tangent_reference(Precision r, Precision s, Precision t, const Precision* state);",
    "solid state helper declarations"
)
text = replace_count(text, "                           Precision*                    state,", "                           const Precision*              state,", 1, "solid linear evaluate const state")
text = replace_count(text, "                           Precision*                       state,", "                           const Precision*                 state,", 1, "solid finite evaluate const state")
write(path, text)

path = "src/model/solid/element_solid_state.ipp"
text = read(path)
insert_marker = "namespace fem::model {\n\n"
helpers = """namespace fem::model {

/**
 * Returns the accepted material-history row for one solid integration point.
 * Stateless materials use no history storage and therefore return `nullptr`.
 */
template<Index N>
const Precision* SolidElement<N>::material_state_source(Index ip) {
    auto* section = get_section();
    logging::error(section->material_ && section->material_->has_elasticity(),
                   "SolidElement: material elasticity is not available for element ", this->elem_id);

    const auto elasticity = section->material_->elasticity();
    if (elasticity->state_size() == 0) return nullptr;

    logging::error(this->_model_data && this->_model_data->material_state,
                   "SolidElement: stateful material requires accepted source history for element ", this->elem_id);
    logging::error(this->_model_data->material_state->components >= elasticity->state_size(),
                   "SolidElement: source material-state field is too narrow for element ", this->elem_id);

    const Index row = this->mp_index(ip);
    logging::error(row < this->_model_data->material_state->rows,
                   "SolidElement: source material-state row is out of range for element ", this->elem_id);
    return &(*this->_model_data->material_state)(row, 0);
}

/**
 * Returns the working target-history row for one solid integration point.
 * The nonlinear state manager owns this field and promotes it only after an
 * accepted physical increment.
 */
template<Index N>
Precision* SolidElement<N>::material_state_target(Index ip) {
    auto* section = get_section();
    logging::error(section->material_ && section->material_->has_elasticity(),
                   "SolidElement: material elasticity is not available for element ", this->elem_id);

    const auto elasticity = section->material_->elasticity();
    if (elasticity->state_size() == 0) return nullptr;

    logging::error(this->_model_data && this->_model_data->material_state_target,
                   "SolidElement: stateful material requires target history for element ", this->elem_id);
    logging::error(this->_model_data->material_state_target->components >= elasticity->state_size(),
                   "SolidElement: target material-state field is too narrow for element ", this->elem_id);

    const Index row = this->mp_index(ip);
    logging::error(row < this->_model_data->material_state_target->rows,
                   "SolidElement: target material-state row is out of range for element ", this->elem_id);
    return &(*this->_model_data->material_state_target)(row, 0);
}

"""
text = replace_once(text, insert_marker, helpers, "solid state helper definitions")
write(path, text)

path = "src/model/solid/element_solid.ipp"
text = read(path)
text = replace_count(text, "                                        Precision*                    state,", "                                        const Precision*              state,", 1, "solid linear evaluate definition")
text = replace_count(text, "                                        Precision*                       state,", "                                        const Precision*                 state,", 1, "solid finite evaluate definition")
text = replace_once(
    text,
    "auto SolidElement<N>::material_tangent_reference(Precision r, Precision s, Precision t, Precision* state)",
    "auto SolidElement<N>::material_tangent_reference(Precision r, Precision s, Precision t, const Precision* state)",
    "solid material tangent const state"
)
text = replace_once(
    text,
    "            Precision* state = &(*this->_model_data->material_state)(this->mp_index(ip++), 0);\n            evaluate_material(r, s, t, strain, state, stress, C);",
    "            const Precision* state = this->material_state_source(ip++);\n            evaluate_material(r, s, t, strain, state, stress, C);",
    "solid state-neutral stiffness"
)
write(path, text)

path = "src/model/solid/element_solid_compute.ipp"
text = read(path)

# Recovery is strictly read-only.
text = replace_once(
    text,
    "        Precision* state = &(*this->_model_data->material_state)(this->mp_index(state_ip), 0);",
    "        const Precision* state = this->material_state_source(state_ip);",
    "solid output source state"
)

# compute_stress_state: only the nonlinear path advances history.
start, end, fn = block(text, "void SolidElement<N>::compute_stress_state", "/**\n * Assembles material tangent")
fn = replace_once(
    fn,
    "        Precision*      state = &(*this->_model_data->material_state)(this->mp_index(static_cast<Index>(n)), 0);",
    "        const Index material_ip = static_cast<Index>(n);\n"
    "        const Precision* state = this->material_state_source(material_ip);\n"
    "        Precision* target_state = use_green_lagrange_nl\n"
    "            ? this->material_state_target(material_ip)\n"
    "            : nullptr;",
    "solid stress-state source target"
)
fn = replace_once(
    fn,
    "        evaluate_material(r, s, t, green_lagrange, state, second_pk, tangent);",
    "        integrate_material(r, s, t, green_lagrange, state, target_state, second_pk, tangent);",
    "solid nonlinear stress-state integration"
)
text = text[:start] + fn + text[end:]

# Combined nonlinear tangent: integrate exactly once per point.
start, end, fn = block(text, "MapMatrix SolidElement<N>::stiffness_tangent", "/**\n * Assembles the Total-Lagrangian internal force")
fn = replace_once(
    fn,
    "        Precision* state = &(*this->_model_data->material_state)(this->mp_index(ip), 0);",
    "        const Precision* state = this->material_state_source(ip);\n"
    "        Precision* target_state = this->material_state_target(ip);",
    "solid tangent source target"
)
fn = replace_once(
    fn,
    "        evaluate_material(point.r, point.s, point.t, strain, state, stress, material_tangent);",
    "        integrate_material(point.r, point.s, point.t, strain, state, target_state, stress, material_tangent);",
    "solid tangent integration"
)
text = text[:start] + fn + text[end:]

# Sensitivities inspect accepted history only.
text = text.replace(
    "        Precision* state = &(*this->_model_data->material_state)(this->mp_index(n), 0);",
    "        const Precision* state = this->material_state_source(n);"
)
if "_model_data->material_state)(" in text:
    raise RuntimeError("solid compute still contains direct material-state row access")
write(path, text)


# -----------------------------------------------------------------------------
# Truss residual-only path must integrate just like tangent assembly.
# -----------------------------------------------------------------------------
path = "src/model/truss/truss.cpp"
text = read(path)
start, end, fn = block(text, "void T3::compute_stress_state", "/**\n * Assembles Total-Lagrangian axial internal force")
fn = replace_once(
    fn,
    "        elasticity->evaluate(strain, source_state(*this, elasticity), stress, tangent);",
    "        elasticity->integrate(\n"
    "            strain,\n"
    "            source_state(*this, elasticity),\n"
    "            target_state(*this, elasticity),\n"
    "            stress,\n"
    "            tangent\n"
    "        );",
    "truss nonlinear stress-state integration"
)
text = text[:start] + fn + text[end:]
write(path, text)


# -----------------------------------------------------------------------------
# Shell section interface: read-only evaluate + explicit integrate.
# -----------------------------------------------------------------------------
path = "src/section/section_shell.h"
text = read(path)
text = replace_once(text, "        Precision*                    material_state,\n        Index                         material_state_stride,", "        const Precision*              material_state,\n        Index                         material_state_stride,", "shell evaluate const source")
needle = """    ) const = 0;

    // Recover generalized resultants for output."""
integrate_decl = """    ) const = 0;

    // Advance constitutive history from an immutable source block into a
    // separate target block. Stateful integrated sections override this path;
    // result recovery never calls it.
    virtual void integrate(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        const Precision*              material_state,
        Precision*                    material_state_target,
        Index                         material_state_stride,
        bool                          use_green_lagrange,
        ShellStressResultants&        resultants_shell,
        Mat8&                         tangent_shell
    ) const = 0;

    // Recover generalized resultants for output."""
text = replace_once(text, needle, integrate_decl, "shell integrate declaration")
text = replace_count(text, "        Precision*                    material_state,", "        const Precision*              material_state,", 2, "shell output const state")
write(path, text)

path = "src/section/section_shell.cpp"
text = read(path)
text = replace_once(text, "    Precision*                    material_state,", "    const Precision*              material_state,", "shell output resultants const state")
write(path, text)


# -----------------------------------------------------------------------------
# ABD shell section is section-stateless: integrate delegates to evaluate.
# -----------------------------------------------------------------------------
path = "src/section/section_shell_abd.h"
text = read(path)
text = replace_once(text, "        Precision*                    material_state,", "        const Precision*              material_state,", "ABD evaluate const state")
needle = """    ) const override;

    // Reconstruct one equivalent homogeneous-layer stress distribution"""
insert = """    ) const override;

    // The prescribed ABD law has no constitutive history. Integration therefore
    // produces the same response as read-only evaluation and ignores the target
    // state block.
    void integrate(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        const Precision*              material_state,
        Precision*                    material_state_target,
        Index                         material_state_stride,
        bool                          use_green_lagrange,
        ShellStressResultants&        resultants_shell,
        Mat8&                         tangent_shell
    ) const override;

    // Reconstruct one equivalent homogeneous-layer stress distribution"""
text = replace_once(text, needle, insert, "ABD integrate declaration")
text = replace_once(text, "        Precision*                    material_state,", "        const Precision*              material_state,", "ABD output const state")
write(path, text)

path = "src/section/section_shell_abd.cpp"
text = read(path)
text = replace_once(text, "    Precision*                    material_state,", "    const Precision*              material_state,", "ABD evaluate definition const")
start, end, eval_fn = block(text, "void ABDShellSection::evaluate(", "/**\n * Reconstructs an equivalent physical Cauchy stress")
integrate_fn = """
/**
 * Applies the prescribed ABD response without changing constitutive history.
 *
 * ABD matrices are already a complete section law and contain no material-point
 * internal variables. The target pointer is therefore intentionally unused.
 */
void ABDShellSection::integrate(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_shell,
    const Precision*              material_state,
    Precision*                    material_state_target,
    Index                         material_state_stride,
    bool                          use_green_lagrange,
    ShellStressResultants&        resultants_shell,
    Mat8&                         tangent_shell
) const {
    (void) material_state_target;
    evaluate(
        position_reference,
        shell_basis_global,
        strain_shell,
        material_state,
        material_state_stride,
        use_green_lagrange,
        resultants_shell,
        tangent_shell
    );
}

"""
text = text[:end] + integrate_fn + text[end:]
text = replace_once(text, "    Precision*                    material_state,", "    const Precision*              material_state,", "ABD output definition const")
write(path, text)


# -----------------------------------------------------------------------------
# Integrated shell section: evaluate is frozen-history; integrate advances all MPs.
# -----------------------------------------------------------------------------
path = "src/section/section_shell_integrated.h"
text = read(path)
text = replace_once(text, "        Precision*                    material_state,", "        const Precision*              material_state,", "integrated shell evaluate const")
needle = """    ) const override;

    // Recover physical Cauchy stress at arbitrary z"""
insert = """    ) const override;

    // Advance all five through-thickness material points from accepted source
    // history into the corresponding target block.
    void integrate(
        const Vec3&                   position_reference,
        const Mat3&                   shell_basis_global,
        const ShellGeneralizedStrain& strain_shell,
        const Precision*              material_state,
        Precision*                    material_state_target,
        Index                         material_state_stride,
        bool                          use_green_lagrange,
        ShellStressResultants&        resultants_shell,
        Mat8&                         tangent_shell
    ) const override;

    // Recover physical Cauchy stress at arbitrary z"""
text = replace_once(text, needle, insert, "integrated shell integrate declaration")
text = replace_once(text, "        Precision*                    material_state,", "        const Precision*              material_state,", "integrated shell output const")
write(path, text)

path = "src/section/section_shell_integrated.cpp"
text = read(path)
start, end, eval_fn = block(text, "void IntegratedShellSection::evaluate(", "/**\n * Evaluates physical Cauchy stress")

eval_fn = replace_once(eval_fn, "    Precision*                    material_state,", "    const Precision*              material_state,", "integrated evaluate const signature")
validation = """    logging::error(material_->elasticity()->supports_shell_integration_linearized() || use_green_lagrange,
        "IntegratedShellSection material does not support linearized shell evaluation");
"""
validation_new = validation + """

    const auto elasticity = material_->elasticity();
    const bool stateful = elasticity->state_size() > 0;
    if (stateful) {
        logging::error(material_state != nullptr,
            "IntegratedShellSection: stateful material requires source history");
        logging::error(material_state_stride >= elasticity->state_size(),
            "IntegratedShellSection: material-state stride is smaller than the constitutive state");
    }
"""
eval_fn = replace_once(eval_fn, validation, validation_new, "integrated evaluate state validation")
eval_fn = replace_once(
    eval_fn,
    "        Precision*          state = material_state + mp * material_state_stride;",
    "        const Precision*    state = stateful\n"
    "            ? material_state + mp * material_state_stride\n"
    "            : nullptr;",
    "integrated evaluate state row"
)
eval_fn = replace_count(eval_fn, "            material_->elasticity()->evaluate(", "            elasticity->evaluate(", 2, "integrated evaluate dispatch")
text = text[:start] + eval_fn + text[end:]

# Clone the now-read-only generalized section integration and turn the clone into
# source->target constitutive integration.
start, end, eval_fn = block(text, "void IntegratedShellSection::evaluate(", "/**\n * Evaluates physical Cauchy stress")
integrate_fn = eval_fn
integrate_fn = replace_once(integrate_fn, "void IntegratedShellSection::evaluate(", "void IntegratedShellSection::integrate(", "integrated clone name")
integrate_fn = replace_once(
    integrate_fn,
    "    const Precision*              material_state,\n    Index                         material_state_stride,",
    "    const Precision*              material_state,\n    Precision*                    material_state_target,\n    Index                         material_state_stride,",
    "integrated clone target parameter"
)
integrate_fn = replace_once(
    integrate_fn,
    "        logging::error(material_state_stride >= elasticity->state_size(),\n            \"IntegratedShellSection: material-state stride is smaller than the constitutive state\");",
    "        logging::error(material_state_target != nullptr,\n"
    "            \"IntegratedShellSection: stateful material requires target history\");\n"
    "        logging::error(material_state != material_state_target,\n"
    "            \"IntegratedShellSection: source and target history must not alias\");\n"
    "        logging::error(material_state_stride >= elasticity->state_size(),\n"
    "            \"IntegratedShellSection: material-state stride is smaller than the constitutive state\");",
    "integrated clone target validation"
)
integrate_fn = replace_once(
    integrate_fn,
    "        const Precision*    state = stateful\n            ? material_state + mp * material_state_stride\n            : nullptr;",
    "        const Precision* state = stateful\n"
    "            ? material_state + mp * material_state_stride\n"
    "            : nullptr;\n"
    "        Precision* target_state = stateful\n"
    "            ? material_state_target + mp * material_state_stride\n"
    "            : nullptr;",
    "integrated clone target row"
)
integrate_fn = replace_count(
    integrate_fn,
    "            elasticity->evaluate(\n                material_strain_gl,\n                state,",
    "            elasticity->integrate(\n                material_strain_gl,\n                state,\n                target_state,",
    1,
    "integrated finite integrate dispatch"
)
integrate_fn = replace_count(
    integrate_fn,
    "            elasticity->evaluate(\n                material_strain_linearized,\n                state,",
    "            elasticity->integrate(\n                material_strain_linearized,\n                state,\n                target_state,",
    1,
    "integrated linear integrate dispatch"
)
integrate_doc = """
/**
 * Integrates all five through-thickness constitutive points from source to target
 * history while assembling the matching generalized resultants and tangent.
 */
"""
# Drop the evaluate function's long leading documentation from the clone by
# starting directly at the function definition; prepend a compact dedicated doc.
integrate_fn = integrate_doc + integrate_fn
text = text[:end] + integrate_fn + text[end:]

# Physical stress recovery remains strictly read-only and null-safe for stateless materials.
text = replace_once(text, "    Precision*                    material_state,", "    const Precision*              material_state,", "integrated output const signature")
output_validation = """    logging::error(material_->elasticity()->supports_shell_integration_linearized() || use_green_lagrange,
        "IntegratedShellSection material does not support linearized shell evaluation");
"""
output_validation_new = output_validation + """

    const auto elasticity = material_->elasticity();
    const bool stateful = elasticity->state_size() > 0;
    if (stateful) {
        logging::error(material_state != nullptr,
            "IntegratedShellSection: stateful material requires source history for stress recovery");
        logging::error(material_state_stride >= elasticity->state_size(),
            "IntegratedShellSection: material-state stride is smaller than the constitutive state");
    }
"""
# The validation text appears once in evaluate and once in output after the prior
# evaluate transformation, so only replace the last occurrence in output block.
out_start = text.index("VolumeStressCauchy IntegratedShellSection::evaluate_output_stress(")
out = text[out_start:]
out = replace_once(out, output_validation, output_validation_new, "integrated output state validation")
out = replace_once(
    out,
    "    Precision* state = material_state + state_mp * material_state_stride;",
    "    const Precision* state = stateful\n"
    "        ? material_state + state_mp * material_state_stride\n"
    "        : nullptr;",
    "integrated output state row"
)
out = replace_count(out, "        material_->elasticity()->evaluate(", "        elasticity->evaluate(", 2, "integrated output dispatch")
text = text[:out_start] + out
write(path, text)


# -----------------------------------------------------------------------------
# FRT shells: central null-safe history addressing and explicit integration mode.
# -----------------------------------------------------------------------------
path = "src/model/shell/frt_shell.h"
text = read(path)
needle = """    void compute_material_resultants(EvaluationData& data) const;

    // Section and material stiffness helpers."""
replacement = """    void compute_material_resultants(EvaluationData& data, bool integrate_history) const;

    // Accepted and target material-history addressing. Stateless sections return
    // null pointers and a zero stride, so elastic shell paths never dereference a
    // non-existent history field.
    Index            material_state_stride() const;
    const Precision* material_state_source(Index ip) const;
    Precision*       material_state_target(Index ip) const;

    // Section and material stiffness helpers."""
text = replace_once(text, needle, replacement, "FRT history helper declarations")
old = """        const Field*        ip_stress    = nullptr,
        int                 ip_start_idx = 0
    ) const;"""
new = """        const Field*        ip_stress    = nullptr,
        int                 ip_start_idx = 0,
        bool                integrate_material_history = false
    ) const;"""
text = replace_once(text, old, new, "FRT init evaluation integration flag")
write(path, text)

path = "src/model/shell/frt_shell_reference.inl"
text = read(path)
insert_at = text.index("/**\n * Evaluates the zero-strain generalized shell-section tangent.")
helpers = """/**
 * Returns the scalar row stride of material history used by this shell section.
 */
template<Index N>
Index FRTShell<N>::material_state_stride() const {
    ShellSection* section = shell_section();
    logging::error(section->material_ && section->material_->has_elasticity(),
                   "FRTShell: section material elasticity is not available for element ", this->elem_id);

    const auto elasticity = section->material_->elasticity();
    if (elasticity->state_size() == 0) return 0;

    logging::error(this->_model_data && this->_model_data->material_state,
                   "FRTShell: stateful material requires accepted source history for element ", this->elem_id);
    logging::error(this->_model_data->material_state->components >= elasticity->state_size(),
                   "FRTShell: source material-state field is too narrow for element ", this->elem_id);
    return this->_model_data->material_state->components;
}

/**
 * Returns the accepted source row of the first through-thickness material point.
 */
template<Index N>
const Precision* FRTShell<N>::material_state_source(Index ip) const {
    const Index stride = material_state_stride();
    if (stride == 0) return nullptr;

    const Index row = this->mp_index(ip, 0);
    logging::error(row < this->_model_data->material_state->rows,
                   "FRTShell: source material-state row is out of range for element ", this->elem_id);
    return &(*this->_model_data->material_state)(row, 0);
}

/**
 * Returns the working target row of the first through-thickness material point.
 */
template<Index N>
Precision* FRTShell<N>::material_state_target(Index ip) const {
    ShellSection* section = shell_section();
    const auto elasticity = section->material_->elasticity();
    if (elasticity->state_size() == 0) return nullptr;

    logging::error(this->_model_data && this->_model_data->material_state_target,
                   "FRTShell: stateful material requires target history for element ", this->elem_id);
    logging::error(this->_model_data->material_state_target->components >= elasticity->state_size(),
                   "FRTShell: target material-state field is too narrow for element ", this->elem_id);

    const Index row = this->mp_index(ip, 0);
    logging::error(row < this->_model_data->material_state_target->rows,
                   "FRTShell: target material-state row is out of range for element ", this->elem_id);
    return &(*this->_model_data->material_state_target)(row, 0);
}

"""
text = text[:insert_at] + helpers + text[insert_at:]
# resultant_stiffness is a frozen-history query.
text = replace_once(
    text,
    "    Precision* state = &(*this->_model_data->material_state)(this->mp_index(state_ip, 0), 0);",
    "    const Precision* state = material_state_source(state_ip);",
    "FRT resultant stiffness source state"
)
text = replace_once(
    text,
    "        this->_model_data->material_state->components,",
    "        material_state_stride(),",
    "FRT resultant stiffness stride"
)
write(path, text)

path = "src/model/shell/frt_shell_assembly.inl"
text = read(path)
text = replace_once(
    text,
    "void FRTShell<N>::compute_material_resultants(EvaluationData& data) const {",
    "void FRTShell<N>::compute_material_resultants(EvaluationData& data, bool integrate_history) const {",
    "FRT material resultants signature"
)
text = replace_once(
    text,
    "    const Index     state_stride = this->_model_data->material_state->components;",
    "    const Index     state_stride = material_state_stride();",
    "FRT material resultants stride"
)
text = replace_once(
    text,
    "        Precision*             state = &(*this->_model_data->material_state)(this->mp_index(ip, 0), 0);",
    "        const Precision* state = material_state_source(ip);\n"
    "        Precision* target_state = integrate_history\n"
    "            ? material_state_target(ip)\n"
    "            : nullptr;",
    "FRT material resultants rows"
)
old_call = """        section->evaluate(
            reference_position(point.r, point.s),
            basis,
            strain,
            state,
            state_stride,
            true,
            resultants,
            tangent
        );"""
new_call = """        if (integrate_history) {
            section->integrate(
                reference_position(point.r, point.s),
                basis,
                strain,
                state,
                target_state,
                state_stride,
                true,
                resultants,
                tangent
            );
        } else {
            section->evaluate(
                reference_position(point.r, point.s),
                basis,
                strain,
                state,
                state_stride,
                true,
                resultants,
                tangent
            );
        }"""
text = replace_once(text, old_call, new_call, "FRT evaluate/integrate dispatch")
# Nonlinear tangent must write a target history; ordinary stiffness remains frozen.
old_tangent_init = """    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        true,
        true
    );"""
new_tangent_init = """    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        true,
        true,
        nullptr,
        0,
        true
    );"""
# This exact block occurs once in stiffness() and once in stiffness_tangent().
# Change only the second occurrence by locating stiffness_tangent.
idx = text.index("MapMatrix FRTShell<N>::stiffness_tangent")
head, tail = text[:idx], text[idx:]
tail = replace_once(tail, old_tangent_init, new_tangent_init, "FRT nonlinear tangent integration flag")
text = head + tail
write(path, text)

path = "src/model/shell/frt_shell_kinematics.inl"
text = read(path)
text = replace_once(
    text,
    "    const Field*        ip_stress,\n    int                 ip_start_idx\n) const {",
    "    const Field*        ip_stress,\n    int                 ip_start_idx,\n    bool                integrate_material_history\n) const {",
    "FRT init evaluation definition flag"
)
# Drill stiffness is frozen-history and must also support stateless sections.
text = replace_once(
    text,
    "            const Index state_stride = this->_model_data->material_state->components;",
    "            const Index state_stride = material_state_stride();",
    "FRT drill state stride"
)
text = replace_once(
    text,
    "                Precision*             material_state =\n                    &(*this->_model_data->material_state)(this->mp_index(ip, 0), 0);",
    "                const Precision* material_state = material_state_source(ip);",
    "FRT drill source state"
)
text = replace_once(
    text,
    "            compute_material_resultants(data);",
    "            compute_material_resultants(data, integrate_material_history);",
    "FRT material resultant integration flag forwarding"
)
write(path, text)

path = "src/model/shell/frt_shell_output.inl"
text = read(path)
# Read-only arbitrary resultant recovery.
text = replace_once(
    text,
    "    Precision* state = &(*this->_model_data->material_state)(this->mp_index(state_ip, 0), 0);",
    "    const Precision* state = material_state_source(state_ip);",
    "FRT output resultant source state"
)
text = replace_once(text, "        this->_model_data->material_state->components,", "        material_state_stride(),", "FRT output resultant stride")
# Read-only physical stress recovery.
text = replace_once(
    text,
    "    Precision* state = &(*this->_model_data->material_state)(this->mp_index(state_ip, 0), 0);",
    "    const Precision* state = material_state_source(state_ip);",
    "FRT physical output source state"
)
text = replace_once(text, "        this->_model_data->material_state->components,", "        material_state_stride(),", "FRT physical output stride")
# Residual-only nonlinear stress state must advance target history.
marker = "void FRTShell<N>::compute_stress_state"
idx = text.index(marker)
head, tail = text[:idx], text[idx:]
old_init = """        const EvaluationData data = init_evaluation(
            state,
            true,
            false,
            false,
            true
        );"""
new_init = """        const EvaluationData data = init_evaluation(
            state,
            true,
            false,
            false,
            true,
            nullptr,
            0,
            true
        );"""
tail = replace_once(tail, old_init, new_init, "FRT residual integration flag")
text = head + tail
# Nodal section-force output remains read-only and null-safe.
text = replace_once(
    text,
    "        Precision* material_state =\n            &(*this->_model_data->material_state)(this->mp_index(state_ip, 0), 0);",
    "        const Precision* material_state = material_state_source(state_ip);",
    "FRT section force source state"
)
text = replace_once(text, "            this->_model_data->material_state->components,", "            material_state_stride(),", "FRT section force stride")
if "_model_data->material_state)(" in text:
    raise RuntimeError("FRT output still contains direct source-state row access")
write(path, text)


# -----------------------------------------------------------------------------
# J2 tests: verify source immutability, target update and read-only recovery.
# -----------------------------------------------------------------------------
path = "tests/test_material_j2.cpp"
text = read(path)
text = replace_once(text, "#include \"../src/material/isotropic_j2_elasticity.h\"", "#include \"../src/material/isotropic_j2_elasticity.h\"\n#include \"../src/material/strain/axial_strain_linearized.h\"\n#include \"../src/material/stress/axial_stress_cauchy.h\"", "J2 test strain includes")
append = r'''

TEST(Material_J2, SourceTargetIntegrationAndReadOnlyRecovery) {
    material::IsotropicJ2Elasticity j2(Precision(210000), Precision(0.3));
    j2.add_yield_point(Precision(250), Precision(0));
    j2.add_yield_point(Precision(275), Precision(0.01));
    j2.add_yield_point(Precision(300), Precision(0.02));

    Precision source[7];
    Precision target[7];
    j2.initialize_state(source);
    j2.initialize_state(target);

    const Precision source_before[7] {
        source[0], source[1], source[2], source[3], source[4], source[5], source[6]
    };

    AxialStressCauchy stress;
    Precision tangent = Precision(0);
    j2.integrate(
        AxialStrainLinearized(Precision(0.005)),
        source,
        target,
        stress,
        tangent
    );

    for (Index i = 0; i < 7; ++i) {
        EXPECT_EQ(source[i], source_before[i]);
    }
    EXPECT_GT(target[6], Precision(0));
    EXPECT_TRUE(std::isfinite(stress.value()));
    EXPECT_TRUE(std::isfinite(tangent));

    const Precision target_before[7] {
        target[0], target[1], target[2], target[3], target[4], target[5], target[6]
    };

    AxialStressCauchy recovered_stress_1;
    AxialStressCauchy recovered_stress_2;
    Precision recovered_tangent_1 = Precision(0);
    Precision recovered_tangent_2 = Precision(0);

    j2.evaluate(
        AxialStrainLinearized(Precision(0.005)),
        target,
        recovered_stress_1,
        recovered_tangent_1
    );
    j2.evaluate(
        AxialStrainLinearized(Precision(0.005)),
        target,
        recovered_stress_2,
        recovered_tangent_2
    );

    for (Index i = 0; i < 7; ++i) {
        EXPECT_EQ(target[i], target_before[i]);
    }
    EXPECT_NEAR(recovered_stress_1.value(), recovered_stress_2.value(), Precision(1e-10));
    EXPECT_NEAR(recovered_tangent_1, recovered_tangent_2, Precision(1e-8));
}
'''
if "SourceTargetIntegrationAndReadOnlyRecovery" in text:
    raise RuntimeError("J2 source-target test already exists")
text += append
write(path, text)


# -----------------------------------------------------------------------------
# Sanity checks for the completed architectural split.
# -----------------------------------------------------------------------------
checks = {
    "src/material/isotropic_j2_elasticity.cpp": [
        "IsotropicJ2Elasticity::integrate(",
        "store_state(target_state, update.state);",
    ],
    "src/material/isotropic_j2_state.cpp": [
        "IsotropicJ2Elasticity::evaluate(",
    ],
    "src/model/solid/element_solid_compute.ipp": [
        "integrate_material(point.r, point.s, point.t",
    ],
    "src/model/truss/truss.cpp": [
        "target_state(*this, elasticity)",
    ],
    "src/model/shell/frt_shell_assembly.inl": [
        "compute_material_resultants(EvaluationData& data, bool integrate_history)",
        "section->integrate(",
    ],
}
for file, needles in checks.items():
    current = read(file)
    for needle in needles:
        if needle not in current:
            raise RuntimeError(f"{file}: missing expected text {needle!r}")

print("Source/target constitutive refactor applied successfully")
