/**
 * @file register_functions.h
 * @brief Declares all native and Abaqus parser command registration functions.
 *
 * The command implementations are split into one translation unit per input
 * command. This header is the single declaration point used by the native and
 * Abaqus parsers as well as by those implementation files.
 *
 * Registration only constructs the DSL grammar. The parser owns semantic
 * execution order, while the registered callbacks translate parsed command
 * data into model and load-case state.
 *
 * @see Parser
 * @see ParserAbq
 * @see io::dsl::Registry
 *
 * @author Finn Eggers
 * @date 29.08.2026
 */

#pragma once

namespace fem::io::dsl {
struct Registry;
}

namespace fem::model {
struct Model;
}

namespace fem::io::reader {
class Parser;
class ParserAbq;
}

namespace fem::io::reader::commands {

// Native model and topology commands
void register_amplitude        (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_assembly         (fem::io::dsl::Registry& registry);
void register_beam_section     (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_cload            (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_connector        (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_contact          (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_convection       (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_conductivity     (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_coupling         (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_density          (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_dload            (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_elastic          (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_element          (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_elset            (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_end_assembly     (fem::io::dsl::Registry& registry);
void register_end_instance     (fem::io::dsl::Registry& registry);
void register_end_part         (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_equation         (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_field            (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_heading          (fem::io::dsl::Registry& registry);
void register_heat_flux        (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_hyperelastic     (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_inertialload     (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_instance         (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_mass             (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_material         (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_node             (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_normal           (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_nset             (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_orientation      (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_overview         (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_part             (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_pload            (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_point_mass       (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_profile          (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_rbm              (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_rotary_inertia   (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_sfset            (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_shell_section    (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_solid_section    (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_spring           (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_support          (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_surface          (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_temperature      (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_thermal_expansion(fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_tie              (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_tload            (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_truss_section    (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_vload            (fem::io::dsl::Registry& registry, fem::model::Model& model);

// Native load-case commands operating on parser-owned analysis state
void register_loadcase_begin            (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_constraintmethod (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_constraintsummary(fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_damping          (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_frequency        (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_inertiarelief    (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_initialvelocity  (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_loads            (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_newmark          (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_nonlinear        (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_numeigenvalues   (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_rebalance        (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_request_stgeom   (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_request_stiffness(fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_sigma            (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_solver           (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_supports         (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_time             (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_topodensity      (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_topoexponent     (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_topoorient       (fem::io::dsl::Registry& registry, Parser& parser);
void register_loadcase_write_every      (fem::io::dsl::Registry& registry, Parser& parser);

} // namespace fem::io::reader::commands

namespace fem::io::reader::commands_abq {

// Abaqus syntax translations backed directly by the FEMaster model
void register_amplitude    (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_element      (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_expansion    (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_orientation  (fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_shell_section(fem::io::dsl::Registry& registry, fem::model::Model& model);
void register_solid_section(fem::io::dsl::Registry& registry, fem::model::Model& model);

// Abaqus commands operating on parser-owned step and load state
void register_boundary (fem::io::dsl::Registry& registry, ParserAbq& parser);
void register_cload    (fem::io::dsl::Registry& registry, ParserAbq& parser);
void register_dload    (fem::io::dsl::Registry& registry, ParserAbq& parser);
void register_dsload   (fem::io::dsl::Registry& registry, ParserAbq& parser);
void register_step     (fem::io::dsl::Registry& registry, ParserAbq& parser);
void register_transform(fem::io::dsl::Registry& registry, ParserAbq& parser);

} // namespace fem::io::reader::commands_abq
