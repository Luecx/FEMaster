/**
 * @file model_build.cpp
 * @brief Implements model-level field construction and structural assembly.
 *
 * The routines operate on the dense assembly created by `Model::compile()`.
 * They bind sections to compiled elements, construct shared shell reference
 * normals, enumerate active generalized DOFs, assemble loads and constraint
 * equations, and combine regular and post-compile point elements into global
 * sparse operators and element-nodal internal-force fields. Contact forces and
 * tangents are assembled separately.
 *
 * FEMaster `POINTMASS` commands create auxiliary `PointElement` objects after
 * compilation because they target compiled NSETs. Those point elements are kept
 * outside dense ELSET enumeration, but own element-nodal rows behind the dense
 * rows and participate in the same DOF, stiffness, force and mass operators as
 * regular structural elements.
 *
 * Element-local matrix evaluation remains the responsibility of each structural
 * formulation. `mattools::assemble_matrix()` owns local-to-global DOF mapping,
 * sparsity construction and atomic matrix accumulation. Internal forces are
 * written by every element into its own ELEMENT_NODAL rows and reduced onto
 * global nodes by the caller.
 *
 * @see Model
 * @see Model::compile
 * @see mattools::assemble_matrix
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "../core/parallel.h"
#include "../mattools/assemble.h"
#include "element/element_structural.h"
#include "model.h"
#include "solid/element_solid.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

namespace fem {
namespace model {

/**
 * Binds every compiled section assignment to its target elements.
 *
 * Section regions already contain dense assembly element identifiers. Each
 * represented non-null element receives the shared section pointer so later
 * stiffness, mass, constitutive and result-recovery operations can access the
 * assigned material and section properties. Auxiliary point elements are bound
 * directly when the post-compile POINTMASS command creates them.
 */
void Model::assign_sections() {
    // Apply each compiled section to all non-null elements in its dense region
    for (Section::Ptr ptr : _data->sections) {
        for (ID elem_id : *ptr->region_) {
            if (_data->elements[elem_id] != nullptr) {
                _data->elements[elem_id]->set_section(ptr);
            }
        }
    }
}

/**
 * Constructs consistent reference normals at every compiled shell element node.
 *
 * Existing finite non-zero entries in `SHELL_ELEMENT_NODAL_NORMALS` are treated
 * as prescribed directions and normalized without being overwritten. Missing
 * entries are evaluated from the shell reference surface at its natural nodal
 * coordinates. Contributions sharing a global node are clustered by angular
 * similarity and each cluster is averaged independently, preserving geometric
 * creases whose normal angle exceeds `equalize_angle_degrees`.
 *
 * The resulting three-component `ELEMENT_NODAL` field remains aligned with the
 * compiled element-nodal offsets. Prescribed rows participate in cluster
 * directions but retain their individual normalized values.
 *
 * @param equalize_angle_degrees Maximum angle in degrees between normals that
 *                               may be averaged into one nodal cluster.
 */
void Model::build_shell_element_normals(Precision equalize_angle_degrees) {
    // Validate the reference geometry and compiled element-nodal enumeration
    logging::error(_data->positions_reference != nullptr,
        "Model: POSITION_REFERENCE field is not set");
    logging::error(_data->element_nodal_offsets != nullptr,
        "Model: element nodal offsets are not initialized");

    // Convert the angular clustering threshold and establish numerical tolerances
    constexpr Precision normal_tolerance = Precision(1e-12);

    const Precision pi             = std::acos(Precision(-1));
    const Precision equalize_angle = equalize_angle_degrees * pi / Precision(180);
    const Precision cos_equalize   = std::cos(equalize_angle);
    const Index element_count      = static_cast<Index>(_data->elements.size());
    const Index expected_rows      = _data->field_rows(FieldDomain::ELEMENT_NODAL);

    // Models without element-local rows do not require shell-normal storage
    if (expected_rows == 0) {
        _data->shell_element_nodal_normals = nullptr;
        return;
    }

    // Create the result field or validate a pre-existing field with prescribed data
    Field::Ptr normals = _data->shell_element_nodal_normals;
    if (normals == nullptr) {
        normals = _data->create_field(
            "SHELL_ELEMENT_NODAL_NORMALS", FieldDomain::ELEMENT_NODAL, 3, false);
        normals->set_zero();
    } else {
        logging::error(normals->domain == FieldDomain::ELEMENT_NODAL,
            "Model: shell normal field '", normals->name,
            "' must use ELEMENT_NODAL domain");
        logging::error(normals->components == 3,
            "Model: shell normal field '", normals->name,
            "' must have exactly three components");
        logging::error(normals->rows == expected_rows,
            "Model: shell normal field '", normals->name,
            "' has ", normals->rows, " rows, expected ", expected_rows);
    }

    // Prepare global-node adjacency and one normal value per element-nodal row
    const Index node_count = _data->positions_reference->rows;
    std::vector<std::vector<Index>> node_rows(static_cast<std::size_t>(node_count));
    std::vector<Vec3>               row_normals(static_cast<std::size_t>(normals->rows), Vec3::Zero());
    std::vector<bool>               prescribed_rows(static_cast<std::size_t>(normals->rows), false);

    // Evaluate or normalize the reference normal at every shell element node
    for (const ElementPtr& element : _data->elements) {
        if (element == nullptr) continue;

        const auto structural = element->as<StructuralElement>();
        if (structural == nullptr || !structural->is_shell()) continue;

        const auto surface = structural->surface(1);
        logging::error(surface != nullptr,
            "Model: shell element ", structural->elem_id,
            " does not provide a reference surface");

        const Index offset  = static_cast<Index>(structural->elem_nodal_offset);
        const Index n_nodes = static_cast<Index>(structural->n_nodes());

        // Natural nodal coordinates must follow the shell connectivity ordering
        const DynamicMatrix natural_coords = surface->node_coords_natural();
        logging::error(static_cast<Index>(natural_coords.rows()) == n_nodes && static_cast<Index>(natural_coords.cols()) == 2,
            "Model: natural surface node coordinates do not match shell element ",
            structural->elem_id);

        for (Index local_node = 0; local_node < n_nodes; ++local_node) {
            const ID    node_id = structural->nodes()[local_node];
            const Index row     = offset + local_node;

            Vec3 normal = normals->row_vec3(row);
            const bool prescribed = normal.allFinite() && normal.norm() > normal_tolerance;

            // Preserve supplied directions and reconstruct only missing rows
            if (prescribed) {
                normal.normalize();
                prescribed_rows[static_cast<std::size_t>(row)] = true;
            } else {
                const Vec2 local = natural_coords.row(local_node).transpose();
                normal = surface->normal(*_data->positions_reference, local);

                logging::error(normal.allFinite() && normal.norm() > normal_tolerance,
                    "Model: invalid shell normal in element ", structural->elem_id,
                    " at local node ", local_node);
                normal.normalize();
            }

            row_normals[static_cast<std::size_t>(row)] = normal;
            node_rows[static_cast<std::size_t>(node_id)].push_back(row);
        }
    }

    // Cluster adjacent element normals independently at every shared global node
    for (const auto& rows : node_rows) {
        std::vector<Vec3>               cluster_sums;
        std::vector<std::vector<Index>> cluster_rows;

        for (Index row : rows) {
            const Vec3 normal = row_normals[static_cast<std::size_t>(row)];
            bool added = false;

            for (Index cluster = 0; cluster < static_cast<Index>(cluster_sums.size()); ++cluster) {
                const Vec3 cluster_normal =
                    cluster_sums[static_cast<std::size_t>(cluster)].normalized();

                if (normal.dot(cluster_normal) < cos_equalize) continue;

                cluster_sums[static_cast<std::size_t>(cluster)] += normal;
                cluster_rows[static_cast<std::size_t>(cluster)].push_back(row);
                added = true;
                break;
            }

            // Start a separate cluster when the normal crosses a geometric crease
            if (!added) {
                cluster_sums.push_back(normal);
                cluster_rows.push_back({row});
            }
        }

        // Write each averaged cluster direction to its non-prescribed rows
        for (Index cluster = 0; cluster < static_cast<Index>(cluster_sums.size()); ++cluster) {
            const Vec3 equalized_normal =
                cluster_sums[static_cast<std::size_t>(cluster)].normalized();

            for (Index row : cluster_rows[static_cast<std::size_t>(cluster)]) {
                if (prescribed_rows[static_cast<std::size_t>(row)]) continue;
                (*normals)(row, 0) = equalized_normal(0);
                (*normals)(row, 1) = equalized_normal(1);
                (*normals)(row, 2) = equalized_normal(2);
            }
        }
    }

    // Publish the completed field through the dedicated solver-facing handle
    _data->shell_element_nodal_normals = normals;
}

/**
 * Enumerates every unconstrained generalized DOF required by the model.
 *
 * Dense elements and post-compile point elements activate the components used by
 * their formulations. Couplings and connectors may additionally activate
 * generalized nodal directions. Manual equations activate every component with
 * a nonzero coefficient, including components of nodes without elements. These
 * components must remain in the algebraic system until constraint reduction;
 * dropping them would change the equations into unintended supports.
 * The activated components are numbered contiguously from zero.
 *
 * @return Node-by-six numbering of the unconstrained system DOFs.
 */
SystemDofs Model::build_structural_dofs() {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");

    SystemDofs system_dofs{_data->positions->rows, 6};

    // Activate generalized DOFs used by regular compiled elements
    for (auto& element : _data->elements) {
        if (element == nullptr) continue;

        const auto dofs = element->dofs();
        for (ID local_node = 0; local_node < element->n_nodes(); ++local_node) {
            system_dofs.activate(element->nodes()[local_node], dofs);
        }
    }

    // Native POINTMASS objects are point elements created after compilation.
    for (auto& element : _data->point_elements) {
        if (element == nullptr) continue;

        const auto dofs = element->dofs();
        for (ID local_node = 0; local_node < element->n_nodes(); ++local_node) {
            system_dofs.activate(element->nodes()[local_node], dofs);
        }
    }

    // Preserve manual equation terms even when their nodes have no elements.
    // Only the referenced directions are activated. Do this before deriving
    // coupling and MPC master DOFs so their masks include these directions.
    for (const auto& equation : _data->equations) {
        for (const auto& entry : equation.entries) {
            logging::error(entry.node_id >= 0 && entry.node_id < system_dofs.nodes()
                        && entry.dof < 6 && std::isfinite(entry.coeff),
                "Model: invalid manual equation term (node=", entry.node_id,
                ", dof=", int(entry.dof), ", coefficient=", entry.coeff, ")");

            if (entry.coeff != Precision(0)) {
                system_dofs.activate(entry.node_id, entry.dof);
            }
        }
    }

    // Add master-node components required by coupling and rigid-body kinematics.
    // Both may use each other's reference nodes as slaves, so repeat until no
    // further reference-node component is activated.
    const auto activate_master = [&system_dofs](ID master_id, const Dofs& dofs) {
        bool changed = false;
        for (Index dof = 0; dof < system_dofs.width(); ++dof) {
            changed = changed || (dofs(dof) && !system_dofs.active(master_id, dof));
        }
        system_dofs.activate(master_id, dofs);
        return changed;
    };

    for (bool changed = true; changed;) {
        changed = false;
        for (auto& coupling : _data->couplings) {
            changed = activate_master(coupling.master_node, coupling.master_dofs(system_dofs, *_data)) || changed;
        }
        for (auto& rigid_body : _data->rigid_bodies) {
            changed = activate_master(rigid_body.ref_node, rigid_body.master_dofs(system_dofs, *_data)) || changed;
        }
    }

    // Activate constrained components at both nodes of every connector
    for (auto& connector : _data->connectors) {
        const ID node1_id = connector.node_1();
        const ID node2_id = connector.node_2();
        const auto dofs   = connector.dofs();

        system_dofs.activate(node1_id, dofs);
        system_dofs.activate(node2_id, dofs);
    }

    // Add independent-node components required by two-node MPC kinematics
    for (const auto& mpc : _data->mpcs) {
        system_dofs.activate(mpc.node_2(), mpc.master_dofs(system_dofs));
    }

    system_dofs.numerate();
    return system_dofs;
}

/**
 * Assembles the current structural condition state into global nodal forces.
 *
 * Active structural conditions are traversed by family. The parser has already
 * resolved input history and selected named collectors for this analysis.
 * Conditions shared by both states contribute their nominal force exactly
 * once. Removed conditions decay linearly, while added conditions follow their
 * named step-time amplitude (or a linear ramp when none is assigned). Dynamic
 * callers provide an empty initial state and evaluate target amplitudes at time.
 *
 * Structural loads write the six-component load field. If the caller requests
 * a load tangent, every load appends its analytic derivative df/du with the
 * same multiplier as its load contribution; dead loads contribute nothing.
 * Follower loads require the tangent and report an error without it.
 *
 * @param initial Frozen initial condition state.
 * @param target Prescribed target condition state.
 * @param time Evaluation time supplied to load amplitudes.
 * @param step_progress Normalized progress for a static load transition.
 * @param system_dof_ids Global equation numbering of the load tangent.
 * @param load_tangent Optional triplets receiving df/du; null requests none.
 * @return Six-component global nodal load field for the current model state.
 */
Field Model::build_load_matrix(const bc::ConditionState& initial,
                               const bc::ConditionState& target,
                               Precision                 time,
                               Precision                 step_progress,
                               const SystemDofs&         system_dof_ids,
                               TripletList*              load_tangent) {
    // Initialize global generalized forces and the outputs unused by loads
    Field load = _data->create_field_("LOAD_MATRIX", FieldDomain::NODE, 6, false);
    load.set_zero();

    constraint::Equations constraints{};
    TripletList           system_matrix{};

    // Partition each family by pointer identity without creating new load
    // objects. A shared pointer represents the same physical definition in
    // both states, so its initial ramp-down and target ramp-up cancel:
    //
    //     (1 - lambda) * f + lambda * f = f.
    //
    // Assemble this intersection only once at its constant nominal magnitude.
    // The parser ensures newly activated static amplitudes start at zero and
    // finish at one so the next step can inherit the nominal load unchanged.
    for (const auto family : bc::STRUCTURAL_LOAD_FAMILIES) {
        // Initial-only definitions are removed linearly. Definitions also
        // present in the target remain constant, regardless of any amplitude
        // attached to their original physical definition.
        for (const auto& condition : initial.get(family)) {
            const bool      continued = target.contains(family, condition);
            const Precision scale     = continued
                ? Precision(1) : Precision(1) - step_progress;

            condition->apply(*_data, load, load_tangent, constraints, system_dof_ids, system_matrix, scale);
        }

        // Target-only definitions are newly activated. Their own amplitude
        // supplies the step-time multiplier; without one, use the standard
        // linear ramp. The intersection has already been assembled above.
        for (const auto& condition : target.get(family)) {
            if (initial.contains(family, condition)) continue;

            const Precision scale = condition->amplitude_
                ? condition->amplitude_->evaluate(time) : step_progress;

            condition->apply(*_data, load, load_tangent, constraints, system_dof_ids, system_matrix, scale);
        }
    }

    return load;
}

/**
 * Decomposes current structural loads into unscaled spatial amplitude bases.
 *
 * Current structural families supply precisely the same active definitions
 * as build_load_matrix(). Loads sharing an amplitude pointer
 * contribute to one six-component nodal field; conditions without amplitude
 * share the null entry. Evaluation through Condition::apply bypasses amplitude
 * scaling, so the transient procedure can multiply each spatial basis by its
 * scalar time history without repeating geometric integration.
 *
 * These fields are derived numerical data for the running analysis. They retain
 * no condition definitions or condition activations.
 *
 * @return Amplitude pointers paired with unscaled global nodal load fields.
 */
std::vector<std::pair<bc::Amplitude::Ptr, Field>> Model::build_load_basis(const bc::ConditionState& target) {
    // Prepare the numerical spatial bases and outputs unused by structural loads
    std::vector<std::pair<bc::Amplitude::Ptr, Field>> basis;
    constraint::Equations constraints{};
    SystemDofs            system_dof_ids{};
    TripletList           system_matrix{};

    // Transient load bases are assembled once and reused over time. A follower
    // load changes with displacement, so it must be rejected before any load
    // integration (also for zero scaling or an empty surface region).
    for (const auto family : bc::STRUCTURAL_LOAD_FAMILIES) {
        for (const auto& condition : target.get(family)) {
            logging::error(!condition->follower_,
                "LOADS: FOLLOWER requires a nonlinear load tangent; transient load bases are unsupported");
        }
    }

    // Reuse one field for each shared scalar amplitude, including the null history
    auto add_condition = [&](const bc::Condition::Ptr& condition) {
        auto entry = std::find_if(basis.begin(), basis.end(), [&](const auto& candidate) {
            return candidate.first == condition->amplitude_;
        });
        if (entry == basis.end()) {
            Field load_matrix = _data->create_field_("LOAD_BASIS", FieldDomain::NODE, 6, false);
            load_matrix.set_zero();
            basis.emplace_back(condition->amplitude_, std::move(load_matrix));
            entry = std::prev(basis.end());
        }

        // Retain nominal magnitudes; the transient solver evaluates time scaling.
        // Linear dynamics requests no load tangent, so follower loads report
        // an error here.
        condition->apply(*_data, entry->second, nullptr, constraints, system_dof_ids, system_matrix, Precision(1));
    };

    // Use the same independent structural input families as direct RHS assembly

    for (const auto family : bc::STRUCTURAL_LOAD_FAMILIES) {
        for (const auto& condition : target.get(family)) {
            add_condition(condition);
        }
    }

    return basis;
}

/**
 * Generates and categorizes the current model's linear constraint equations.
 *
 * Target SUPPORT conditions contribute prescribed displacements and rotations.
 * Inactive reusable collector definitions do not impose constraints. Each prescription is
 * expanded through the common Condition::apply interface and annotated with its
 * origin for rank diagnostics and constraint reporting.
 *
 * Connector, MPC, coupling, tie and rigid-body-mode objects then generate their
 * equations in model order. Manually stored equations enter the fallback group.
 * Algebraic reduction and enforcement remain with the constraint transformer.
 *
 * @param system_dof_ids Active global DOF identifiers used by constraint builders.
 * @return Constraint equations grouped and annotated by their origin.
 */
constraint::ConstraintGroups Model::collect_constraints(SystemDofs& system_dof_ids,
                                                           const bc::ConditionState& target) {
    // Prepare grouped equations and temporary data required by support prescriptions
    constraint::ConstraintGroups groups{};
    Field                       load{};
    TripletList                 system_matrix{};

    // Collect equations from a sequence of constraint objects. Each generated
    // equation retains the source type and index of its originating object.
    auto collect = [&](auto& items, auto source, auto& target, auto get_equations) {
        Index index = 0;
        for (auto& item : items) {
            // Generate all equations associated with the current constraint object
            auto equations = get_equations(item);

            // Assign the same source index to all equations from this object
            // and move them into their corresponding constraint group
            for (auto& equation : equations) {
                equation.source       = source;
                equation.source_index = index;
                target.push_back(std::move(equation));
            }

            // Advance once per constraint object, not per generated equation
            ++index;
        }
    };

    // Generate equations from topology-dependent constraints using the active
    // system DOF numbering and the current model data
    auto get_equations = [&](auto& item) {
        return item.get_equations(system_dof_ids, *_data);
    };

    // Expand prescribed supports into equations using the condition interface.
    // The temporary load and system matrix are required by apply(), while
    // the generated equations are returned for grouping.
    auto get_supports = [&](const bc::Condition::Ptr& condition) {
        constraint::Equations equations{};
        condition->apply(*_data, load, nullptr, equations, system_dof_ids, system_matrix, Precision(1));
        return equations;
    };

    // Collect active supports and topology-dependent constraints in model order.
    // Each constraint family receives its own source indexing and output group.
    using Source = constraint::EquationSourceKind;

    collect(target.get(bc::SUPPORT), Source::Support,   groups.supports,   get_supports);
    collect(_data->connectors,                  Source::Connector, groups.connectors, get_equations);
    collect(_data->mpcs,                        Source::Mpc,       groups.mpcs,       get_equations);
    collect(_data->couplings,                   Source::Coupling,  groups.couplings,  get_equations);
    collect(_data->rigid_bodies,                Source::RigidBody, groups.rigid_bodies, get_equations);
    collect(_data->ties,                        Source::Tie,       groups.ties,       get_equations);
    collect(_data->rbms,                        Source::Rbm,       groups.rbms,       get_equations);

    // Append manually defined equations without regenerating their coefficients.
    // Preserve explicitly assigned source types and mark unknown origins as manual.
    Index manual_idx = 0;
    for (auto equation : _data->equations) {
        if (equation.source == Source::Unknown) {
            equation.source = Source::Manual;
        }

        // Preserve the original ordering for diagnostics and source identification
        equation.source_index = manual_idx++;
        groups.others.push_back(std::move(equation));
    }

    return groups;
}

/**
 * Assembles the linear global structural stiffness matrix.
 *
 * Regular compiled elements and post-compile native point elements both use the
 * common sparse element assembly path. The optional ELEMENT-domain stiffness
 * scale applies only to dense compiled elements because auxiliary point elements
 * have no element-domain row.
 *
 * @param indices Active global DOF identifiers used for local-to-global mapping.
 * @param stiffness_scalar Optional one-component element stiffness scale.
 * @return Assembled sparse linear stiffness matrix.
 */
SparseMatrix Model::build_stiffness_matrix(
    SystemDofs&   indices,
    const Field*  base_temperature,
    const Field*  stiffness_scalar
) {
    logging::error(_data->contacts.empty(),
        "CONTACT requires NONLINEARSTATIC; linear stiffness assembly cannot include contact");

    auto lambda = [&](const ElementPtr& element, Precision* storage) {
        if (auto structural = element->as<StructuralElement>()) {
            MapMatrix stiffness = structural->evaluate(
                storage, nullptr, nullptr,
                nullptr, nullptr,
                nullptr, base_temperature,
                false);
            if (stiffness_scalar) {
                logging::error(stiffness_scalar->domain == FieldDomain::ELEMENT,
                    "stiffness scale field must use ELEMENT domain");
                logging::error(stiffness_scalar->components == 1,
                    "stiffness scale field must have 1 component");
                stiffness *= (*stiffness_scalar)(static_cast<Index>(structural->elem_id), 0);
            }
            return stiffness;
        }

        MapMatrix matrix{storage, 0, 0};
        return matrix;
    };

    SparseMatrix matrix = mattools::assemble_matrix(_data->elements, indices, lambda);

    if (!_data->point_elements.empty()) {
        auto point_lambda = [&](const ElementPtr& element, Precision* storage) {
            if (auto structural = element->as<StructuralElement>()) {
                return structural->evaluate(
                    storage, nullptr, nullptr,
                    nullptr, nullptr,
                    nullptr, base_temperature,
                    false);
            }
            MapMatrix matrix{storage, 0, 0};
            return matrix;
        };
        matrix += mattools::assemble_matrix(_data->point_elements, indices, point_lambda);
    }

    return matrix;
}

/**
 * Assembles the nonlinear tangent and matching element-nodal internal force.
 *
 * Every structural element performs one physical trial evaluation through
 * `evaluate()`. The same constitutive response supplies its internal
 * force and, when matrix storage is present, its consistent tangent. No global
 * integration-point stress scratch field is allocated during assembly.
 *
 * Native POINTMASS elements use the same path: their linear spring law produces
 * the matching force and constant tangent directly.
 *
 * A supplied pattern is reused for regular elements and must match their
 * connectivity and active DOF numbering. Point elements are assembled
 * separately without modifying the supplied pattern. Contact is not included;
 * see build_contact_force().
 *
 * @param indices Active global DOF identifiers used for sparse assembly.
 * @param element_forces Six-component ELEMENT_NODAL field overwritten with the
 *                       element internal forces.
 * @param target_displacement Requested global nodal displacement field.
 * @param target_temperature Requested global nodal temperature field.
 * @param base_displacement Base displacement state used for linearization.
 * @param base_temperature Base temperature state used for linearization.
 * @param stiffness_scalar Optional one-component element tangent scale.
 * @param update_state Whether element evaluations may write trial constitutive history.
 * @param pattern Optional compressed sparsity pattern for regular elements.
 * @return Assembled sparse nonlinear tangent matrix.
 */
SparseMatrix Model::build_tangent_stiffness_matrix(
    SystemDofs&         indices,
    Field&              element_forces,
    const Field&        target_displacement,
    const Field*        target_temperature,
    const Field&        base_displacement,
    const Field*        base_temperature,
    const Field*        stiffness_scalar,
    bool                update_state,
    const SparseMatrix* pattern
) {
    logging::error(element_forces.domain == FieldDomain::ELEMENT_NODAL,
        "tangent internal force output must use ELEMENT_NODAL domain");
    logging::error(element_forces.rows == _data->field_rows(FieldDomain::ELEMENT_NODAL),
        "tangent internal force output has wrong row count");
    logging::error(element_forces.components >= 6,
        "tangent internal force output requires at least 6 components");
    logging::error(!stiffness_scalar || stiffness_scalar->domain == FieldDomain::ELEMENT,
        "stiffness scale field must use ELEMENT domain");
    logging::error(!stiffness_scalar || stiffness_scalar->components == 1,
        "stiffness scale field must have 1 component");

    element_forces.set_zero();

    // Element rows are disjoint, so every callback writes its forces directly
    auto lambda = [&](const ElementPtr& element, Precision* storage) -> MapMatrix {
        auto* structural = element->as<StructuralElement>();
        if (!structural) {
            MapMatrix matrix{storage, 0, 0};
            return matrix;
        }

        MapMatrix tangent = structural->evaluate(
            storage,
            nullptr,
            &element_forces,
            &target_displacement,
            target_temperature,
            &base_displacement,
            base_temperature,
            update_state
        );

        // Point elements have negative ids and are never scaled
        if (stiffness_scalar && structural->elem_id >= 0) {
            tangent *= (*stiffness_scalar)(static_cast<Index>(structural->elem_id), 0);
        }

        return tangent;
    };

    SparseMatrix global_matrix = mattools::assemble_matrix(_data->elements, indices, lambda, pattern);

    if (!_data->point_elements.empty()) {
        global_matrix += mattools::assemble_matrix(_data->point_elements, indices, lambda);
    }

    element_forces.check_finite("Internal force");
    return global_matrix;
}

/**
 * Evaluates element-nodal internal forces without retaining a tangent.
 *
 * Structural elements use the same physical `evaluate()` path as full
 * Newton assembly, but receive a null matrix buffer. They therefore update trial
 * material state and write their matching internal force without constructing
 * tangent-only terms. Every element owns disjoint ELEMENT_NODAL rows, so regular
 * and point elements are evaluated in one parallel loop without private fields.
 * Contact is not included; see build_contact_force().
 *
 * @param element_forces Six-component ELEMENT_NODAL field overwritten with the
 *                       element internal forces.
 * @param target_displacement Requested global nodal displacement field.
 * @param target_temperature Requested global nodal temperature field.
 * @param base_displacement Optional base displacement state; nullptr denotes zero.
 * @param base_temperature Optional base temperature state; nullptr denotes the
 *                         material stress-free temperature.
 * @param update_state If the material state shall be updated or not
 */
void Model::build_internal_force(
    Field&        element_forces,
    const Field&  target_displacement,
    const Field*  target_temperature,
    const Field*  base_displacement,
    const Field*  base_temperature,
    bool          update_state
) {
    logging::error(element_forces.domain == FieldDomain::ELEMENT_NODAL,
        "internal force output must use ELEMENT_NODAL domain");
    logging::error(element_forces.rows == _data->field_rows(FieldDomain::ELEMENT_NODAL),
        "internal force output has wrong row count");
    logging::error(element_forces.components >= 6,
        "internal force output requires at least 6 components");

    element_forces.set_zero();

    const std::size_t regular_count = _data->elements.size();
    const std::size_t element_count = regular_count + _data->point_elements.size();

    parallel::for_index(
        element_count,
        static_cast<int>(global_config.max_threads),
        [&](std::size_t element_index, int /*worker*/) {
            const ElementPtr& element =
                element_index < regular_count
                    ? _data->elements[element_index]
                    : _data->point_elements[element_index - regular_count];

            if (!element) return;

            auto* structural = element->as<StructuralElement>();
            if (!structural) return;

            structural->evaluate(
                nullptr, nullptr,
                &element_forces,
                &target_displacement, target_temperature,
                base_displacement, base_temperature,
                update_state
            );
        },
        std::size_t(64)
    );

    element_forces.check_finite("Internal force");
}

/**
 * Adds contact forces and, optionally, the matching contact tangent.
 *
 * Contact couples changing node pairs and is therefore not part of the
 * element-nodal layout. Its forces are added directly to global nodes, usually
 * onto the reduced element internal force from Field::nodal_sum().
 *
 * @param indices Active global DOF identifiers used by the contact tangent.
 * @param nodal_forces Six-component NODE field receiving contact forces.
 * @param tangent Optional active-numbering matrix receiving the contact tangent.
 */
void Model::build_contact_force(SystemDofs&   indices,
                                NodeData&     nodal_forces,
                                SparseMatrix* tangent) {
    logging::error(nodal_forces.domain == FieldDomain::NODE,
        "contact force output must use NODE domain");
    logging::error(nodal_forces.rows == _data->field_rows(FieldDomain::NODE),
        "contact force output has wrong node count");
    logging::error(nodal_forces.components >= 6,
        "contact force output requires at least 6 components");

    if (_data->contacts.empty()) return;

    TripletList triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, triplets);
    }

    if (tangent && !triplets.empty()) {
        SparseMatrix contact_matrix(tangent->rows(), tangent->cols());
        contact_matrix.insertFromTriplets(triplets.begin(), triplets.end());
        *tangent += contact_matrix;
    }

    nodal_forces.check_finite("Contact force");
}

/**
 * Assembles the global perturbation geometric stiffness between two states.
 *
 * Target displacement and temperature define (u,T), while base displacement
 * and temperature define (u0,T0). Each structural element forms the stress or
 * stress-resultant increment relative to the base state and assembles the
 * corresponding geometric stiffness locally. Stress already present at
 * (u0,T0) belongs to the complete tangent and is not repeated here.
 *
 * No global integration-point stress field is created or consumed. Auxiliary
 * point elements are excluded because they carry no stress-dependent geometric
 * stiffness.
 *
 * @param indices Active global DOF identifiers used for sparse assembly.
 * @param target_displacement Requested global nodal displacement state u.
 * @param target_temperature Requested global nodal temperature state T.
 * @param base_displacement Optional base displacement state u0.
 * @param base_temperature Optional base temperature state T0.
 * @param stiffness_scalar Optional one-component element stiffness scale.
 * @return Assembled sparse perturbation geometric stiffness matrix.
 */
SparseMatrix Model::build_geom_stiffness_matrix(
    SystemDofs&   indices,
    const Field&  target_displacement,
    const Field*  target_temperature,
    const Field*  base_displacement,
    const Field*  base_temperature,
    const Field*  stiffness_scalar
) {
    auto lambda = [&](const ElementPtr& element, Precision* storage) -> MapMatrix {
        if (auto structural = element->as<StructuralElement>()) {
            MapMatrix geometric_stiffness = structural->evaluate(
                nullptr,
                storage,
                nullptr,
                &target_displacement,
                target_temperature,
                base_displacement,
                base_temperature,
                false
            );

            if (stiffness_scalar) {
                logging::error(stiffness_scalar->domain == FieldDomain::ELEMENT,
                    "stiffness scale field must use ELEMENT domain");
                logging::error(stiffness_scalar->components == 1,
                    "stiffness scale field must have 1 component");
                geometric_stiffness *=
                    (*stiffness_scalar)(static_cast<Index>(structural->elem_id), 0);
            }

            return geometric_stiffness;
        }

        MapMatrix matrix{storage, 0, 0};
        return matrix;
    };

    return mattools::assemble_matrix(_data->elements, indices, lambda);
}

/**
 * Assembles the global lumped mass matrix.
 *
 * Regular structural elements and post-compile native point elements both
 * provide local translational/rotational mass matrices through the common sparse
 * assembly path.
 *
 * @param indices Active global DOF identifiers used for sparse assembly.
 * @return Assembled sparse lumped mass matrix.
 */
SparseMatrix Model::build_lumped_mass_matrix(SystemDofs& indices) {
    auto lambda = [&](const ElementPtr& element, Precision* storage) {
        if (auto structural = element->as<StructuralElement>()) {
            return structural->mass(storage);
        }

        MapMatrix matrix{storage, 0, 0};
        return matrix;
    };

    SparseMatrix matrix = mattools::assemble_matrix(_data->elements, indices, lambda);

    if (!_data->point_elements.empty()) {
        matrix += mattools::assemble_matrix(_data->point_elements, indices, lambda);
    }

    return matrix;
}

} // namespace model
} // namespace fem
