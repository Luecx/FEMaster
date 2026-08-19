/**
 * @file model_build.cpp
 * @brief Implements model-level field construction and structural assembly.
 *
 * The routines operate on the dense assembly created by `Model::compile()`.
 * They bind sections to compiled elements, construct shared shell reference
 * normals, enumerate active generalized DOFs, assemble loads and constraint
 * equations, and combine element, contact and feature contributions into global
 * sparse operators and nodal internal-force fields.
 *
 * Element-local matrix evaluation remains the responsibility of each structural
 * formulation. `mattools::assemble_matrix()` owns local-to-global DOF mapping,
 * sparse triplet accumulation and parallel reduction; this file supplies the
 * formulation-specific callbacks and adds model-level contributions that are not
 * represented by regular elements.
 *
 * @see Model
 * @see Model::compile
 * @see mattools::assemble_matrix
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../mattools/assemble.h"
#include "../mattools/numerate_dofs.h"
#include "element/element_structural.h"
#include "model.h"
#include "solid/element_solid.h"

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
 * assigned material and section properties.
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
    const Index expected_rows      = static_cast<Index>((*_data->element_nodal_offsets)(element_count));

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
        logging::error(static_cast<Index>(natural_coords.rows()) == n_nodes
                    && static_cast<Index>(natural_coords.cols()) == 2,
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
 * Element DOF masks activate the translational and rotational components used by
 * their formulations. Couplings may additionally require master-node DOFs, and
 * connectors activate their constrained components at both endpoint nodes. The
 * final boolean mask is converted to contiguous zero-based system indices;
 * inactive components receive negative identifiers.
 *
 * @return Node-by-six matrix of global unconstrained system DOF identifiers.
 */
SystemDofIds Model::build_unconstrained_index_matrix() {
    // Validate the compiled nodal dimension and initialize an inactive DOF mask
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");

    SystemDofs mask{_data->positions->rows, 6};
    mask.fill(false);

    // Activate every generalized DOF used by compiled element formulations
    for (auto& element : _data->elements) {
        if (element == nullptr) continue;

        const auto dofs = element->dofs();
        for (ID local_node = 0; local_node < element->n_nodes(); ++local_node) {
            const ID node_id = element->nodes()[local_node];
            for (ID dof = 0; dof < 6; ++dof) {
                mask(node_id, dof) |= dofs(0, dof);
            }
        }
    }

    // Add master-node components required to represent coupling kinematics
    for (auto& coupling : _data->couplings) {
        const ID master_id = coupling.master_node;
        const auto master_dofs = coupling.master_dofs(mask, *_data);
        for (ID dof = 0; dof < 6; ++dof) {
            mask(master_id, dof) |= master_dofs(0, dof);
        }
    }

    // Activate constrained components at both nodes of every connector
    for (auto& connector : _data->connectors) {
        const ID node1_id = connector.node_1();
        const ID node2_id = connector.node_2();
        const auto dofs   = connector.dofs();

        for (ID dof = 0; dof < 6; ++dof) {
            mask(node1_id, dof) |= dofs(0, dof);
            mask(node2_id, dof) |= dofs(0, dof);
        }
    }

    // Convert the boolean activity mask to contiguous global equation numbers
    return mattools::numerate_dofs(mask);
}

/**
 * Assembles selected load collectors into one global nodal load field.
 *
 * Each named collector evaluates its loads at `time`, including any attached
 * amplitude and coordinate-system transformation, into a six-component nodal
 * field. Couplings then redistribute generalized master-node loads to their
 * slave topology while preserving the supported resultant.
 *
 * @param load_sets Names of load collectors applied to the field.
 * @param time Evaluation time supplied to load amplitudes.
 * @return Six-component global nodal load field.
 */
Field Model::build_load_matrix(std::vector<std::string> load_sets, Precision time) {
    // Allocate a zeroed generalized load row for every compiled node
    Field load_matrix{"LOAD_MATRIX", FieldDomain::NODE, _data->field_rows(FieldDomain::NODE), 6};
    load_matrix.set_zero();

    // Accumulate all requested collectors at the supplied physical time
    for (auto& key : load_sets) {
        auto data = _data->load_cols.get(key);
        data->apply(*_data, load_matrix, time);
    }

    // Replace supported master loads by their statically equivalent slave loads
    for (auto& coupling : _data->couplings) {
        coupling.apply_loads(*_data, load_matrix);
    }

    return load_matrix;
}

/**
 * Decomposes selected loads into amplitude-independent spatial basis fields.
 *
 * Loads sharing the same amplitude pointer contribute to one six-component
 * nodal basis field; loads without an amplitude share the null-amplitude entry.
 * Each load is evaluated with amplitude scaling disabled so a later transient or
 * harmonic procedure can multiply the returned spatial field by the associated
 * scalar history. Coupling load redistribution is applied independently to every
 * completed basis field.
 *
 * @param load_sets Names of load collectors included in the decomposition.
 * @return Pairs of amplitude pointers and their unscaled global nodal load fields.
 */
std::vector<std::pair<bc::Amplitude::Ptr, Field>>
Model::build_load_basis(std::vector<std::string> load_sets) {
    std::vector<std::pair<bc::Amplitude::Ptr, Field>> basis;

    // Reuse one spatial field for every distinct shared amplitude definition
    auto field_for = [this, &basis](const bc::Amplitude::Ptr& amplitude) -> Field& {
        for (auto& entry : basis) {
            if (entry.first == amplitude) return entry.second;
        }

        Field load_matrix{"LOAD_BASIS", FieldDomain::NODE, _data->field_rows(FieldDomain::NODE), 6};
        load_matrix.set_zero();
        basis.emplace_back(amplitude, std::move(load_matrix));
        return basis.back().second;
    };

    // Assemble nominal load magnitudes without applying their scalar histories
    for (auto& key : load_sets) {
        auto data = _data->load_cols.get(key);
        for (const auto& load : data->entries()) {
            if (!load) continue;
            auto& load_matrix = field_for(load->amplitude_);
            load->apply(*_data, load_matrix, Precision(0), true);
        }
    }

    // Redistribute master-node loads independently in every amplitude basis
    for (auto& entry : basis) {
        for (auto& coupling : _data->couplings) {
            coupling.apply_loads(*_data, entry.second);
        }
    }

    return basis;
}

/**
 * Generates and categorizes all linear constraint equations for the model.
 *
 * When no support collector names are supplied, the aggregate support collection
 * is used when available; otherwise only the named collectors participate.
 * Connector, coupling, tie and rigid-body-mode objects generate their equations
 * in model order. Manually stored equations are copied into the fallback group.
 *
 * Every generated row receives source-kind and source-index metadata before it
 * is moved into `ConstraintGroups`. Contact is excluded because its unilateral
 * residual and tangent are assembled directly during nonlinear evaluation.
 *
 * @param system_dof_ids Active global DOF identifiers used by constraint builders.
 * @param supp_sets Optional names of support collectors to include.
 * @return Constraint equations grouped and annotated by their origin.
 */
constraint::ConstraintGroups Model::collect_constraints(
    SystemDofIds&                   system_dof_ids,
    const std::vector<std::string>& supp_sets
) {
    constraint::ConstraintGroups groups{};

    // Generate support equations from the aggregate or explicitly selected sets
    Index support_idx = 0;
    if (supp_sets.empty() && _data->supp_cols.has_all()) {
        if (auto all = _data->supp_cols.all()) {
            auto eqs = all->get_equations(*_data);
            for (auto& eq : eqs) {
                eq.source       = constraint::EquationSourceKind::Support;
                eq.source_index = support_idx;
                groups.supports.push_back(std::move(eq));
            }
            ++support_idx;
        }
    }

    for (const auto& key : supp_sets) {
        if (auto data = _data->supp_cols.get(key)) {
            auto eqs = data->get_equations(*_data);
            for (auto& eq : eqs) {
                eq.source       = constraint::EquationSourceKind::Support;
                eq.source_index = support_idx;
                groups.supports.push_back(std::move(eq));
            }
            ++support_idx;
        }
    }

    // Generate connector equations and retain their originating object indices
    Index connector_idx = 0;
    for (auto& connector : _data->connectors) {
        auto eqs = connector.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source       = constraint::EquationSourceKind::Connector;
            eq.source_index = connector_idx;
            groups.connectors.push_back(std::move(eq));
        }
        ++connector_idx;
    }

    // Generate kinematic or structural coupling equations
    Index coupling_idx = 0;
    for (auto& coupling : _data->couplings) {
        auto eqs = coupling.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source       = constraint::EquationSourceKind::Coupling;
            eq.source_index = coupling_idx;
            groups.couplings.push_back(std::move(eq));
        }
        ++coupling_idx;
    }

    // Generate line/surface tie equations from current compiled geometry
    Index tie_idx = 0;
    for (auto& tie : _data->ties) {
        auto eqs = tie.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source       = constraint::EquationSourceKind::Tie;
            eq.source_index = tie_idx;
            groups.ties.push_back(std::move(eq));
        }
        ++tie_idx;
    }

    // Generate rigid-body-mode suppression equations
    Index rbm_idx = 0;
    for (auto& rbm : _data->rbms) {
        auto eqs = rbm.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source       = constraint::EquationSourceKind::Rbm;
            eq.source_index = rbm_idx;
            groups.rbms.push_back(std::move(eq));
        }
        ++rbm_idx;
    }

    // Preserve explicit equations as manual or otherwise categorized rows
    if (!_data->equations.empty()) {
        Index manual_idx = 0;
        for (auto eq : _data->equations) {
            if (eq.source == constraint::EquationSourceKind::Unknown) {
                eq.source = constraint::EquationSourceKind::Manual;
            }
            eq.source_index = manual_idx++;
            groups.others.push_back(std::move(eq));
        }
    }

    return groups;
}

/**
 * Assembles the linear global structural stiffness matrix.
 *
 * Every structural element contributes its formulation-specific local stiffness
 * through the common sparse assembly pipeline. An optional scalar
 * `ELEMENT`-domain field scales each local matrix before global insertion.
 * Non-element features append their stiffness triplets afterwards.
 *
 * Unilateral contact is intentionally rejected because its active set, residual
 * and tangent depend on the current nonlinear configuration.
 *
 * @param indices Active global DOF identifiers used for local-to-global mapping.
 * @param stiffness_scalar Optional one-component element stiffness scale.
 * @return Assembled sparse linear stiffness matrix.
 */
SparseMatrix Model::build_stiffness_matrix(SystemDofIds& indices, const Field* stiffness_scalar) {
    // Reject configuration-dependent contact from the linear assembly path
    logging::error(_data->contacts.empty(),
        "CONTACT requires NONLINEARSTATIC; linear stiffness assembly cannot include contact");

    // Define the structural element callback and optional element-wise scaling
    auto lambda = [&](const ElementPtr& element, Precision* storage) {
        if (auto structural = element->as<StructuralElement>()) {
            MapMatrix stiffness = structural->stiffness(storage);
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

    // Assemble all element-local matrices over the active global DOFs
    SparseMatrix matrix = mattools::assemble_matrix(_data->elements, indices, lambda);

    // Append stiffness contributions represented by non-element features
    if (!_data->features.empty()) {
        TripletList triplets;
        for (const auto& feature : _data->features) {
            if (feature) feature->assemble_stiffness(indices, triplets);
        }
        if (!triplets.empty()) matrix.insertFromTriplets(triplets.begin(), triplets.end());
    }

    return matrix;
}

/**
 * Assembles the nonlinear tangent and matching nodal internal-force vector.
 *
 * Structural elements evaluate their consistent tangent in the current
 * displacement configuration while accumulating generalized internal forces in
 * the supplied nodal field. A temporary integration-point stress-state field
 * provides shared scratch storage required by the element interfaces. Optional
 * element stiffness scaling applies to the structural tangent matrices.
 *
 * Contact contributes its current dual-mortar residual directly to
 * `nodal_forces` and supplies frozen-geometry tangent triplets. Non-element
 * features append their stiffness contributions after contact assembly. The
 * completed force field is checked for finite values.
 *
 * @param indices Active global DOF identifiers used for sparse assembly.
 * @param nodal_forces Six-component nodal field receiving internal and contact
 *                     force contributions.
 * @param displacement Current global nodal displacement field.
 * @param stiffness_scalar Optional one-component element tangent scale.
 * @return Assembled sparse nonlinear tangent matrix.
 */
SparseMatrix Model::build_tangent_stiffness_matrix(
    SystemDofIds& indices,
    NodeData&     nodal_forces,
    const Field&  displacement,
    const Field*  stiffness_scalar
) {
    // Validate the nodal force output field before parallel accumulation
    logging::error(nodal_forces.domain == FieldDomain::NODE,
        "tangent internal force output must use NODE domain");
    logging::error(nodal_forces.rows == _data->field_rows(FieldDomain::NODE),
        "tangent internal force output has wrong node count");
    logging::error(nodal_forces.components >= 6,
        "tangent internal force output requires at least 6 components");

    // Allocate zeroed integration-point scratch storage for element evaluation
    Field ip_stress_state{
        "IP_STRESS_STATE", FieldDomain::ELEMENT_IP,
        _data->field_rows(FieldDomain::ELEMENT_IP), 8};
    ip_stress_state.set_zero();

    // Define the structural tangent/internal-force callback used by sparse assembly
    auto lambda = [&](const ElementPtr& element,
                      Precision*        local_matrix_storage,
                      NodeData&         local_nodal_forces) -> MapMatrix {
        auto* structural = element->as<StructuralElement>();
        if (!structural) {
            MapMatrix matrix{local_matrix_storage, 0, 0};
            return matrix;
        }

        MapMatrix tangent = structural->stiffness_tangent(
            local_matrix_storage,
            ip_stress_state,
            local_nodal_forces,
            displacement
        );

        if (stiffness_scalar) {
            logging::error(stiffness_scalar->domain == FieldDomain::ELEMENT,
                "stiffness scale field must use ELEMENT domain");
            logging::error(stiffness_scalar->components == 1,
                "stiffness scale field must have 1 component");
            tangent *= (*stiffness_scalar)(static_cast<Index>(structural->elem_id), 0);
        }

        return tangent;
    };

    // Assemble structural tangents and nodal forces with the common parallel path
    SparseMatrix global_matrix = mattools::assemble_matrix(
        _data->elements,
        indices,
        lambda,
        &nodal_forces
    );

    // Assemble configuration-dependent contact residuals and tangent triplets
    TripletList contact_triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, contact_triplets);
    }

    if (!contact_triplets.empty()) {
        SparseMatrix contact_matrix(global_matrix.rows(), global_matrix.cols());
        contact_matrix.insertFromTriplets(contact_triplets.begin(), contact_triplets.end());
        global_matrix += contact_matrix;
    }

    // Append stiffness contributions represented by non-element features
    if (!_data->features.empty()) {
        TripletList feature_triplets;
        for (const auto& feature : _data->features) {
            if (feature) feature->assemble_stiffness(indices, feature_triplets);
        }
        if (!feature_triplets.empty()) {
            global_matrix.insertFromTriplets(feature_triplets.begin(), feature_triplets.end());
        }
    }

    // Reject invalid structural or contact force contributions
    nodal_forces.check_finite("Internal force");

    return global_matrix;
}

/**
 * Evaluates nonlinear internal and contact forces without retaining a tangent.
 *
 * The output field is reset before structural elements accumulate their current
 * generalized internal forces. Contact then evaluates the same residual path as
 * tangent assembly, but its generated tangent triplets are deliberately
 * discarded. This supports residual-only evaluations used by nonlinear control
 * and line-search algorithms.
 *
 * @param indices Active global DOF identifiers required by contact assembly.
 * @param nodal_forces Six-component nodal field overwritten with the resulting
 *                     internal and contact forces.
 * @param displacement Current global nodal displacement field.
 */
void Model::build_internal_force_nonlinear(
    SystemDofIds& indices,
    NodeData&     nodal_forces,
    const Field&  displacement
) {
    // Validate the nodal force output field
    logging::error(nodal_forces.domain == FieldDomain::NODE,
        "nonlinear internal force output must use NODE domain");
    logging::error(nodal_forces.rows == _data->field_rows(FieldDomain::NODE),
        "nonlinear internal force output has wrong node count");
    logging::error(nodal_forces.components >= 6,
        "nonlinear internal force output requires at least 6 components");

    // Initialize element scratch state and clear all previous force contributions
    Field ip_stress_state{
        "IP_STRESS_STATE", FieldDomain::ELEMENT_IP,
        _data->field_rows(FieldDomain::ELEMENT_IP), 8};
    ip_stress_state.set_zero();
    nodal_forces.set_zero();

    // Accumulate nonlinear internal forces from every structural element
    for (const auto& element : _data->elements) {
        if (!element) continue;

        auto* structural = element->as<StructuralElement>();
        if (!structural) continue;

        structural->internal_force_nonlinear(
            ip_stress_state,
            nodal_forces,
            displacement
        );
    }

    // Evaluate contact residuals while intentionally discarding tangent entries
    TripletList discarded_contact_triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, discarded_contact_triplets);
    }

    // Reject invalid structural or contact force contributions
    nodal_forces.check_finite("Internal force");
}

/**
 * Assembles the global geometric stiffness matrix from integration-point stress.
 *
 * The compiled element-IP offset field identifies the first stress row belonging
 * to each structural element. Formulations use that stress state to compute
 * their local initial-stress/geometric matrix. An optional scalar
 * `ELEMENT`-domain field scales each local contribution before sparse assembly.
 *
 * Contact and non-element features do not contribute to this dedicated
 * geometric operator.
 *
 * @param indices Active global DOF identifiers used for sparse assembly.
 * @param ip_stress Integration-point stress field driving geometric stiffness.
 * @param stiffness_scalar Optional one-component element stiffness scale.
 * @return Assembled sparse geometric stiffness matrix.
 */
SparseMatrix Model::build_geom_stiffness_matrix(
    SystemDofIds& indices,
    const Field&  ip_stress,
    const Field*  stiffness_scalar
) {
    // Validate and access the compiled integration-point enumeration
    logging::error(_data->element_ip_offsets != nullptr,
        "element IP offset field has not been initialized");

    const Field& ip_enum = *_data->element_ip_offsets;

    // Define element-local geometric stiffness and optional scaling
    auto lambda = [&](const ElementPtr& element, Precision* storage) -> MapMatrix {
        if (auto structural = element->as<StructuralElement>()) {
            const ID element_id = structural->elem_id;
            const ID ip_start = static_cast<ID>(ip_enum(static_cast<Index>(element_id), 0));

            MapMatrix geometric_stiffness =
                structural->stiffness_geom(storage, ip_stress, ip_start);

            if (stiffness_scalar) {
                logging::error(stiffness_scalar->domain == FieldDomain::ELEMENT,
                    "stiffness scale field must use ELEMENT domain");
                logging::error(stiffness_scalar->components == 1,
                    "stiffness scale field must have 1 component");
                geometric_stiffness *=
                    (*stiffness_scalar)(static_cast<Index>(element_id), 0);
            }

            return geometric_stiffness;
        }

        MapMatrix matrix{storage, 0, 0};
        return matrix;
    };

    // Assemble all element-local geometric matrices over active global DOFs
    return mattools::assemble_matrix(_data->elements, indices, lambda);
}

/**
 * Assembles the global lumped mass matrix.
 *
 * Structural elements provide their local lumped translational and rotational
 * mass matrices through the common sparse assembly path. Non-element features,
 * including concentrated point masses and rotary inertias, append additional
 * diagonal mass triplets afterwards.
 *
 * @param indices Active global DOF identifiers used for sparse assembly.
 * @return Assembled sparse lumped mass matrix.
 */
SparseMatrix Model::build_lumped_mass_matrix(SystemDofIds& indices) {
    // Define the structural element mass callback
    auto lambda = [&](const ElementPtr& element, Precision* storage) {
        if (auto structural = element->as<StructuralElement>()) {
            return structural->mass(storage);
        }

        MapMatrix matrix{storage, 0, 0};
        return matrix;
    };

    // Assemble all structural element mass contributions
    SparseMatrix matrix = mattools::assemble_matrix(_data->elements, indices, lambda);

    // Append concentrated mass contributions represented by model features
    if (!_data->features.empty()) {
        TripletList triplets;
        for (const auto& feature : _data->features) {
            if (feature) feature->assemble_mass(indices, triplets);
        }
        if (!triplets.empty()) matrix.insertFromTriplets(triplets.begin(), triplets.end());
    }

    return matrix;
}

} // namespace model
} // namespace fem
