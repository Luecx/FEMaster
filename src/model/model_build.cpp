/**
 * @file model_build.cpp
 * @brief Implements model-level field construction and structural assembly.
 *
 * This file contains build operations that combine element-local information
 * into model-wide fields, constraints, forces and sparse matrices. Element
 * formulations remain responsible for local kinematics, constitutive response
 * and element matrices, while `Model` coordinates global enumeration, field
 * ownership and sparse assembly.
 *
 * Nonlinear structural assembly combines element contributions with feature and
 * surface-to-surface mortar contact contributions. Contact adds both nodal
 * internal forces and tangent triplets after the structural element assembly so
 * the nonlinear residual and tangent use the same current contact state.
 *
 * Shell reference directors are represented by an `ELEMENT_NODAL` field. The
 * shell-normal build routine completes missing field rows from reference geometry
 * and equalises compatible element-local normals at shared nodes while preserving
 * explicitly prescribed values.
 *
 * @see Model
 * @see ModelData
 * @see StructuralElement
 * @see constraint::Contact
 *
 * @author Finn Eggers
 * @date 10.08.2026
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
 * Assigns every configured section to all elements in the section's region.
 *
 * Section objects own their target element regions. Null element slots are
 * ignored so sparse element id spaces remain valid.
 */
void Model::assign_sections() {
    for (Section::Ptr ptr : this->_data->sections) {
        for (ID elem_id : *ptr->region_) {
            if (this->_data->elements[elem_id] != nullptr) {
                this->_data->elements[elem_id]->set_section(ptr);
            }
        }
    }
}

/**
 * Builds and completes the element-nodal reference-normal field for shells.
 *
 * If `ModelData::shell_element_nodal_normals` already references a field, every
 * finite row with magnitude above the numerical tolerance is interpreted as a
 * prescribed element-node normal. Prescribed values remain unchanged in field
 * storage. A normalized copy is used only for angular comparison and averaging.
 *
 * Missing rows are evaluated from the physical reference surface at the
 * corresponding natural node coordinate. At each global node, the resulting
 * unit vectors are divided into angularly compatible clusters according to
 *
 * \f[
 *     \hat{\boldsymbol n}_i \cdot \hat{\boldsymbol n}_c
 *     \ge \cos(\theta_{\mathrm{eq}}).
 * \f]
 *
 * Prescribed and geometrically evaluated vectors contribute equally to the
 * normalized cluster mean
 *
 * \f[
 *     \bar{\boldsymbol n}
 *     = \frac{\sum_i \hat{\boldsymbol n}_i}
 *            {\left\|\sum_i \hat{\boldsymbol n}_i\right\|}.
 * \f]
 *
 * The cluster mean is written only to previously missing rows. Consequently,
 * two prescribed normals and one missing normal at the same smooth node cause
 * the missing row to follow the mean of all three directions, while both
 * prescribed field entries remain untouched. Angularly incompatible clusters
 * remain separate and therefore preserve sharp folds.
 *
 * If no field was selected, the routine creates the default
 * `SHELL_ELEMENT_NODAL_NORMALS` field and reproduces the fully automatic
 * construction used previously.
 *
 * @param equalize_angle_degrees Maximum angle between normals that may belong
 *                               to one smooth cluster at a global node.
 */
void Model::build_shell_element_normals(Precision equalize_angle_degrees) {
    // ---------------------------------------------------------------------
    // Validate model state and configure angular equalisation
    // ---------------------------------------------------------------------
    logging::error(_data->positions_reference != nullptr,
        "Model: POSITION_REFERENCE field is not set");
    logging::error(_data->element_nodal_offsets != nullptr,
        "Model: element nodal offsets are not initialized");

    constexpr Precision normal_tolerance = Precision(1e-12);

    const Precision pi             = std::acos(Precision(-1));
    const Precision equalize_angle = equalize_angle_degrees * pi / Precision(180);
    const Precision cos_equalize   = std::cos(equalize_angle);
    const Index expected_rows =
        static_cast<Index>((*_data->element_nodal_offsets)(static_cast<Index>(_data->max_elems)));

    // ---------------------------------------------------------------------
    // Reuse a selected field or create the default automatic normal field
    // ---------------------------------------------------------------------
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

    // Store the element-nodal rows attached to each global node. The working
    // vectors are normalized independently of the original field values, while
    // the prescribed mask controls which rows may be written during completion.
    std::vector<std::vector<Index>> node_rows(static_cast<std::size_t>(_data->max_nodes));
    std::vector<Vec3>               row_normals(static_cast<std::size_t>(normals->rows), Vec3::Zero());
    std::vector<bool>               prescribed_rows(static_cast<std::size_t>(normals->rows), false);

    // ---------------------------------------------------------------------
    // Gather prescribed normals and evaluate missing geometric normals
    // ---------------------------------------------------------------------
    for (const ElementPtr& element : _data->elements) {
        if (element == nullptr) {
            continue;
        }

        const auto structural = element->as<StructuralElement>();
        if (structural == nullptr || !structural->is_shell()) {
            continue;
        }

        const auto surface = structural->surface(1);
        logging::error(surface != nullptr,
            "Model: shell element ", structural->elem_id,
            " does not provide a reference surface");

        const Index offset  = static_cast<Index>(structural->elem_nodal_offset);
        const Index n_nodes = static_cast<Index>(structural->n_nodes());

        const DynamicMatrix natural_coords = surface->node_coords_natural();
        logging::error(static_cast<Index>(natural_coords.rows()) == n_nodes
                    && static_cast<Index>(natural_coords.cols()) == 2,
            "Model: natural surface node coordinates do not match shell element ",
            structural->elem_id);

        for (Index local_node = 0; local_node < n_nodes; ++local_node) {
            const ID    node_id = structural->nodes()[local_node];
            const Index row     = offset + local_node;

            // A valid field vector is prescribed. Normalize only the local copy
            // so the user-provided field row remains exactly unchanged.
            Vec3 normal = normals->row_vec3(row);
            const bool prescribed = normal.allFinite() && normal.norm() > normal_tolerance;

            if (prescribed) {
                normal.normalize();
                prescribed_rows[static_cast<std::size_t>(row)] = true;
            } else {
                // Evaluate the physical reference normal at the natural shell
                // node only when no usable field value was supplied.
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

    // ---------------------------------------------------------------------
    // Cluster compatible normals and complete only non-prescribed rows
    // ---------------------------------------------------------------------
    for (const auto& rows : node_rows) {
        std::vector<Vec3>               cluster_sums;
        std::vector<std::vector<Index>> cluster_rows;

        // Build angular clusters using the current normalized sum as the
        // representative direction. Every member contributes unit weight.
        for (Index row : rows) {
            const Vec3 normal = row_normals[static_cast<std::size_t>(row)];
            bool added = false;

            for (Index cluster = 0; cluster < static_cast<Index>(cluster_sums.size()); ++cluster) {
                const Vec3 cluster_normal =
                    cluster_sums[static_cast<std::size_t>(cluster)].normalized();

                if (normal.dot(cluster_normal) < cos_equalize) {
                    continue;
                }

                cluster_sums[static_cast<std::size_t>(cluster)] += normal;
                cluster_rows[static_cast<std::size_t>(cluster)].push_back(row);
                added = true;
                break;
            }

            if (!added) {
                cluster_sums.push_back(normal);
                cluster_rows.push_back({row});
            }
        }

        // The normalized sum is the equalised direction of one smooth cluster.
        // Prescribed rows participated in this sum but are deliberately skipped
        // during the write-back step.
        for (Index cluster = 0; cluster < static_cast<Index>(cluster_sums.size()); ++cluster) {
            const Vec3 equalized_normal =
                cluster_sums[static_cast<std::size_t>(cluster)].normalized();

            for (Index row : cluster_rows[static_cast<std::size_t>(cluster)]) {
                if (prescribed_rows[static_cast<std::size_t>(row)]) {
                    continue;
                }

                (*normals)(row, 0) = equalized_normal(0);
                (*normals)(row, 1) = equalized_normal(1);
                (*normals)(row, 2) = equalized_normal(2);
            }
        }
    }

    // Publish either the completed selected field or the newly created default
    // field as the semantic shell-normal storage consumed by shell elements.
    _data->shell_element_nodal_normals = normals;
}

/**
 * Builds the global active-degree-of-freedom index matrix.
 *
 * Structural element DOFs form the initial mask. Coupling master DOFs and both
 * connector-node DOF sets are then added because those algebraic constraints may
 * introduce active unknowns that are not directly requested by an element.
 * `numerate_dofs()` converts the final boolean mask to contiguous system ids.
 */
SystemDofIds Model::build_unconstrained_index_matrix() {
    SystemDofs mask{this->_data->max_nodes, 6};
    mask.fill(false);

    // Gather DOFs requested by all existing elements.
    for (auto& element : _data->elements) {
        if (element == nullptr) {
            continue;
        }

        const auto dofs = element->dofs();
        for (ID local_node = 0; local_node < element->n_nodes(); ++local_node) {
            const ID node_id = element->nodes()[local_node];
            for (ID dof = 0; dof < 6; ++dof) {
                mask(node_id, dof) |= dofs(0, dof);
            }
        }
    }

    // Coupling masters may require additional translational or rotational DOFs.
    for (auto& coupling : this->_data->couplings) {
        const ID master_id = coupling.master_node;
        const auto master_dofs = coupling.master_dofs(mask, *_data);
        for (ID dof = 0; dof < 6; ++dof) {
            mask(master_id, dof) |= master_dofs(0, dof);
        }
    }

    // Connector equations operate on both connector nodes.
    for (auto& connector : this->_data->connectors) {
        const ID node1_id = connector.node_1();
        const ID node2_id = connector.node_2();
        const auto dofs   = connector.dofs();

        for (ID dof = 0; dof < 6; ++dof) {
            mask(node1_id, dof) |= dofs(0, dof);
            mask(node2_id, dof) |= dofs(0, dof);
        }
    }

    return mattools::numerate_dofs(mask);
}

/**
 * Builds the total nodal load field at one analysis time.
 *
 * Selected load collectors are evaluated first. Coupling transformations are
 * then applied so loads acting on coupled slave regions are represented on the
 * algebraic master/slave DOFs consistently with the constraint equations.
 *
 * @param load_sets Load-collector names to include.
 * @param time Analysis time supplied to amplitude interpolation.
 * @return Six-component nodal load field.
 */
Field Model::build_load_matrix(std::vector<std::string> load_sets, Precision time) {
    Field load_matrix{"LOAD_MATRIX", FieldDomain::NODE, static_cast<Index>(this->_data->max_nodes), 6};
    load_matrix.set_zero();

    for (auto& key : load_sets) {
        auto data = this->_data->load_cols.get(key);
        data->apply(*_data, load_matrix, time);
    }

    for (auto& coupling : this->_data->couplings) {
        coupling.apply_loads(*_data, load_matrix);
    }

    return load_matrix;
}

/**
 * Builds amplitude-separated nodal load basis fields.
 *
 * Loads sharing the same amplitude pointer are accumulated into one basis field.
 * The stored field represents the unscaled spatial load pattern; time-dependent
 * amplitude factors are applied by the transient/nonlinear caller. Coupling load
 * transformations are applied independently to each basis field.
 *
 * @param load_sets Load-collector names to include.
 * @return Pairs of amplitude handles and corresponding nodal load patterns.
 */
std::vector<std::pair<bc::Amplitude::Ptr, Field>>
Model::build_load_basis(std::vector<std::string> load_sets) {
    std::vector<std::pair<bc::Amplitude::Ptr, Field>> basis;

    auto field_for = [this, &basis](const bc::Amplitude::Ptr& amplitude) -> Field& {
        for (auto& entry : basis) {
            if (entry.first == amplitude) {
                return entry.second;
            }
        }

        Field load_matrix{"LOAD_BASIS", FieldDomain::NODE, static_cast<Index>(this->_data->max_nodes), 6};
        load_matrix.set_zero();
        basis.emplace_back(amplitude, std::move(load_matrix));
        return basis.back().second;
    };

    for (auto& key : load_sets) {
        auto data = this->_data->load_cols.get(key);
        for (const auto& load : data->entries()) {
            if (!load) {
                continue;
            }

            auto amplitude = load->amplitude_;
            auto& load_matrix = field_for(amplitude);
            load->apply(*_data, load_matrix, Precision(0), true);
        }
    }

    for (auto& entry : basis) {
        for (auto& coupling : this->_data->couplings) {
            coupling.apply_loads(*_data, entry.second);
        }
    }

    return basis;
}

/**
 * Collects all active linear constraints while preserving their semantic source.
 *
 * Support, connector, coupling, tie, rigid-body and manual equations are gathered
 * into separate groups and annotated with source indices. This grouped form is
 * used for diagnostics before `ConstraintGroups::flatten()` produces the algebraic
 * equation list consumed by the constraint transformer.
 *
 * @param system_dof_ids Active global DOF ids required by connector/coupling
 *                       equation generation.
 * @param supp_sets Explicit support collectors; when empty, the global `ALL`
 *                  support collector is used when available.
 * @return Semantic constraint groups in model order.
 */
constraint::ConstraintGroups Model::collect_constraints(
    SystemDofIds&                   system_dof_ids,
    const std::vector<std::string>& supp_sets
) {
    constraint::ConstraintGroups groups{};

    Index support_idx = 0;
    if (supp_sets.empty() && this->_data->supp_cols.has_all()) {
        if (auto all = this->_data->supp_cols.all()) {
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
        if (auto data = this->_data->supp_cols.get(key)) {
            auto eqs = data->get_equations(*_data);
            for (auto& eq : eqs) {
                eq.source       = constraint::EquationSourceKind::Support;
                eq.source_index = support_idx;
                groups.supports.push_back(std::move(eq));
            }
            ++support_idx;
        }
    }

    Index connector_idx = 0;
    for (auto& connector : this->_data->connectors) {
        auto eqs = connector.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source       = constraint::EquationSourceKind::Connector;
            eq.source_index = connector_idx;
            groups.connectors.push_back(std::move(eq));
        }
        ++connector_idx;
    }

    Index coupling_idx = 0;
    for (auto& coupling : this->_data->couplings) {
        auto eqs = coupling.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source       = constraint::EquationSourceKind::Coupling;
            eq.source_index = coupling_idx;
            groups.couplings.push_back(std::move(eq));
        }
        ++coupling_idx;
    }

    Index tie_idx = 0;
    for (auto& tie : this->_data->ties) {
        auto eqs = tie.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source       = constraint::EquationSourceKind::Tie;
            eq.source_index = tie_idx;
            groups.ties.push_back(std::move(eq));
        }
        ++tie_idx;
    }

    Index rbm_idx = 0;
    for (auto& rbm : this->_data->rbms) {
        auto eqs = rbm.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source       = constraint::EquationSourceKind::Rbm;
            eq.source_index = rbm_idx;
            groups.rbms.push_back(std::move(eq));
        }
        ++rbm_idx;
    }

    if (!this->_data->equations.empty()) {
        Index manual_idx = 0;
        for (auto eq : this->_data->equations) {
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
 * Builds the flat linear equation set used by algebraic constraint handling.
 *
 * This is the compact compatibility entry point around `collect_constraints()`;
 * semantic grouping is discarded only after all source metadata has been added.
 */
constraint::Equations Model::build_constraints(
    SystemDofIds&             system_dof_ids,
    std::vector<std::string> supp_sets
) {
    return collect_constraints(system_dof_ids, supp_sets).flatten();
}

/**
 * Assembles the linear structural stiffness matrix.
 *
 * Contact is rejected explicitly because unilateral mortar contact requires the
 * current nonlinear configuration and is supported only by nonlinear tangent
 * assembly. Optional element scalar fields scale structural element stiffness
 * before global sparse assembly. Feature stiffness contributions are appended as
 * sparse triplets afterwards.
 *
 * @param indices Active global DOF ids.
 * @param stiffness_scalar Optional one-component element stiffness scale field.
 * @return Global sparse linear stiffness matrix.
 */
SparseMatrix Model::build_stiffness_matrix(SystemDofIds& indices, const Field* stiffness_scalar) {
    logging::error(_data->contacts.empty(),
        "CONTACT requires NONLINEARSTATIC; linear stiffness assembly cannot include contact");

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

    SparseMatrix matrix = mattools::assemble_matrix(_data->elements, indices, lambda);

    if (!_data->features.empty()) {
        TripletList triplets;
        for (const auto& feature : _data->features) {
            if (feature) {
                feature->assemble_stiffness(indices, triplets);
            }
        }
        if (!triplets.empty()) {
            matrix.insertFromTriplets(triplets.begin(), triplets.end());
        }
    }

    return matrix;
}

/**
 * Assembles the nonlinear structural/contact tangent and current internal force.
 *
 * Structural elements first evaluate their material and geometric tangent while
 * scattering element internal forces through the common parallel assembly path.
 * The element callback also fills the current integration-point stress/resultant
 * field used internally by nonlinear elements.
 *
 * Surface-to-surface mortar contact is assembled afterwards from the same current
 * `ModelData::positions`. Contact contributes directly to `nodal_forces` and
 * returns sparse tangent triplets which are converted to a matrix and added to
 * the structural tangent. Feature stiffness contributions are appended last.
 *
 * @param indices Active global DOF ids used by structural/contact scattering.
 * @param nodal_forces Nodal internal-force field receiving structural and contact
 *                     contributions.
 * @param displacement Current total displacement field used by nonlinear elements.
 * @param stiffness_scalar Optional one-component element stiffness scale field.
 * @return Global sparse tangent matrix for the current nonlinear configuration.
 */
SparseMatrix Model::build_tangent_stiffness_matrix(
    SystemDofIds& indices,
    NodeData&     nodal_forces,
    const Field&  displacement,
    const Field*  stiffness_scalar
) {
    // ---------------------------------------------------------------------
    // Validate nonlinear force output and prepare current IP result storage
    // ---------------------------------------------------------------------
    logging::error(nodal_forces.domain == FieldDomain::NODE,
        "tangent internal force output must use NODE domain");
    logging::error(nodal_forces.rows == static_cast<Index>(_data->max_nodes),
        "tangent internal force output has wrong node count");
    logging::error(nodal_forces.components >= 6,
        "tangent internal force output requires at least 6 components");

    const Index max_ip = this->_data->max_integration_points;

    Field ip_stress_state{"IP_STRESS_STATE", FieldDomain::ELEMENT_IP, max_ip, 8};
    ip_stress_state.set_zero();

    // ---------------------------------------------------------------------
    // Assemble nonlinear structural element tangent and internal force
    // ---------------------------------------------------------------------
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

    SparseMatrix global_matrix = mattools::assemble_matrix(
        _data->elements,
        indices,
        lambda,
        &nodal_forces
    );

    // ---------------------------------------------------------------------
    // Add current surface-to-surface mortar residual and tangent
    // ---------------------------------------------------------------------
    TripletList contact_triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, contact_triplets);
    }

    if (!contact_triplets.empty()) {
        SparseMatrix contact_matrix(global_matrix.rows(), global_matrix.cols());
        contact_matrix.insertFromTriplets(contact_triplets.begin(), contact_triplets.end());
        global_matrix += contact_matrix;
    }

    // ---------------------------------------------------------------------
    // Add feature stiffness and validate assembled nodal internal forces
    // ---------------------------------------------------------------------
    if (!_data->features.empty()) {
        TripletList feature_triplets;
        for (const auto& feature : _data->features) {
            if (feature) {
                feature->assemble_stiffness(indices, feature_triplets);
            }
        }
        if (!feature_triplets.empty()) {
            global_matrix.insertFromTriplets(feature_triplets.begin(), feature_triplets.end());
        }
    }

    for (Index i = 0; i < nodal_forces.rows; ++i) {
        for (Index j = 0; j < nodal_forces.components; ++j) {
            const bool bad = std::isnan(nodal_forces(i, j)) || std::isinf(nodal_forces(i, j));
            logging::error(!bad,
                "Internal force at node ", i, " has invalid value at col ", j);
        }
    }

    return global_matrix;
}

/**
 * Assembles the nonlinear internal force at the current trial displacement
 * without assembling the structural tangent stiffness.
 *
 * The nonlinear Newton line search only needs the residual norm at temporary
 * trial states. This routine therefore evaluates each structural element's
 * current stress or shell resultant state and scatters `f_int` directly into the
 * supplied nodal field. Contact uses the same slave discretisation as tangent
 * assembly so line-search residuals remain identical to the residual associated
 * with the assembled tangent.
 *
 * @param indices Active global degree-of-freedom ids used by contact assembly.
 * @param nodal_forces Nodal internal-force field to overwrite and fill.
 * @param displacement Trial displacement defining the current configuration.
 */
void Model::build_internal_force_nonlinear(
    SystemDofIds& indices,
    NodeData&     nodal_forces,
    const Field&  displacement
) {
    // ---------------------------------------------------------------------
    // Validate output field and initialize current residual state
    // ---------------------------------------------------------------------
    logging::error(nodal_forces.domain == FieldDomain::NODE,
        "nonlinear internal force output must use NODE domain");
    logging::error(nodal_forces.rows == static_cast<Index>(_data->max_nodes),
        "nonlinear internal force output has wrong node count");
    logging::error(nodal_forces.components >= 6,
        "nonlinear internal force output requires at least 6 components");

    const Index max_ip = this->_data->max_integration_points;

    Field ip_stress_state{"IP_STRESS_STATE", FieldDomain::ELEMENT_IP, max_ip, 8};
    ip_stress_state.set_zero();
    nodal_forces.set_zero();

    // ---------------------------------------------------------------------
    // Assemble structural nonlinear internal force without element tangents
    // ---------------------------------------------------------------------
    for (const auto& element : _data->elements) {
        if (!element) {
            continue;
        }

        auto* structural = element->as<StructuralElement>();
        if (!structural) {
            continue;
        }

        structural->internal_force_nonlinear(
            ip_stress_state,
            nodal_forces,
            displacement
        );
    }

    // ---------------------------------------------------------------------
    // Add mortar contact residual and validate assembled force field
    // ---------------------------------------------------------------------
    TripletList discarded_contact_triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, discarded_contact_triplets);
    }

    for (Index i = 0; i < nodal_forces.rows; ++i) {
        for (Index j = 0; j < nodal_forces.components; ++j) {
            const bool bad =
                std::isnan(nodal_forces(i, j)) ||
                std::isinf(nodal_forces(i, j));

            logging::error(!bad,
                "Internal force at node ", i, " has invalid value at col ", j);
        }
    }
}

/**
 * Assembles the structural geometric stiffness from a supplied IP stress state.
 *
 * Element integration-point offsets map global `ELEMENT_IP` rows back to each
 * structural element. Optional element scalar values scale the resulting local
 * geometric stiffness before sparse assembly.
 */
SparseMatrix Model::build_geom_stiffness_matrix(
    SystemDofIds& indices,
    const Field&  ip_stress,
    const Field*  stiffness_scalar
) {
    logging::error(_data->element_ip_offsets != nullptr,
        "element IP offset field has not been initialized");

    const Field& ip_enum = *_data->element_ip_offsets;

    auto lambda = [&](const ElementPtr& element, Precision* storage) -> MapMatrix {
        if (auto structural = element->as<StructuralElement>()) {
            const ID element_id = structural->elem_id;
            const ID ip_start =
                static_cast<ID>(ip_enum(static_cast<Index>(element_id), 0));

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

    return mattools::assemble_matrix(_data->elements, indices, lambda);
}

/**
 * Assembles the global lumped structural mass matrix.
 *
 * Element-local lumped masses are assembled through the common sparse assembly
 * path. Feature mass contributions, such as point masses, are appended as sparse
 * triplets afterwards.
 */
SparseMatrix Model::build_lumped_mass_matrix(SystemDofIds& indices) {
    auto lambda = [&](const ElementPtr& element, Precision* storage) {
        if (auto structural = element->as<StructuralElement>()) {
            return structural->mass(storage);
        }

        MapMatrix matrix{storage, 0, 0};
        return matrix;
    };

    SparseMatrix matrix = mattools::assemble_matrix(_data->elements, indices, lambda);

    if (!_data->features.empty()) {
        TripletList triplets;
        for (const auto& feature : _data->features) {
            if (feature) {
                feature->assemble_mass(indices, triplets);
            }
        }
        if (!triplets.empty()) {
            matrix.insertFromTriplets(triplets.begin(), triplets.end());
        }
    }

    return matrix;
}

} // namespace model
} // namespace fem
