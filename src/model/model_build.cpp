/**
 * @file model_build.cpp
 * @brief Implements model-level field construction and structural assembly.
 *
 * This file contains build operations that combine element-local information
 * into model-wide fields, constraints, forces and sparse matrices. Element
 * formulations remain responsible for their local mechanics, while `Model`
 * coordinates enumeration, field ownership and global assembly.
 *
 * Shell reference directors are represented by an `ELEMENT_NODAL` field. The
 * shell-normal build routine completes missing field rows from the reference
 * geometry and equalises compatible element-local normals at shared nodes while
 * preserving explicitly prescribed values.
 *
 * @see Model
 * @see ModelData
 *
 * @author Finn Eggers
 * @date 30.07.2026
 */

#include "../mattools/assemble.h"
#include "../mattools/numerate_dofs.h"
#include "element/element_structural.h"
#include "model.h"
#include "solid/element_solid.h"

#include <cmath>
#include <iterator>
#include <vector>
#include <string>
#include <utility>

namespace fem {
namespace model {
void Model::assign_sections() {
    for (Section::Ptr ptr: this->_data->sections) {
        for (ID elem_id: *ptr->region_) {
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

SystemDofIds Model::build_unconstrained_index_matrix() {
    // first build a boolean matrix of which dofs are present in the system
    SystemDofs mask{this->_data->max_nodes, 6};
    mask.fill(false);

    // go through each elements and mask the dofs for all the nodes
    for (auto &e: _data->elements) {
        if (e != nullptr) {
            auto dofs = e->dofs();
            for (ID node_local_id = 0; node_local_id < e->n_nodes(); node_local_id++) {
                ID node_id = e->nodes()[node_local_id];
                for (ID dof = 0; dof < 6; dof++) {
                    mask(node_id, dof) |= dofs(0, dof);
                }
            }
        }
    }

    // go through all _couplings and mask the master dof
    for (auto &c: this->_data->couplings) {
        ID master_id = c.master_node;
        auto master_dofs = c.master_dofs(mask, *_data);
        for (ID dof = 0; dof < 6; dof++) {
            mask(master_id, dof) |= master_dofs(0, dof);
        }
    }
    // go through all connectors and mask the dofs of both nodes
    for (auto &c: this->_data->connectors) {
        ID node1_id = c.node_1();
        ID node2_id = c.node_2();
        auto dofs   = c.dofs();

        for (ID dof = 0; dof < 6; dof++) {
            mask(node1_id, dof) |= dofs(0, dof);
            mask(node2_id, dof) |= dofs(0, dof);
        }
    }

    auto res = mattools::numerate_dofs(mask);

    return res;
}

Field Model::build_load_matrix(std::vector<std::string> load_sets, Precision time) {
    Field load_matrix{"LOAD_MATRIX", FieldDomain::NODE, static_cast<Index>(this->_data->max_nodes), 6};
    load_matrix.set_zero();

    for (auto &key: load_sets) {
        auto data = this->_data->load_cols.get(key);
        data->apply(*_data, load_matrix, time);
    }

    // apply constrained loads
    for (auto &c: this->_data->couplings) {
        c.apply_loads(*_data, load_matrix);
    }

    return load_matrix;
}

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

    for (auto& key: load_sets) {
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
        for (auto& c: this->_data->couplings) {
            c.apply_loads(*_data, entry.second);
        }
    }

    return basis;
}

constraint::ConstraintGroups Model::collect_constraints(SystemDofIds& system_dof_ids,
                                                   const std::vector<std::string>& supp_sets) {
    constraint::ConstraintGroups groups{};

    Index support_idx = 0;
    if (supp_sets.empty() && this->_data->supp_cols.has_all()) {
        if (auto all = this->_data->supp_cols.all()) {
            auto eqs = all->get_equations(*_data);
            for (auto& eq : eqs) {
                eq.source = constraint::EquationSourceKind::Support;
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
                eq.source = constraint::EquationSourceKind::Support;
                eq.source_index = support_idx;
                groups.supports.push_back(std::move(eq));
            }
            ++support_idx;
        }
    }

    Index connector_idx = 0;
    for (auto& c : this->_data->connectors) {
        auto eqs = c.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source = constraint::EquationSourceKind::Connector;
            eq.source_index = connector_idx;
            groups.connectors.push_back(std::move(eq));
        }
        ++connector_idx;
    }

    Index coupling_idx = 0;
    for (auto& c : this->_data->couplings) {
        auto eqs = c.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source = constraint::EquationSourceKind::Coupling;
            eq.source_index = coupling_idx;
            groups.couplings.push_back(std::move(eq));
        }
        ++coupling_idx;
    }

    Index tie_idx = 0;
    for (auto& t : this->_data->ties) {
        auto eqs = t.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source = constraint::EquationSourceKind::Tie;
            eq.source_index = tie_idx;
            groups.ties.push_back(std::move(eq));
        }
        ++tie_idx;
    }

    Index rbm_idx = 0;
    for (auto& r : this->_data->rbms) {
        auto eqs = r.get_equations(system_dof_ids, *_data);
        for (auto& eq : eqs) {
            eq.source = constraint::EquationSourceKind::Rbm;
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

constraint::Equations Model::build_constraints(SystemDofIds& system_dof_ids,
                                               std::vector<std::string> supp_sets) {
    return collect_constraints(system_dof_ids, supp_sets).flatten();
}

SparseMatrix Model::build_stiffness_matrix(SystemDofIds &indices, const Field* stiffness_scalar) {
    logging::error(_data->contacts.empty(),
                   "CONTACT requires NONLINEARSTATIC; linear stiffness assembly cannot include contact");

    auto lambda = [&](const ElementPtr &el, Precision* storage) {
        if (auto sel = el->as<StructuralElement>()) {
            MapMatrix stiff = sel->stiffness(storage);
            if (stiffness_scalar) {
                logging::error(stiffness_scalar->domain == FieldDomain::ELEMENT,
                               "stiffness scale field must use ELEMENT domain");
                logging::error(stiffness_scalar->components == 1,
                               "stiffness scale field must have 1 component");
                stiff *= (*stiffness_scalar)(static_cast<Index>(sel->elem_id), 0);
            }
            return stiff;
        } else {
            MapMatrix mat{storage, 0, 0};
            return mat;
        }
    };
    auto res = mattools::assemble_matrix(_data->elements, indices, lambda);
    // Add feature-based stiffness contributions (diagonal triplets etc.)
    if (!_data->features.empty()) {
        TripletList trips;
        for (const auto& f : _data->features) if (f) f->assemble_stiffness(indices, trips);
        if (!trips.empty()) res.insertFromTriplets(trips.begin(), trips.end());
    }
    return res;
}

SparseMatrix Model::build_tangent_stiffness_matrix(SystemDofIds& indices,
                                                   NodeData& nodal_forces,
                                                   const Field& displacement,
                                                   const Field* stiffness_scalar) {
    logging::error(nodal_forces.domain == FieldDomain::NODE,
        "tangent internal force output must use NODE domain");
    logging::error(nodal_forces.rows == static_cast<Index>(_data->max_nodes),
        "tangent internal force output has wrong node count");
    logging::error(nodal_forces.components >= 6,
        "tangent internal force output requires at least 6 components");

    const Index max_ip = this->_data->max_integration_points;

    // Store the current integration-point resultants produced during tangent
    // evaluation. Element IP rows are element-local, so parallel element calls
    // write disjoint rows.
    Field ip_stress_state{"IP_STRESS_STATE", FieldDomain::ELEMENT_IP, max_ip, 8};
    ip_stress_state.set_zero();

    // Evaluate the element tangent and scatter internal force contributions into
    // the nodal force field owned by the active assembly worker.
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

    // Assemble structural element tangents through the common sparse-matrix
    // assembly path. With OpenMP enabled, nodal internal forces are accumulated in
    // thread-local fields and reduced after the element loop.
    SparseMatrix global_matrix = mattools::assemble_matrix(
        _data->elements,
        indices,
        lambda,
        &nodal_forces
    );

    // Add nonlinear contact contributions after the structural element matrix has
    // been assembled. Contact assembly owns additional search state and is kept
    // outside the generic element loop.
    TripletList contact_triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, contact_triplets);
    }

    if (!contact_triplets.empty()) {
        SparseMatrix contact_matrix(global_matrix.rows(), global_matrix.cols());
        contact_matrix.insertFromTriplets(contact_triplets.begin(), contact_triplets.end());
        global_matrix += contact_matrix;
    }

    // Add feature-based stiffness contributions such as point springs.
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

    // Validate the reduced internal-force assembly before the nonlinear solver
    // forms the residual.
    for (Index i = 0; i < nodal_forces.rows; ++i) {
        for (Index j = 0; j < nodal_forces.components; ++j) {
            const bool bad = std::isnan(nodal_forces(i, j)) || std::isinf(nodal_forces(i, j));
            logging::error(!bad, "Internal force at node ", i, " has invalid value at col ", j);
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
 * supplied nodal field. The same contact force contribution as tangent assembly
 * is retained so the residual equation remains unchanged. Contact is explicitly
 * invoked in force-only mode, avoiding construction of unused sparse tangent
 * entries during line-search evaluations.
 *
 * @param indices Active global degree-of-freedom ids used by contact assembly.
 * @param nodal_forces Nodal internal-force field to overwrite and fill.
 * @param displacement Trial displacement defining the current configuration.
 */
void Model::build_internal_force_nonlinear(SystemDofIds& indices,
                                           NodeData&      nodal_forces,
                                           const Field&   displacement) {
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

    // Element callbacks evaluate the current nonlinear stress/resultant state
    // and scatter the corresponding element internal force directly. This is a
    // residual assembly, so no material or geometric tangent matrix is built.
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

    // Contact contributes to the same nonlinear residual without constructing
    // the Gauss-Newton matrix required only by a Newton tangent evaluation.
    TripletList discarded_contact_triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(
            indices,
            *_data,
            nodal_forces,
            discarded_contact_triplets,
            false);
    }

    // Match the validation performed by tangent assembly before the nonlinear
    // solver projects the residual.
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

SparseMatrix Model::build_geom_stiffness_matrix(SystemDofIds &indices,
                                                const Field& ip_stress,
                                                const Field* stiffness_scalar)
{
    logging::error(_data->element_ip_offsets != nullptr,
                   "element IP offset field has not been initialized");
    const Field& ip_enum = *_data->element_ip_offsets;

    auto lambda = [&](const ElementPtr &e, Precision* storage) -> MapMatrix {
        if (auto sel = e->as<StructuralElement>()) {
            const ID eid = sel->elem_id;
            // start index of this element’s IPs
            const ID ip_start = static_cast<ID>(ip_enum(static_cast<Index>(eid), 0));

            MapMatrix Kg = sel->stiffness_geom(storage, ip_stress, ip_start);
            if (stiffness_scalar) {
                logging::error(stiffness_scalar->domain == FieldDomain::ELEMENT,
                               "stiffness scale field must use ELEMENT domain");
                logging::error(stiffness_scalar->components == 1,
                               "stiffness scale field must have 1 component");
                Kg *= (*stiffness_scalar)(static_cast<Index>(eid), 0);
            }
            return Kg;
        } else {
            MapMatrix mat{storage, 0, 0};
            return mat;
        }
    };

    return mattools::assemble_matrix(_data->elements, indices, lambda);
}

SparseMatrix Model::build_lumped_mass_matrix(SystemDofIds& indices) {
    auto lambda = [&](const ElementPtr &el, Precision* storage) {
        if (auto sel = el->as<StructuralElement>()) {
            MapMatrix element_mass = sel->mass(storage);
            return element_mass;
        } else {
            MapMatrix mat{storage, 0, 0};
            return mat;
        }
    };
    auto res = mattools::assemble_matrix(_data->elements, indices, lambda);
    // Add feature-based mass contributions
    if (!_data->features.empty()) {
        TripletList trips;
        for (const auto& f : _data->features) if (f) f->assemble_mass(indices, trips);
        if (!trips.empty()) res.insertFromTriplets(trips.begin(), trips.end());
    }
    return res;
}
}
}
