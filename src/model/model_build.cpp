/**
 * @file model_build.cpp
 * @brief Implements model-level field construction and structural assembly.
 *
 * The routines operate on the dense assembly created by `Model::compile()`.
 * They bind sections to compiled elements, construct shared shell reference
 * normals, enumerate active generalized DOFs, assemble loads and constraint
 * equations, and combine regular elements, post-compile point elements, contact
 * and generic feature contributions into global sparse operators and nodal
 * internal-force fields.
 *
 * FEMaster `POINTMASS` commands create auxiliary `PointElement` objects after
 * compilation because they target compiled NSETs. Those point elements are kept
 * outside dense ELSET and element-field enumeration, but participate in the same
 * DOF, stiffness and mass operators as regular structural elements.
 *
 * Element-local matrix evaluation remains the responsibility of each structural
 * formulation. `mattools::assemble_matrix()` owns local-to-global DOF mapping,
 * sparse triplet accumulation and parallel reduction. Stress and stress
 * resultants needed for geometric or nonlinear operators remain element-local;
 * model assembly exchanges only displacement, matrices and nodal forces.
 *
 * @see Model
 * @see Model::compile
 * @see mattools::assemble_matrix
 *
 * @author Finn Eggers
 * @date 25.08.2026
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
    const Index expected_rows      = static_cast<Index>((*_data->element_nodal_offsets)(element_count));

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
 * their formulations. Generic non-element features, couplings and connectors may
 * additionally activate generalized nodal directions. The final boolean mask is
 * converted to contiguous zero-based system indices.
 *
 * @return Node-by-six matrix of global unconstrained system DOF identifiers.
 */
SystemDofIds Model::build_unconstrained_index_matrix() {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");

    SystemDofs mask{_data->positions->rows, 6};
    mask.fill(false);

    // Activate generalized DOFs used by regular compiled elements
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

    // Native POINTMASS objects are point elements created after compilation.
    for (auto& element : _data->point_elements) {
        if (element == nullptr) continue;

        const auto dofs = element->dofs();
        for (ID local_node = 0; local_node < element->n_nodes(); ++local_node) {
            const ID node_id = element->nodes()[local_node];
            for (ID dof = 0; dof < 6; ++dof) {
                mask(node_id, dof) |= dofs(0, dof);
            }
        }
    }

    // Activate generalized components carrying other non-element feature data
    for (const auto& feature : _data->features) {
        if (feature) feature->activate_dofs(mask);
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
    Field load_matrix{"LOAD_MATRIX", FieldDomain::NODE, _data->field_rows(FieldDomain::NODE), 6};
    load_matrix.set_zero();

    for (auto& key : load_sets) {
        auto data = _data->load_cols.get(key);
        data->apply(*_data, load_matrix, time);
    }

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

    auto field_for = [this, &basis](const bc::Amplitude::Ptr& amplitude) -> Field& {
        for (auto& entry : basis) {
            if (entry.first == amplitude) return entry.second;
        }

        Field load_matrix{"LOAD_BASIS", FieldDomain::NODE, _data->field_rows(FieldDomain::NODE), 6};
        load_matrix.set_zero();
        basis.emplace_back(amplitude, std::move(load_matrix));
        return basis.back().second;
    };

    for (auto& key : load_sets) {
        auto data = _data->load_cols.get(key);
        for (const auto& load : data->entries()) {
            if (!load) continue;
            auto& load_matrix = field_for(load->amplitude_);
            load->apply(*_data, load_matrix, Precision(0), true);
        }
    }

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
 * @param system_dof_ids Active global DOF identifiers used by constraint builders.
 * @param supp_sets Optional names of support collectors to include.
 * @return Constraint equations grouped and annotated by their origin.
 */
constraint::ConstraintGroups Model::collect_constraints(
    SystemDofIds&                   system_dof_ids,
    const std::vector<std::string>& supp_sets
) {
    constraint::ConstraintGroups groups{};

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
 * Regular compiled elements and post-compile native point elements both use the
 * common sparse element assembly path. The optional ELEMENT-domain stiffness
 * scale applies only to dense compiled elements because auxiliary point elements
 * have no element-domain row.
 *
 * @param indices Active global DOF identifiers used for local-to-global mapping.
 * @param stiffness_scalar Optional one-component element stiffness scale.
 * @return Assembled sparse linear stiffness matrix.
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

    if (!_data->point_elements.empty()) {
        auto point_lambda = [](const ElementPtr& element, Precision* storage) {
            if (auto structural = element->as<StructuralElement>()) {
                return structural->stiffness(storage);
            }
            MapMatrix matrix{storage, 0, 0};
            return matrix;
        };
        matrix += mattools::assemble_matrix(_data->point_elements, indices, point_lambda);
    }

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
 * Every structural element performs one physical trial evaluation through
 * `stiffness_tangent()`. The same constitutive response supplies its internal
 * force and, when matrix storage is present, its consistent tangent. No global
 * integration-point stress scratch field is allocated during assembly.
 *
 * Native POINTMASS elements use the same path: their linear spring law produces
 * the matching nodal force and constant tangent directly.
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
    logging::error(nodal_forces.domain == FieldDomain::NODE,
        "tangent internal force output must use NODE domain");
    logging::error(nodal_forces.rows == _data->field_rows(FieldDomain::NODE),
        "tangent internal force output has wrong node count");
    logging::error(nodal_forces.components >= 6,
        "tangent internal force output requires at least 6 components");

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

    if (!_data->point_elements.empty()) {
        auto point_lambda = [&](const ElementPtr& element,
                                Precision*        storage,
                                NodeData&         local_nodal_forces) -> MapMatrix {
            if (auto structural = element->as<StructuralElement>()) {
                return structural->stiffness_tangent(storage, local_nodal_forces, displacement);
            }
            MapMatrix matrix{storage, 0, 0};
            return matrix;
        };
        global_matrix += mattools::assemble_matrix(
            _data->point_elements,
            indices,
            point_lambda,
            &nodal_forces
        );
    }

    TripletList contact_triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, contact_triplets);
    }

    if (!contact_triplets.empty()) {
        SparseMatrix contact_matrix(global_matrix.rows(), global_matrix.cols());
        contact_matrix.insertFromTriplets(contact_triplets.begin(), contact_triplets.end());
        global_matrix += contact_matrix;
    }

    if (!_data->features.empty()) {
        TripletList feature_triplets;
        for (const auto& feature : _data->features) {
            if (feature) feature->assemble_stiffness(indices, feature_triplets);
        }
        if (!feature_triplets.empty()) {
            global_matrix.insertFromTriplets(feature_triplets.begin(), feature_triplets.end());
        }
    }

    nodal_forces.check_finite("Internal force");
    return global_matrix;
}

/**
 * Evaluates nonlinear internal and contact forces without retaining a tangent.
 *
 * Structural elements use the same physical `stiffness_tangent()` path as full
 * Newton assembly, but receive a null matrix buffer. They therefore update trial
 * material state and assemble their matching internal force without constructing
 * tangent-only terms. Contact follows with tangent triplets discarded.
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
    logging::error(nodal_forces.domain == FieldDomain::NODE,
        "nonlinear internal force output must use NODE domain");
    logging::error(nodal_forces.rows == _data->field_rows(FieldDomain::NODE),
        "nonlinear internal force output has wrong node count");
    logging::error(nodal_forces.components >= 6,
        "nonlinear internal force output requires at least 6 components");

    nodal_forces.set_zero();

    for (const auto& element : _data->elements) {
        if (!element) continue;
        auto* structural = element->as<StructuralElement>();
        if (!structural) continue;
        structural->stiffness_tangent(nullptr, nodal_forces, displacement);
    }

    for (const auto& element : _data->point_elements) {
        if (!element) continue;
        auto* structural = element->as<StructuralElement>();
        if (!structural) continue;
        structural->stiffness_tangent(nullptr, nodal_forces, displacement);
    }

    TripletList discarded_contact_triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, discarded_contact_triplets);
    }

    nodal_forces.check_finite("Internal force");
}

/**
 * Assembles the global geometric stiffness from a supplied prestress displacement.
 *
 * Each dense structural element evaluates the stress or stress resultants needed
 * by its geometric operator locally from `displacement` and committed material
 * history. No global integration-point stress field is created or consumed.
 * Auxiliary point elements are excluded because they carry no stress-dependent
 * geometric stiffness.
 *
 * @param indices Active global DOF identifiers used for sparse assembly.
 * @param displacement Global nodal displacement field defining the prestress.
 * @param stiffness_scalar Optional one-component element stiffness scale.
 * @return Assembled sparse geometric stiffness matrix.
 */
SparseMatrix Model::build_geom_stiffness_matrix(
    SystemDofIds& indices,
    const Field&  displacement,
    const Field*  stiffness_scalar
) {
    auto lambda = [&](const ElementPtr& element, Precision* storage) -> MapMatrix {
        if (auto structural = element->as<StructuralElement>()) {
            MapMatrix geometric_stiffness = structural->stiffness_geom(storage, displacement);

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
 * assembly path. Other generic features may still append independent mass
 * triplets afterwards.
 *
 * @param indices Active global DOF identifiers used for sparse assembly.
 * @return Assembled sparse lumped mass matrix.
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

    if (!_data->point_elements.empty()) {
        matrix += mattools::assemble_matrix(_data->point_elements, indices, lambda);
    }

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
