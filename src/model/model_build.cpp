/**
 * @file model_build.cpp
 * @brief Implements model-level field construction and structural assembly.
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

void Model::assign_sections() {
    for (Section::Ptr ptr : _data->sections) {
        for (ID elem_id : *ptr->region_) {
            if (_data->elements[elem_id] != nullptr) {
                _data->elements[elem_id]->set_section(ptr);
            }
        }
    }
}

void Model::build_shell_element_normals(Precision equalize_angle_degrees) {
    logging::error(_data->positions_reference != nullptr,
        "Model: POSITION_REFERENCE field is not set");
    logging::error(_data->element_nodal_offsets != nullptr,
        "Model: element nodal offsets are not initialized");

    constexpr Precision normal_tolerance = Precision(1e-12);

    const Precision pi             = std::acos(Precision(-1));
    const Precision equalize_angle = equalize_angle_degrees * pi / Precision(180);
    const Precision cos_equalize   = std::cos(equalize_angle);
    const Index element_count      = static_cast<Index>(_data->elements.size());
    const Index expected_rows      = static_cast<Index>((*_data->element_nodal_offsets)(element_count));

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

    const Index node_count = _data->positions_reference->rows;
    std::vector<std::vector<Index>> node_rows(static_cast<std::size_t>(node_count));
    std::vector<Vec3>               row_normals(static_cast<std::size_t>(normals->rows), Vec3::Zero());
    std::vector<bool>               prescribed_rows(static_cast<std::size_t>(normals->rows), false);

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

            if (!added) {
                cluster_sums.push_back(normal);
                cluster_rows.push_back({row});
            }
        }

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

    _data->shell_element_nodal_normals = normals;
}

SystemDofIds Model::build_unconstrained_index_matrix() {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");

    SystemDofs mask{_data->positions->rows, 6};
    mask.fill(false);

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

    for (auto& coupling : _data->couplings) {
        const ID master_id = coupling.master_node;
        const auto master_dofs = coupling.master_dofs(mask, *_data);
        for (ID dof = 0; dof < 6; ++dof) {
            mask(master_id, dof) |= master_dofs(0, dof);
        }
    }

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
            if (feature) feature->assemble_stiffness(indices, triplets);
        }
        if (!triplets.empty()) matrix.insertFromTriplets(triplets.begin(), triplets.end());
    }

    return matrix;
}

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

    Field ip_stress_state{
        "IP_STRESS_STATE", FieldDomain::ELEMENT_IP,
        _data->field_rows(FieldDomain::ELEMENT_IP), 8};
    ip_stress_state.set_zero();

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

    for (Index i = 0; i < nodal_forces.rows; ++i) {
        for (Index j = 0; j < nodal_forces.components; ++j) {
            const bool bad = std::isnan(nodal_forces(i, j)) || std::isinf(nodal_forces(i, j));
            logging::error(!bad,
                "Internal force at node ", i, " has invalid value at col ", j);
        }
    }

    return global_matrix;
}

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

    Field ip_stress_state{
        "IP_STRESS_STATE", FieldDomain::ELEMENT_IP,
        _data->field_rows(FieldDomain::ELEMENT_IP), 8};
    ip_stress_state.set_zero();
    nodal_forces.set_zero();

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

    TripletList discarded_contact_triplets;
    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, discarded_contact_triplets);
    }

    for (Index i = 0; i < nodal_forces.rows; ++i) {
        for (Index j = 0; j < nodal_forces.components; ++j) {
            const bool bad = std::isnan(nodal_forces(i, j)) || std::isinf(nodal_forces(i, j));
            logging::error(!bad,
                "Internal force at node ", i, " has invalid value at col ", j);
        }
    }
}

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

    return mattools::assemble_matrix(_data->elements, indices, lambda);
}

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
            if (feature) feature->assemble_mass(indices, triplets);
        }
        if (!triplets.empty()) matrix.insertFromTriplets(triplets.begin(), triplets.end());
    }

    return matrix;
}

} // namespace model
} // namespace fem
