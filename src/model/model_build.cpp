//
// Created by Finn Eggers on 03.09.23.
//

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
 * Builds an element-nodal reference normal field for shell elements.
 *
 * Each shell element first contributes the physical reference normal evaluated
 * at its natural nodes. Contributions that meet at the same global node are
 * then grouped by angular compatibility. Normals inside one group are
 * equalised by averaging, while groups separated by more than the supplied
 * angle remain element-local. This preserves sharp shell folds but provides a
 * smooth director field on gently curved shell patches.
 *
 * @param equalize_angle_degrees Maximum angle for equalising adjacent shell
 *                               normals at one global node.
 */
void Model::build_shell_element_normals(Precision equalize_angle_degrees) {
    logging::error(_data->positions_reference != nullptr,
                   "Model: POSITION_REFERENCE field is not set");
    logging::error(_data->element_nodal_offsets != nullptr,
                   "Model: element nodal offsets are not initialized");

    const Precision pi             = std::acos(Precision(-1));
    const Precision equalize_angle = equalize_angle_degrees * pi / Precision(180);
    const Precision cos_equalize   = std::cos(equalize_angle);

    auto normals = _data->create_field(
        "SHELL_ELEMENT_NODAL_NORMALS", FieldDomain::ELEMENT_NODAL, 3, false);
    normals->set_zero();

    std::vector<std::vector<Index>> node_rows(static_cast<std::size_t>(_data->max_nodes));
    std::vector<Vec3>               row_normals(static_cast<std::size_t>(normals->rows), Vec3::Zero());

    // Evaluate the raw shell normal owned by each element-node pair.
    for (const ElementPtr& element: _data->elements) {
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
        logging::error(static_cast<Index>(natural_coords.rows()) == n_nodes &&
                       static_cast<Index>(natural_coords.cols()) == 2,
                       "Model: natural surface node coordinates do not match shell element ",
                       structural->elem_id);

        for (Index local_node = 0; local_node < n_nodes; ++local_node) {
            const ID    node_id = structural->nodes()[local_node];
            const Index row     = offset + local_node;
            const Vec2  local   = natural_coords.row(local_node).transpose();
            Vec3        normal  = surface->normal(*_data->positions_reference,
                                                  local);

            logging::error(normal.allFinite() && normal.norm() > Precision(0),
                           "Model: invalid shell normal in element ", structural->elem_id);

            normal.normalize();

            row_normals[static_cast<std::size_t>(row)] = normal;
            node_rows[static_cast<std::size_t>(node_id)].push_back(row);
        }
    }

    // Equalise compatible normals at each global node. A node may produce
    // several clusters when it lies on a fold sharper than the threshold.
    for (const auto& rows: node_rows) {
        std::vector<Vec3>          cluster_normals;
        std::vector<Precision>     cluster_weights;
        std::vector<std::vector<Index>> cluster_rows;

        for (Index row: rows) {
            const Vec3 normal = row_normals[static_cast<std::size_t>(row)];
            bool       added  = false;

            for (Index cluster = 0; cluster < static_cast<Index>(cluster_normals.size()); ++cluster) {
                if (normal.dot(cluster_normals[static_cast<std::size_t>(cluster)]) < cos_equalize) {
                    continue;
                }

                Vec3&     cluster_normal = cluster_normals[static_cast<std::size_t>(cluster)];
                Precision& weight        = cluster_weights[static_cast<std::size_t>(cluster)];

                cluster_normal = (weight * cluster_normal + normal).normalized();
                weight += Precision(1);
                cluster_rows[static_cast<std::size_t>(cluster)].push_back(row);
                added = true;
                break;
            }

            if (!added) {
                cluster_normals.push_back(normal);
                cluster_weights.push_back(Precision(1));
                cluster_rows.push_back({row});
            }
        }

        for (Index cluster = 0; cluster < static_cast<Index>(cluster_normals.size()); ++cluster) {
            const Vec3 normal = cluster_normals[static_cast<std::size_t>(cluster)];
            for (Index row: cluster_rows[static_cast<std::size_t>(cluster)]) {
                (*normals)(row, 0) = normal(0);
                (*normals)(row, 1) = normal(1);
                (*normals)(row, 2) = normal(2);
            }
        }
    }

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

    nodal_forces.set_zero();

    const int   global_size = indices.maxCoeff() + 1;
    const Index max_ip      = this->_data->max_integration_points;
    TripletList triplets;
    SparseMatrix global_matrix(global_size, global_size);

    // create a stress state which is optional to be used by the elements and can hold up
    // to 8 entries
    Field ip_stress_state{ "IP_STRESS_STATE", FieldDomain::ELEMENT_IP, max_ip, 8};
    ip_stress_state.set_zero();

    for (const auto& element : _data->elements) {
        if (!element) {
            continue;
        }

        auto* structural = element->as<StructuralElement>();
        if (!structural) {
            continue;
        }

        constexpr int MAX_LOCAL_MATRIX_SIZE = 128;
        alignas(64) Precision local_matrix_storage[MAX_LOCAL_MATRIX_SIZE * MAX_LOCAL_MATRIX_SIZE]{};

        MapMatrix tangent = structural->stiffness_tangent(
            local_matrix_storage,
            ip_stress_state,
            nodal_forces,
            displacement
        );

        if (stiffness_scalar) {
            logging::error(stiffness_scalar->domain == FieldDomain::ELEMENT,
                           "stiffness scale field must use ELEMENT domain");
            logging::error(stiffness_scalar->components == 1,
                           "stiffness scale field must have 1 component");
            tangent *= (*stiffness_scalar)(static_cast<Index>(structural->elem_id), 0);
        }

        const int num_nodes = element->n_nodes();
        const int local_matrix_size = static_cast<int>(tangent.rows());
        const int dofs_per_node = local_matrix_size / num_nodes;

        for (int i = 0; i < num_nodes; ++i) {
            for (int j = 0; j < num_nodes; ++j) {
                for (int idof = 0; idof < dofs_per_node; ++idof) {
                    for (int jdof = 0; jdof < dofs_per_node; ++jdof) {
                        const int global_row = indices(element->nodes()[i], idof);
                        const int global_col = indices(element->nodes()[j], jdof);

                        if (global_row < 0 || global_col < 0) {
                            continue;
                        }

                        triplets.emplace_back(
                            global_row,
                            global_col,
                            tangent(i * dofs_per_node + idof, j * dofs_per_node + jdof)
                        );
                    }
                }
            }
        }
    }

    for (const auto& contact : _data->contacts) {
        contact.assemble(indices, *_data, nodal_forces, triplets);
    }

    global_matrix.insertFromTriplets(triplets.begin(), triplets.end());

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
            logging::error(!bad, "Internal force at node ", i, " has invalid value at col ", j);
        }
    }

    return global_matrix;
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
