//
// Created by Finn Eggers on 03.09.23.
//

#include "../mattools/assemble.h"
#include "../mattools/numerate_dofs.h"
#include "element/element_structural.h"
#include "model.h"
#include "solid/element_solid.h"

#include <iterator>
#include <utility>
#include <string>
#include <cmath>

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

            load->amplitude_ = nullptr;
            load->apply(*_data, load_matrix, Precision(0));
            load->amplitude_ = amplitude;
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
