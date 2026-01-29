//
// Created by Finn Eggers on 03.09.23.
//

#include "../mattools/assemble.h"
#include "../mattools/lump_matrix.h"
#include "../mattools/numerate_dofs.h"
#include "element/element_structural.h"
#include "model.h"
#include "solid/element_solid.h"

#include <iterator>
#include <utility>
#include <string>

namespace fem {

namespace model {

void Model::assign_sections() {
    for (Section::Ptr ptr: this->_data->sections) {
        for (ID elem_id: *ptr->region) {
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

Field Model::build_integration_point_numeration() {
    // +1 row to store the total number of integration points
    Field ip_numeration{"IP_ENUM", FieldDomain::ELEMENT,
                        static_cast<Index>(this->_data->max_elems + 1), 1};

    ID next_id = 0;
    for (auto &e : _data->elements) {
        if (e != nullptr) {
            if (auto el = e->as<StructuralElement>()) {
                ip_numeration(e->elem_id, 0) = next_id;
                next_id += el->n_integration_points();
            }
        }
    }

    // last entry = total number of IPs
    ip_numeration(static_cast<Index>(this->_data->max_elems), 0) = next_id;

    return ip_numeration;
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
    return res;
}

SparseMatrix Model::build_geom_stiffness_matrix(SystemDofIds &indices,
                                                const Field& ip_stress,
                                                const Field* stiffness_scalar)
{
    Field ip_enum = build_integration_point_numeration();

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
    return res;
}


}
}
