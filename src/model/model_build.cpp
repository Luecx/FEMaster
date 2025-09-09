//
// Created by Finn Eggers on 03.09.23.
//

#include "../mattools/assemble.h"
#include "../mattools/assemble_bc.h"
#include "../mattools/lump_matrix.h"
#include "../mattools/numerate_dofs.h"
#include "element/element_structural.h"
#include "model.h"
#include "solid/element_solid.h"

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

    auto res = mattools::numerate_dofs(mask);

    return res;
}

ElementData Model::build_integration_point_numeration() {
    // +1 row to store the total number of integration points
    ElementData ip_numeration{this->_data->max_elems + 1, 1};

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
    ip_numeration(this->_data->max_elems, 0) = next_id;

    return ip_numeration;
}


std::tuple<NodeData, constraint::Equations> Model::build_support_matrix(std::vector<std::string> supp_sets) {
    NodeData disp_matrix{this->_data->max_nodes, 6};
    disp_matrix.fill(std::numeric_limits<Precision>::quiet_NaN());

    // todo: add support for constraints
    constraint::Equations equations;

    for (auto &key: supp_sets) {
        auto supp_col = this->_data->supp_cols.get(key);
        supp_col->apply(*_data, disp_matrix, equations);
    }
    return {disp_matrix, equations};
}

NodeData Model::build_load_matrix(std::vector<std::string> load_sets) {
    NodeData load_matrix{this->_data->max_nodes, 6};
    load_matrix.setZero();

    for (auto &key: load_sets) {
        auto data = this->_data->load_cols.get(key);
        data->apply(*_data, load_matrix);
    }

    // apply constrained loads
    for (auto &c: this->_data->couplings) {
        c.apply_loads(*_data, load_matrix);
    }

    return load_matrix;
}

SparseMatrix Model::build_constraint_matrix   (SystemDofIds& indices, constraint::Equations& bc_equations, Precision characteristic_stiffness) {
    constraint::Equations equations;

    // add all bc equations
    for (auto &eq: bc_equations) {
        equations.push_back(eq);
    }

    for (auto &c: this->_data->couplings) {
        auto coupling_equations = c.get_equations(indices, *_data);
        for (auto &eq: coupling_equations) {
            equations.push_back(eq);
        }
    }

    for (auto &t: this->_data->ties) {
        auto tie_equations = t.get_equations(indices, *_data);
        for (auto &eq: tie_equations) {
            equations.push_back(eq);
        }
    }

    for (auto &t : this->_data->connectors) {
        auto connector_equations = t.get_equations(indices, *_data);
        for (auto &eq: connector_equations) {
            equations.push_back(eq);
        }
    }

    for (auto &eq: this->_data->equations) {
        equations.push_back(eq);
    }

    if (equations.empty()) {
        return SparseMatrix{0, indices.maxCoeff() + 1};
    }
    auto triplets = constraint::Equation::get_triplets(equations, indices, 0);

    SparseMatrix matrix{(Eigen::Index) equations.size(), indices.maxCoeff() + 1};
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix *= characteristic_stiffness;
    return matrix;
}



SparseMatrix Model::build_stiffness_matrix(SystemDofIds &indices, ElementData stiffness_scalar) {
    auto lambda = [&](const ElementPtr &el, Precision* storage) {
        if (auto sel = el->as<StructuralElement>()) {
            MapMatrix stiff = sel->stiffness(storage);
            stiff *= (stiffness_scalar.rows() > 0) ? stiffness_scalar(sel->elem_id, 0) : 1.0;
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
                                                const IPData& ip_stress,
                                                ElementData stiffness_scalar)
{
    ElementData ip_enum = build_integration_point_numeration();

    auto lambda = [&](const ElementPtr &e, Precision* storage) -> MapMatrix {
        if (auto sel = e->as<StructuralElement>()) {
            const ID eid = sel->elem_id;
            // start index of this element’s IPs
            const ID ip_start = ip_enum(eid, 0);

            MapMatrix Kg = sel->stiffness_geom(storage, const_cast<IPData&>(ip_stress), ip_start);
            Kg *= (stiffness_scalar.rows() > 0) ? stiffness_scalar(eid, 0) : 1.0;
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