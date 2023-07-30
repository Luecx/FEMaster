
#include "model.h"
namespace fem {

namespace model {

IndexMatrix Model::build_dof_index_vector(DynamicVector& constraints){
    // first build a boolean matrix of which dofs are present in the system
    BooleanMatrix mask{this->max_nodes, 6};

    // go through each elements and mask the dofs for all the nodes
    for(auto& e:elements){
        if(e != nullptr){
            auto dofs = e->dofs();
            for(ID node_local_id = 0; node_local_id < e->n_nodes(); node_local_id++){
                ID node_id = e->nodes()[node_local_id];
                for(ID dof = 0; dof < 6; dof++){
                    mask(node_id, dof) |= dofs(dof);
                }
            }
        }
    }

    // numerate the entries in the mask
    IndexMatrix res{this->max_nodes, 6};

    ID c = 0;
    for(ID n = 0; n < this->max_nodes; n++){
        for(ID d = 0; d < 6; d++){
            if(mask(n,d)){
                res(n,d) = c++;
            }else{
                res(n,d) = NAN;
            }
        }
    }

    return res;
}

DynamicVector Model::solve(DynamicVector& constraints, DynamicVector& loads) {
    //        SparseMatrixBuilder tripplets{};
}

SparseMatrix Model::stiffness(){

}

void Model::apply_constraints(DynamicVector& constraints, SparseMatrix& matrix, DynamicVector& loads){

}

}    // namespace model
}    // namespace fem
