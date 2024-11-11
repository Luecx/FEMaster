/******************************************************************************
 * @file element_solid.ipp
 * @brief Implementation of the SolidElement class template. This file contains
 * the definitions for methods declared in the SolidElement class, including
 * the computation of the strain-displacement matrix, Jacobian, stiffness and
 * mass matrices, and other element-level calculations.
 *
 * @date Created on 12.06.2023
 * @author Finn Eggers
 ******************************************************************************/

#pragma once

#include "element_solid.h"

namespace fem::model {


//-----------------------------------------------------------------------------
// dofs
//-----------------------------------------------------------------------------
template<Index N>
ElDofs
SolidElement<N>::dofs() {
    return ElDofs {true, true, true, false, false, false};
}

//-----------------------------------------------------------------------------
// dimensions
//-----------------------------------------------------------------------------
template<Index N>
Dim
SolidElement<N>::dimensions() {
    return D;
}

//-----------------------------------------------------------------------------
// n_nodes
//-----------------------------------------------------------------------------
template<Index N>
Dim
SolidElement<N>::n_nodes() {
    return node_ids.size();
}

//-----------------------------------------------------------------------------
// nodes
//-----------------------------------------------------------------------------
template<Index N>
ID*
SolidElement<N>::nodes() {
    return &node_ids[0];
}

//-----------------------------------------------------------------------------
// n_integration_points
//-----------------------------------------------------------------------------
template<Index N>
Dim
SolidElement<N>::n_integration_points() {
    return integration_scheme().count();
}


}  // namespace fem::model
