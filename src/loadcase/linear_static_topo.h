/**
 * @file LinearStaticTopo.h
 * @brief LinearStaticTopo.h defines a load case that extends linear static 
 * analysis to include topology optimization parameters such as element density 
 * and penalization exponent. This class runs the analysis considering 
 * these optimization effects.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 27.08.2024
 *
 */

#pragma once

#include "loadcase.h"
#include "linear_static.h"
#include "../solve/solver.h"

namespace fem {
namespace loadcase {

/**
 * LinearStaticTopo class
 * This class represents a load case for linear static analysis with topology
 * optimization considerations. It extends the LinearStatic class, adding
 * element density and a penalization exponent to account for _material
 * interpolation during topology optimization.
 */
struct LinearStaticTopo : public LinearStatic {

    //-------------------------------------------------------------------------
    // Constructor
    //-------------------------------------------------------------------------
    /**
     * @brief Constructs a LinearStaticTopo load case with a given ID, writer, 
     * and model.
     *
     * @param id The unique ID of the load case.
     * @param writer Pointer to a Writer object for output handling.
     * @param model Pointer to the finite element model.
     */
    explicit LinearStaticTopo(ID id, reader::Writer* writer, model::Model* model);

    //-------------------------------------------------------------------------
    // Data Members
    //-------------------------------------------------------------------------
    ElementData density;  /**< Element density used for topology optimization. */
    ElementData orientation;  /**< Element orientation used for topology optimization. */
    Precision exponent = 1;  /**< Exponent used in SIMP-like penalization schemes. */

public:

    //-------------------------------------------------------------------------
    // Run Analysis
    //-------------------------------------------------------------------------
    /**
     * @brief Executes the linear static analysis with topology optimization 
     * considerations. This function overrides the base run() method.
     */
    void run() override;
};

} // namespace loadcase
} // namespace fem
