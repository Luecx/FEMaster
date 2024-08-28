/******************************************************************************
 * @file EigenFrequencyAnalysis.h
 * @brief EigenFrequencyAnalysis.h defines a load case for performing modal
 * (eigenfrequency) analysis. The class computes natural frequencies and mode
 * shapes of a structure based on the finite element model.
 *
 * @details The analysis returns the first N eigenfrequencies, as defined by
 * the user, along with the corresponding mode shapes.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 27.08.2024
 *
 ******************************************************************************/

#pragma once

#include "loadcase.h"
#include "../solve/solver.h"

namespace fem {
namespace loadcase {

/******************************************************************************
 * EigenFrequencyAnalysis class
 * This class performs modal analysis on a finite element model to compute the
 * first N natural frequencies and associated mode shapes. It extends the
 * LoadCase class, implementing a run function for eigenfrequency computations.
 ******************************************************************************/
struct EigenFrequencyAnalysis : public LoadCase {

    //-------------------------------------------------------------------------
    // Constructor
    //-------------------------------------------------------------------------
    /**
     * @brief Constructs an EigenFrequencyAnalysis load case with a given ID,
     * writer, model, and number of eigenvalues.
     *
     * @param id The unique ID of the load case.
     * @param writer Pointer to a Writer object for output handling.
     * @param model Pointer to the finite element model.
     * @param numEigenvalues The number of eigenfrequencies to compute.
     */
    EigenFrequencyAnalysis(ID id, reader::Writer* writer, model::Model* model, int numEigenvalues=10)
        : LoadCase(id, writer, model), numEigenvalues(numEigenvalues) {};

    //-------------------------------------------------------------------------
    // Data Members
    //-------------------------------------------------------------------------
    int numEigenvalues;  /**< The number of eigenfrequencies to compute. */

public:

    //-------------------------------------------------------------------------
    // Run Analysis
    //-------------------------------------------------------------------------
    /**
     * @brief Executes the eigenfrequency analysis. This function overrides
     * the base run() method and computes the first `numEigenvalues` natural
     * frequencies and corresponding mode shapes.
     */
    void run() override;
};

} // namespace loadcase
} // namespace fem
