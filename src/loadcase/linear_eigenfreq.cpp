#include "linear_eigenfreq.h"
#include "../solve/eigen.h"

// fem::loadcase::EigenFrequencyAnalysis::EigenFrequencyAnalysis(ID id, reader::Writer* writer, model::Model* model, int numEigenvalues)
//     : LoadCase(id, writer, model), numEigenvalues(numEigenvalues) {}

void fem::loadcase::EigenFrequencyAnalysis::run() {
    // logging::info(true, "");
    // logging::info(true, "");
    // logging::info(true, "================================================================================================");
    // logging::info(true, "EIGENFREQUENCY ANALYSIS");
    // logging::info(true, "================================================================================================");
    // logging::info(true, "");
    //
    // auto unconstrained = Timer::measure(
    //     [&]() { return this->m_model->build_unconstrained_index_matrix(); },
    //     "generating unconstrained index matrix"
    // );
    //
    // auto supp_vec = Timer::measure(
    //     [&]() { return this->m_model->build_support_vector(unconstrained); },
    //     "building support vector"
    // );
    //
    // auto stiffness = Timer::measure(
    //     [&]() { return this->m_model->build_stiffness_matrix(unconstrained); },
    //     "constructing stiffness matrix"
    // );
    //
    // auto mass = Timer::measure(
    //     [&]() { return this->m_model->build_mass_matrix(unconstrained); },
    //     "constructing mass matrix"
    // );
    //
    // auto constrained = Timer::measure(
    //     [&]() { return this->m_model->build_constrained_index_matrix(supp_vec); },
    //     "creating constrained index matrix"
    // );
    //
    // auto mapping_vec = Timer::measure(
    //     [&]() { return this->m_model->build_mapping_vector(unconstrained, constrained); },
    //     "mapping unconstrained indices to constrained ones"
    // );
    //
    // auto reduced_stiffness = Timer::measure(
    //     [&]() { return this->m_model->build_reduced_stiffness(mapping_vec, stiffness); },
    //     "reducing stiffness matrix"
    // );
    //
    // // auto reduced_mass = Timer::measure(
    // //     [&]() { return this->m_model->build_reduced_mass(mapping_vec, mass); },
    // //     "reducing mass matrix"
    // // );
    //
    // auto [eigenvalues, eigenvectors] = Timer::measure(
    //     [&]() { return solver::compute_eigenvalues(device, reduced_stiffness, numEigenvalues, true); },
    //     "computing eigenvalues and eigenvectors"
    // );
    //
    // // Optionally, you can log or display the eigenvalues here
    // logging::info(true, "Computed eigenvalues:");
    // for (size_t i = 0; i < eigenvalues.size(); ++i) {
    //     logging::info(true, "Eigenvalue ", i + 1, ": ", eigenvalues[i]);
    // }

    // Note: No results are written out to files as requested
}
