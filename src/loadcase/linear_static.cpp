//
// Created by Luecx on 04.09.2023.
//
#include "linear_static.h"

fem::loadcase::LinearStatic::LinearStatic(fem::model::Model* model)
    : model(model) {}

void fem::loadcase::LinearStatic::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR STATIC");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

    auto unconstrained = Timer::measure(
        [&]() { return model->build_unconstrained_index_matrix(); },
        "generating unconstrained index matrix"
    );

    auto supp_vec = Timer::measure(
        [&]() { return model->build_support_vector(unconstrained); },
        "building support vector"
    );

    auto load_vec = Timer::measure(
        [&]() { return model->build_load_vector(unconstrained); },
        "formulating load vector"
    );

    auto stiffness = Timer::measure(
        [&]() { return model->build_stiffness_matrix(unconstrained); },
        "constructing stiffness matrix"
    );

    auto impl_load_vec = Timer::measure(
        [&]() { return model->build_implicit_load_vector(stiffness, supp_vec); },
        "computing implicit load vector"
    );

    auto constrained = Timer::measure(
        [&]() { return model->build_constrained_index_matrix(supp_vec); },
        "creating constrained index matrix"
    );

    auto mapping_vec = Timer::measure(
        [&]() { return model->build_mapping_vector(unconstrained, constrained); },
        "mapping unconstrained indices to constrained ones"
    );

    load_vec += impl_load_vec;

    auto reduced_stiffness = Timer::measure(
        [&]() { return model->build_reduced_stiffness(mapping_vec, stiffness); },
        "reducing stiffness matrix"
    );

    auto reduced_load = Timer::measure(
        [&]() { return model->build_reduced_load(mapping_vec, load_vec); },
        "reducing load vector"
    );

    auto displacement = Timer::measure(
        [&]() { return solver::solve(device, method, reduced_stiffness, reduced_load); },
        "solving system of equations"
    );
}
