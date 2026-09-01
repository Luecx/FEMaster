/**
 * @file test_nonlinear_state_manager.cpp
 * @brief Tests nonlinear material-state ownership and lifecycle transitions.
 */

#include "../src/loadcase/tools/nonlinear_state_manager.h"
#include "../src/material/elasticity.h"
#include "../src/model/model.h"
#include "../src/section/section.h"

#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <string>

namespace {

using namespace fem;

struct StatefulElasticity final : material::Elasticity {
    Index state_size() const override {
        return 3;
    }

    void initialize_state(Precision* state) const override {
        state[0] = Precision(1);
        state[1] = Precision(2);
        state[2] = Precision(3);
    }
};

struct TestSection final : Section {
    void info() override {}

    std::string str() const override {
        return "TestSection";
    }
};

struct TestElement final : model::ElementInterface {
    std::array<ID, 1> node_ids{0};

    explicit TestElement(ID id)
        : ElementInterface(id) {}

    model::ElementPtr copy() const override {
        auto copy = std::make_shared<TestElement>(elem_id);
        copy->node_ids = node_ids;
        return copy;
    }

    ElDofs dofs() const override {
        return ElDofs{true, false, false, false, false, false};
    }

    Dim dimensions() const override {
        return 1;
    }

    Dim n_nodes() const override {
        return 1;
    }

    Dim num_ip() const override {
        return 2;
    }

    Index num_mp_per_ip() const override {
        return 2;
    }

    const ID* nodes() const override {
        return node_ids.data();
    }
};

TEST(NonlinearStateManager, SeparatesResetsAndCommitsMaterialState) {
    model::Model model;
    model._data = std::make_shared<model::ModelData>();
    model._data->elements.resize(1);

    auto material = std::make_shared<material::Material>("STATEFUL");
    material->set_elasticity<StatefulElasticity>();

    auto section = std::make_shared<TestSection>();
    section->material_ = material;

    auto element = std::make_shared<TestElement>(0);
    element->set_section(section);
    element->_model_data = model._data.get();
    model._data->elements[0] = element;
    model._data->initialize_element_enumeration();

    auto enumerated_state_old = model._data->material_state_old;
    auto enumerated_state_new = model._data->material_state_new;

    model::Field::Ptr committed_state;
    model::Field::Ptr trial_state;

    {
        loadcase::tools::NonlinearStateManager state(model);

        ASSERT_NE(model._data->material_state_old, nullptr);
        ASSERT_NE(model._data->material_state_new, nullptr);
        EXPECT_NE(model._data->material_state_old, model._data->material_state_new);
        EXPECT_EQ(model._data->material_state_old, enumerated_state_old);
        EXPECT_EQ(model._data->material_state_new, enumerated_state_new);
        EXPECT_EQ(model._data->material_state_old->domain, model::FieldDomain::ELEMENT_MP);
        EXPECT_EQ(model._data->material_state_old->rows, 4);
        EXPECT_EQ(model._data->material_state_old->components, 3);

        for (Index row = 0; row < 4; ++row) {
            EXPECT_EQ((*model._data->material_state_old)(row, 0), Precision(1));
            EXPECT_EQ((*model._data->material_state_old)(row, 1), Precision(2));
            EXPECT_EQ((*model._data->material_state_old)(row, 2), Precision(3));
            EXPECT_EQ((*model._data->material_state_new)(row, 0), Precision(1));
            EXPECT_EQ((*model._data->material_state_new)(row, 1), Precision(2));
            EXPECT_EQ((*model._data->material_state_new)(row, 2), Precision(3));
        }

        const Index row = element->mp_index(1, 0);
        const Precision* old_state = &(*model._data->material_state_old)(row, 0);
        Precision*       new_state = &(*model._data->material_state_new)(row, 0);
        EXPECT_EQ(old_state[1], Precision(2));
        EXPECT_EQ(new_state[1], Precision(2));

        new_state[1] = Precision(99);
        EXPECT_EQ(old_state[1], Precision(2));
        state.reset_material_state();
        EXPECT_EQ((*model._data->material_state_new)(row, 1), Precision(2));

        (*model._data->material_state_new)(row, 1) = Precision(7);
        state.commit_material_state();
        EXPECT_EQ((*model._data->material_state_old)(row, 1), Precision(7));
        EXPECT_EQ((*model._data->material_state_new)(row, 1), Precision(7));

        (*model._data->material_state_new)(row, 1) = Precision(42);
        state.reset_material_state();
        EXPECT_EQ((*model._data->material_state_old)(row, 1), Precision(7));
        EXPECT_EQ((*model._data->material_state_new)(row, 1), Precision(7));

        committed_state = model._data->material_state_old;
        trial_state     = model._data->material_state_new;
    }

    EXPECT_EQ(model._data->material_state_old, committed_state);
    EXPECT_EQ(model._data->material_state_new, trial_state);
}

TEST(NonlinearStateManager, KeepsEnumeratedStateForStatelessMaterials) {
    model::Model model;
    model._data = std::make_shared<model::ModelData>();
    model._data->elements.resize(1);

    auto material = std::make_shared<material::Material>("STATELESS");

    auto section = std::make_shared<TestSection>();
    section->material_ = material;

    auto element = std::make_shared<TestElement>(0);
    element->set_section(section);
    element->_model_data = model._data.get();
    model._data->elements[0] = element;
    model._data->initialize_element_enumeration();

    auto enumerated_state_old = model._data->material_state_old;
    auto enumerated_state_new = model._data->material_state_new;
    ASSERT_NE(enumerated_state_old, nullptr);
    ASSERT_NE(enumerated_state_new, nullptr);

    {
        loadcase::tools::NonlinearStateManager state(model);
        EXPECT_EQ(model._data->material_state_old, enumerated_state_old);
        EXPECT_EQ(model._data->material_state_new, enumerated_state_new);
    }

    EXPECT_EQ(model._data->material_state_old, enumerated_state_old);
    EXPECT_EQ(model._data->material_state_new, enumerated_state_new);
}

} // namespace
