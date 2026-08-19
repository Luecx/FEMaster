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

TEST(NonlinearStateManager, ResetsCommitsAndRestoresMaterialStateBinding) {
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

    auto previous_state = std::make_shared<model::Field>(
        "PREVIOUS_STATE", model::FieldDomain::ELEMENT_MP, 4, 1);
    model._data->material_state = previous_state;

    {
        loadcase::tools::NonlinearStateManager state(model);

        ASSERT_NE(model._data->material_state, nullptr);
        EXPECT_EQ(model._data->material_state->domain, model::FieldDomain::ELEMENT_MP);
        EXPECT_EQ(model._data->material_state->rows, 4);
        EXPECT_EQ(model._data->material_state->components, 3);

        for (Index row = 0; row < 4; ++row) {
            EXPECT_EQ((*model._data->material_state)(row, 0), Precision(1));
            EXPECT_EQ((*model._data->material_state)(row, 1), Precision(2));
            EXPECT_EQ((*model._data->material_state)(row, 2), Precision(3));
        }

        const Index row = element->mp_index(1, 0);
        Precision* material_state = &(*model._data->material_state)(row, 0);
        EXPECT_EQ(material_state[1], Precision(2));

        material_state[1] = Precision(99);
        state.reset_material_state();
        EXPECT_EQ((*model._data->material_state)(row, 1), Precision(2));

        (*model._data->material_state)(row, 1) = Precision(7);
        state.commit_material_state();
        EXPECT_EQ((*model._data->material_state)(row, 1), Precision(7));

        (*model._data->material_state)(row, 1) = Precision(42);
        state.reset_material_state();
        EXPECT_EQ((*model._data->material_state)(row, 1), Precision(7));
    }

    EXPECT_EQ(model._data->material_state, previous_state);
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

    auto enumerated_state = model._data->material_state;
    ASSERT_NE(enumerated_state, nullptr);

    {
        loadcase::tools::NonlinearStateManager state(model);
        EXPECT_EQ(model._data->material_state, enumerated_state);
    }

    EXPECT_EQ(model._data->material_state, enumerated_state);
}

} // namespace
