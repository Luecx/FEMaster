#include "../src/loadcase/tools/material_state_transaction.h"
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
    explicit TestElement(ID id)
        : ElementInterface(id) {}

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

    std::array<ID, 1> node_ids{0};
};

TEST(MaterialStateTransaction, InitializesCommitsRollsBackAndRestoresBinding) {
    model::Model model(1, 1, 0);

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
        loadcase::tools::MaterialStateTransaction transaction(model);

        ASSERT_TRUE(transaction.active());
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
        Precision* state = model._data->material_state->data()
                         + row * model._data->material_state->components;
        EXPECT_EQ(state[1], Precision(2));

        state[1] = Precision(99);
        transaction.begin_evaluation();
        EXPECT_EQ((*model._data->material_state)(row, 1), Precision(2));

        (*model._data->material_state)(row, 1) = Precision(7);
        transaction.commit_increment();
        EXPECT_EQ((*model._data->material_state)(row, 1), Precision(7));

        (*model._data->material_state)(row, 1) = Precision(42);
        transaction.rollback_increment();
        EXPECT_EQ((*model._data->material_state)(row, 1), Precision(7));
    }

    EXPECT_EQ(model._data->material_state, previous_state);
}

TEST(MaterialStateTransaction, RemainsInactiveForStatelessMaterials) {
    model::Model model(1, 1, 0);

    auto material = std::make_shared<material::Material>("STATELESS");

    auto section = std::make_shared<TestSection>();
    section->material_ = material;

    auto element = std::make_shared<TestElement>(0);
    element->set_section(section);
    element->_model_data = model._data.get();
    model._data->elements[0] = element;
    model._data->initialize_element_enumeration();

    loadcase::tools::MaterialStateTransaction transaction(model);

    EXPECT_FALSE(transaction.active());
    EXPECT_EQ(model._data->material_state, nullptr);
}

} // namespace
