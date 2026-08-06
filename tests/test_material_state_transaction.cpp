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

TEST(MaterialStateTransaction, InitializesCommitsRollsBackAndRestoresBindings) {
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

    auto previous_old = std::make_shared<model::Field>(
        "PREVIOUS_OLD", model::FieldDomain::ELEMENT_MP, 4, 1);
    auto previous_new = std::make_shared<model::Field>(
        "PREVIOUS_NEW", model::FieldDomain::ELEMENT_MP, 4, 1);

    model._data->material_state_old = previous_old;
    model._data->material_state_new = previous_new;

    {
        loadcase::tools::MaterialStateTransaction transaction(model);

        ASSERT_TRUE(transaction.active());
        ASSERT_NE(model._data->material_state_old, nullptr);
        ASSERT_NE(model._data->material_state_new, nullptr);
        EXPECT_NE(model._data->material_state_old, model._data->material_state_new);

        EXPECT_EQ(model._data->material_state_old->domain, model::FieldDomain::ELEMENT_MP);
        EXPECT_EQ(model._data->material_state_old->rows, 4);
        EXPECT_EQ(model._data->material_state_old->components, 3);

        for (Index row = 0; row < 4; ++row) {
            EXPECT_EQ((*model._data->material_state_old)(row, 0), Precision(1));
            EXPECT_EQ((*model._data->material_state_old)(row, 1), Precision(2));
            EXPECT_EQ((*model._data->material_state_old)(row, 2), Precision(3));
        }

        ASSERT_NE(element->material_state_old(1, 0), nullptr);
        ASSERT_NE(element->material_state_new(1, 0), nullptr);
        EXPECT_EQ(element->material_state_old(1, 0)[1], Precision(2));

        element->material_state_new(1, 0)[1] = Precision(99);
        transaction.begin_evaluation();
        EXPECT_EQ(element->material_state_new(1, 0)[1], Precision(2));

        element->material_state_new(1, 0)[1] = Precision(7);
        transaction.commit_increment();
        EXPECT_EQ(element->material_state_old(1, 0)[1], Precision(7));

        element->material_state_new(1, 0)[1] = Precision(42);
        transaction.rollback_increment();
        EXPECT_EQ(element->material_state_old(1, 0)[1], Precision(7));
        EXPECT_EQ(element->material_state_new(1, 0)[1], Precision(7));
    }

    EXPECT_EQ(model._data->material_state_old, previous_old);
    EXPECT_EQ(model._data->material_state_new, previous_new);
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
    EXPECT_EQ(model._data->material_state_old, nullptr);
    EXPECT_EQ(model._data->material_state_new, nullptr);
    EXPECT_EQ(element->material_state_old(0, 0), nullptr);
    EXPECT_EQ(element->material_state_new(0, 0), nullptr);
}

} // namespace
