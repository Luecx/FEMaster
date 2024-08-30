#include <gtest/gtest.h>
#include "../src/model/element_solid.h"
#include "../src/model/c3d4.h"
#include "../src/model/c3d6.h"
#include "../src/model/c3d8.h"
#include "../src/model/c3d10.h"
#include "../src/model/c3d15.h"
#include "../src/model/c3d20.h"

using namespace fem::model;

// Test class template for elements with specific number of nodes
template <typename ElementType, size_t N>
class ElementTest : public ::testing::Test {
    protected:
    ElementType* element;

    void SetUp() override {
        std::array<ID, N> node_ids;
        for (size_t i = 0; i < N; ++i) {
            node_ids[i] = i; // Example node IDs
        }
        element = new ElementType(1, node_ids);
    }

    void TearDown() override {
        delete element;
    }
};

// Test cases for each specific element type and number of nodes
class C3D4Test : public ElementTest<C3D4, 4> {};
class C3D6Test : public ElementTest<C3D6, 6> {};
class C3D8Test : public ElementTest<C3D8, 8> {};
class C3D10Test : public ElementTest<C3D10, 10> {};
class C3D15Test : public ElementTest<C3D15, 15> {};
class C3D20Test : public ElementTest<C3D20, 20> {};

TEST_F(C3D4Test, TestImplementation) {
    ASSERT_TRUE(this->element->test_implementation<C3D4>());
}

TEST_F(C3D6Test, TestImplementation) {
    ASSERT_TRUE(this->element->test_implementation<C3D6>());
}

TEST_F(C3D8Test, TestImplementation) {
    ASSERT_TRUE(this->element->test_implementation<C3D8>());
}

TEST_F(C3D10Test, TestImplementation) {
    ASSERT_TRUE(this->element->test_implementation<C3D10>());
}

TEST_F(C3D15Test, TestImplementation) {
    ASSERT_TRUE(this->element->test_implementation<C3D15>());
}

TEST_F(C3D20Test, TestImplementation) {
    ASSERT_TRUE(this->element->test_implementation<C3D20>());
}
