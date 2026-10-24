/**
 * @file test_nonlinear_step_continuation.cpp
 * @brief Regression tests for nonlinear continuation and static load amplitudes.
 */

#include "../src/io/reader/parser.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>

#include <exception>
#include <filesystem>
#include <fstream>
#include <string>

using namespace fem;

namespace {

// Exercise the actual STEP/parser history and nonlinear solver, not only the
// individual load factor functions. The initial regressions failed at the first
// Newton evaluation of the continuation step with a zero load derivative.
void run_continuation(const std::string& label,
                      const std::string& steps,
                      Precision          expected_displacement,
                      bool               expect_failure = false) {
    const std::string base = "tests/TMP_NONLINEAR_CONTINUATION_" + label;
    const std::string input_path = base + ".inp";

    // Always remove generated files, including if a solver assertion throws.
    struct Cleanup {
        std::string base;
        ~Cleanup() {
            std::filesystem::remove(base + ".inp");
            std::filesystem::remove(base + ".res");
            std::filesystem::remove(base + ".frd");
            std::filesystem::remove(base + ".femr");
        }
    } cleanup{base};

    {
        std::ofstream input(input_path);
        ASSERT_TRUE(input.is_open());
        input << R"(*NODE
1, 0., 0., 0.
2, 1000., 0., 0.
*ELEMENT, TYPE=T3, ELSET=BAR
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*TRUSS SECTION, ELSET=BAR, MATERIAL=STEEL
100.
*BOUNDARY
1, 1, 3, 0.
2, 2, 3, 0.
)" << steps;
    }

    io::reader::Parser parser;

    // Invalid load histories must fail before any nonlinear solve. The same
    // helper exercises both successful continuation and rejected step input.
    if (expect_failure) {
        EXPECT_THROW(parser.run(input_path, base + ".res"), std::exception);
        return;
    }

    ASSERT_NO_THROW(parser.run(input_path, base + ".res"));

    const auto& data = *parser.model()._data;
    ASSERT_NE(data.positions, nullptr);
    ASSERT_NE(data.positions_reference, nullptr);

    const ID tip = parser.model().compiled_node_id(2);
    const Precision displacement =
        (*data.positions)(tip, 0) - (*data.positions_reference)(tip, 0);
    EXPECT_NEAR(displacement, expected_displacement, Precision(1e-4));
}

} // namespace

TEST(NonlinearContinuation, UnchangedLoadAfterPerturbation) {
    run_continuation("PERTURBATION", R"(*STEP, NAME=PRELOAD, NLGEOM=YES, INC=40
*STATIC
0.25, 1., 1.e-5, 0.25
*CLOAD
2, 1, 1000.
*OUTPUT, FIELD
*NODE OUTPUT
U, CF
*END STEP
*STEP, NAME=PERTURBATION, PERTURBATION
*STATIC
*CLOAD, OP=MOD
2, 1, 1500.
*OUTPUT, FIELD
*NODE OUTPUT
U, CF
*END STEP
*STEP, NAME=UNCHANGED, NLGEOM=YES, INC=40
*STATIC
0.25, 1., 1.e-5, 0.25
*OUTPUT, FIELD
*NODE OUTPUT
U, CF
*END STEP
)", Precision(0.047616));
}

TEST(NonlinearContinuation, UnchangedLoadAfterRiks) {
    run_continuation("RIKS", R"(*STEP, NAME=RIKS, NLGEOM=YES, INC=80
*STATIC, RIKS
0.1, 1., 1.e-5, 0.2
*CLOAD
2, 1, 1000.
*OUTPUT, FIELD
*NODE OUTPUT
U, CF
*END STEP
*STEP, NAME=UNCHANGED, NLGEOM=YES, INC=40
*STATIC
0.2, 1., 1.e-5, 0.2
*OUTPUT, FIELD
*NODE OUTPUT
U, CF
*END STEP
*STEP, NAME=MODIFIED, NLGEOM=YES, INC=40
*STATIC
0.2, 1., 1.e-5, 0.2
*CLOAD, OP=MOD
2, 1, 1500.
*OUTPUT, FIELD
*NODE OUTPUT
U, CF
*END STEP
)", Precision(0.071421));
}

TEST(NonlinearContinuation, RiksRejectsNewStaticAmplitude) {
    run_continuation("RIKS_AMPLITUDE", R"(*AMPLITUDE, NAME=PULSE
0., 0., 1., 1.
*STEP, NAME=PRELOAD, NLGEOM=YES, INC=40
*STATIC
0.25, 1., 1.e-5, 0.25
*CLOAD, AMPLITUDE=PULSE
2, 1, 1000.
*END STEP
*STEP, NAME=RIKS, NLGEOM=YES, INC=80
*STATIC, RIKS
0.1, 1., 1.e-5, 0.2
*CLOAD, OP=MOD, AMPLITUDE=PULSE
2, 1, 1500.
*END STEP
)", Precision(0), true);
}

TEST(NonlinearContinuation, ReactionOnlyStartLoadsDoNotRelaxNewtonTolerance) {
    // The 1e12 force on the fully fixed first node is carried into step two.
    // It is a reaction load, not a driving force, and must not scale down the
    // reduced Newton residual for the 200000-N load on the free second node.
    // For the GL truss: F = A*E*(lambda^3-lambda)/2, L0 = 1000.
    run_continuation("FIXED_REACTION", R"(*STEP, NAME=PRELOAD, NLGEOM=YES, INC=40
*STATIC
0.25, 1., 1.e-5, 0.25
*CLOAD
1, 1, 1000000000000.
2, 1, 1000.
*END STEP
*STEP, NAME=STRONG_LOAD_CHANGE, NLGEOM=YES, INC=40
*STATIC
0.25, 1., 1.e-5, 0.25
*CLOAD, OP=MOD
2, 1, 200000.
*END STEP
)", Precision(9.39110609885));
}
