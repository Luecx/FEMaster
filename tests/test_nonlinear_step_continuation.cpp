/**
 * @file test_nonlinear_step_continuation.cpp
 * @brief Regressions for Newton residual scaling after an inherited nonlinear load.
 */

#include "../src/io/reader/parser.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

using namespace fem;

namespace {

// Exercise the actual STEP/parser history and nonlinear solver, not only the
// individual load factor functions. Both regressions used to fail at the first
// Newton evaluation of the continuation step with a zero load derivative.
void run_continuation(const std::string& label,
                      const std::string& steps,
                      Precision          expected_displacement) {
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
