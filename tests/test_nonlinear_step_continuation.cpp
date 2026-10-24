/**
 * @file test_nonlinear_step_continuation.cpp
 * @brief Regression tests for nonlinear continuation and static load amplitudes.
 */

#include "../src/io/reader/parser.h"
#include "../src/loadcase/tools/load_control.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>

#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <stdexcept>
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


TEST(NonlinearContinuation, NewlyAmplitudedLoadRemainsConstantInNextStep) {
    // A nonlinear load-control step follows the intermediate amplitude values.
    // Its completed load is then inherited without evaluating the amplitude
    // again during the next step.
    run_continuation("AMPLITUDE_CONTINUED", R"(*AMPLITUDE, NAME=CURVE
0., 0., 0.5, 0.25, 1., 1.
*STEP, NAME=PRELOAD, NLGEOM=YES, INC=40
*STATIC
0.25, 1., 1.e-5, 0.25
*CLOAD, AMPLITUDE=CURVE
2, 1, 1000.
*END STEP
*STEP, NAME=UNCHANGED, NLGEOM=YES, INC=40
*STATIC
0.25, 1., 1.e-5, 0.25
*END STEP
)", Precision(0.047616));
}

TEST(NonlinearContinuation, AmplitudeEndpointIsSolvedWithinLastIncrement) {
    // The old tolerance-based termination accepted lambda=0.9995 with
    // zero load and skipped the full force at lambda=1. With INC=1 the
    // pending step must be snapped to the endpoint before Newton starts.
    run_continuation("AMPLITUDE_EXACT_END", R"(*AMPLITUDE, NAME=CURVE
0., 0., 0.9995, 0., 1., 1.
*STEP, NAME=ENDPOINT, NLGEOM=YES, INC=1
*STATIC, DIRECT
0.9995, 1., 1.e-5, 1.
*CLOAD, AMPLITUDE=CURVE
2, 1, 1000.
*END STEP
)", Precision(0.047616));
}

TEST(NonlinearContinuation, AllowsFinalRemainderBelowMinimumIncrement) {
    // A final remainder smaller than MINIMUM_INCREMENT is valid when it
    // solves the exact endpoint, rather than an intermediate load factor.
    run_continuation("SMALL_FINAL_INCREMENT", R"(*STEP, NAME=ENDPOINT, NLGEOM=YES, INC=2
*STATIC, DIRECT
0.6, 1., 0.41, 0.6
*CLOAD
2, 1, 1000.
*END STEP
)", Precision(0.047616));
}

TEST(NonlinearContinuation, EndpointSnapRespectsMaximumIncrement) {
    // Reaching lambda=1 would exceed the configured maximum step, so a
    // single permitted increment must fail rather than silently overshoot.
    run_continuation("LIMITED_ENDPOINT", R"(*STEP, NAME=ENDPOINT, NLGEOM=YES, INC=1
*STATIC, DIRECT
0.9995, 1., 1.e-5, 0.9995
*CLOAD
2, 1, 1000.
*END STEP
)", Precision(0), true);
}

TEST(NonlinearContinuation, FixedTenthsReachExactEndpoint) {
    // Ten fixed increments of 0.1 can accumulate lambda=0.8999999999999999
    // after the ninth increment. The last target must still snap to exactly
    // one despite the rounded remainder marginally exceeding MAX_INCREMENT.
    run_continuation("FIXED_TENTHS", R"(*STEP, NAME=ENDPOINT, NLGEOM=YES, INC=10
*STATIC, DIRECT
0.1, 1., 0.1, 0.1
*CLOAD
2, 1, 1000.
*END STEP
)", Precision(0.047616));
}

TEST(LoadControl, EightyFixedIncrementsReachExactEndpoint) {
    loadcase::tools::LoadControl control;
    control.maximum_increments = 80;
    control.initial_increment  = Precision(0.0125);
    control.minimum_increment  = Precision(0.0125);
    control.maximum_increment  = Precision(0.0125);
    control.adaptive           = false;

    DynamicVector q = DynamicVector::Zero(1);
    Precision     lambda = Precision(0);
    Index         accepted_increments = 0;

    // Ordinary floating-point summation of eighty increments of 0.0125
    // can leave lambda slightly below one when the increment limit is hit.
    // Compensated accumulation must preserve all eighty increments and
    // converge to the exact endpoint rather than reporting MAXIMUM_INCREMENTS.
    const bool converged = control.solve(
        q, lambda,
        [](const DynamicVector& state, Precision target, DynamicVector& residual, SparseMatrix& tangent) {
            residual = DynamicVector::Constant(1, target - state(0));
            tangent.resize(1, 1);
            tangent.setIdentity();
        },
        [](const SparseMatrix&, const DynamicVector& residual) -> DynamicVector {
            return residual;
        },
        [](const DynamicVector& residual, Precision) -> Precision {
            return residual.norm();
        },
        [](const DynamicVector&, const DynamicVector&) -> Precision {
            return Precision(0);
        },
        {},
        [&](Index, const DynamicVector&, Precision) {
            ++accepted_increments;
        }
    );

    EXPECT_TRUE(converged);
    EXPECT_EQ(accepted_increments, Index(80));
    EXPECT_EQ(lambda, Precision(1));
    EXPECT_NEAR(q(0), Precision(1), Precision(1e-5));
}

TEST(LoadControl, CutbackUsesUnroundedProposedIncrement) {
    loadcase::tools::LoadControl control;
    control.maximum_increments       = 12;
    control.initial_increment        = Precision(0.2);
    control.minimum_increment        = Precision(0.1);
    control.maximum_increment        = Precision(0.2);
    control.growth_factor            = Precision(1);
    control.cutback_factor           = Precision(0.5);
    control.severe_cutback_factor    = Precision(0.5);

    DynamicVector q = DynamicVector::Zero(1);
    Precision     lambda = Precision(0);
    bool          rejected_once = false;

    // Force one retry near lambda=0.8 after three accepted increments.
    // For an unsnapped attempt, the cutback must use the proposed 0.2
    // rather than 0.8 - 0.6000000000000001, which would round 0.1
    // slightly below MINIMUM_INCREMENT and incorrectly abort.
    const bool converged = control.solve(
        q, lambda,
        [&](const DynamicVector& state, Precision target, DynamicVector& residual, SparseMatrix& tangent) {
            if (!rejected_once && std::abs(target - Precision(0.8)) < Precision(1e-5)) {
                rejected_once = true;
                throw std::runtime_error("Deliberate cutback regression");
            }

            residual = DynamicVector::Constant(1, target - state(0));
            tangent.resize(1, 1);
            tangent.setIdentity();
        },
        [](const SparseMatrix&, const DynamicVector& residual) -> DynamicVector {
            return residual;
        },
        [](const DynamicVector& residual, Precision) -> Precision {
            return residual.norm();
        },
        [](const DynamicVector&, const DynamicVector& correction) -> Precision {
            return correction.norm();
        }
    );

    EXPECT_TRUE(rejected_once);
    EXPECT_TRUE(converged);
    EXPECT_EQ(lambda, Precision(1));
    EXPECT_NEAR(q(0), Precision(1), Precision(1e-5));
}

TEST(NonlinearContinuation, RejectsIncompleteStaticAmplitudeRamp) {
    // Without a full end multiplier, nominal pointer continuation would
    // silently change the physical force between consecutive steps.
    run_continuation("AMPLITUDE_ENDPOINT", R"(*AMPLITUDE, NAME=CURVE
0., 0., 1., 0.5
*STEP, NAME=PRELOAD, NLGEOM=YES, INC=40
*STATIC
0.25, 1., 1.e-5, 0.25
*CLOAD, AMPLITUDE=CURVE
2, 1, 1000.
*END STEP
)", Precision(0), true);
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
