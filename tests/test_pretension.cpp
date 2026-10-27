#include "../src/io/reader/parser.h"
#include "../src/model/model.h"
#include "../src/constraints/transformer/constraint_transformer.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace fem;

namespace {
std::string read_file(const std::filesystem::path& path) {
    std::ifstream input(path);
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

void replace(std::string& text, const std::string& from, const std::string& to) {
    const auto position = text.find(from);
    ASSERT_NE(position, std::string::npos) << from;
    text.replace(position, from.size(), to);
}

using Rows = std::vector<std::vector<double>>;
using Frame = std::map<std::string, Rows>;
std::vector<Frame> read_frames(const std::string& text) {
    std::vector<Frame> frames;
    std::istringstream input(text);
    std::string line, field;
    while (std::getline(input, line)) {
        if (line.rfind("LC ", 0) == 0) frames.emplace_back();
        else if (line.rfind("FIELD, NAME=", 0) == 0) {
            const auto end = line.find(',', 12);
            field = line.substr(12, end - 12);
        } else if (line == "END FIELD") field.clear();
        else if (!field.empty() && !frames.empty()) {
            std::istringstream row(line);
            std::string value;
            row >> value;
            std::vector<double> values;
            while (row >> value) values.push_back(std::strtod(value.c_str(), nullptr));
            frames.back()[field].push_back(std::move(values));
        }
    }
    return frames;
}

class Pretension : public testing::Test {
protected:
    std::filesystem::path directory;
    std::string deck;
    io::reader::Parser parser;

    void SetUp() override {
        directory = std::filesystem::temp_directory_path() /
            (std::string("femaster_pretension_") + testing::UnitTest::GetInstance()->current_test_info()->name());
        std::filesystem::create_directories(directory);
        deck = read_file(std::filesystem::path(__FILE__).parent_path() / "data/pretension_face_pair_loadsteps.inp");
        replace(deck, "210000.0, 0.3", "210000.0, 0.0");
    }
    void TearDown() override {
        std::filesystem::remove(directory / "model.inp");
        std::filesystem::remove(directory / "model.res");
        std::filesystem::remove(directory / "model.frd");
        std::filesystem::remove(directory / "model.femr");
        std::filesystem::remove(directory);
    }
    void run() {
        std::ofstream(directory / "model.inp") << deck;
        parser.run((directory / "model.inp").string(), (directory / "model").string(), {true, true, false});
    }
};

TEST_F(Pretension, ForceHoldExternalLoad) {
    run();
    const auto& data = *parser.model()._data;
    ASSERT_EQ(data.positions->rows, 20);
    ASSERT_EQ(data.pretension_sections.size(), 1u);
    const auto& section = *data.pretension_sections.front();
    ASSERT_EQ(section.interface_pairs.size(), 4u);
    const auto frames = read_frames(read_file(directory / "model.res"));
    ASSERT_EQ(frames.size(), 3u);
    const double expected_gap = 1000.0 * 3.0 / (210000.0 * 4.0);
    for (const auto& frame : frames) {
        const auto& gap = frame.at("PTGAP");
        const auto& displacement = frame.at("DISPLACEMENT");
        for (const auto& pair : section.interface_pairs) {
            EXPECT_NEAR(gap.at(pair.side_a).at(0), expected_gap, 1e-8);
            EXPECT_NEAR(gap.at(pair.side_b).at(0), expected_gap, 1e-8);
            Vec3 relative;
            for (Dim c = 0; c < 3; ++c)
                relative(c) = displacement.at(pair.side_b).at(c) - displacement.at(pair.side_a).at(c);
            EXPECT_NEAR((relative - expected_gap * section.axis_direction).norm(), 0, 1e-8);
        }
    }
    const auto cut_force = [&](const Frame& frame) {
        double total = 0;
        for (const auto& pair : section.interface_pairs)
            for (Dim c = 0; c < 3; ++c)
                total += frame.at("PTFORC").at(pair.side_b).at(c) * section.axis_direction(c);
        return total;
    };
    EXPECT_NEAR(cut_force(frames[0]), 1000.0, 1e-3);
    EXPECT_NEAR(cut_force(frames[1]), 1000.0, 1e-3);
    EXPECT_GT(std::abs(cut_force(frames[2]) - cut_force(frames[1])), 0.1);
    EXPECT_NEAR(section.locked_gap, expected_gap, 1e-8);
}

TEST_F(Pretension, RejectsSharedTopology) {
    replace(deck, "2, 17, 18, 19, 20, 9, 10, 11, 12",
        "2, 5, 6, 7, 8, 9, 10, 11, 12");
    EXPECT_THROW(run(), std::exception);
}
TEST_F(Pretension, RejectsMissingMatch) {
    replace(deck, "17, -1, -1, 1", "17, -1, -1, 1.1");
    EXPECT_THROW(run(), std::exception);
}
TEST_F(Pretension, RejectsAmbiguousMatch) {
    replace(deck, "18,  1, -1, 1", "18, -1, -1, 1");
    EXPECT_THROW(run(), std::exception);
}
TEST_F(Pretension, RejectsHoldWithoutCompletedStep) {
    replace(deck, "*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=1000.0",
        "*PRETENSION, SECTION=PT1, MODE=HOLD");
    EXPECT_THROW(run(), std::exception);
}
TEST(PretensionPostCheck, RejectsUnbalancedAndNonfiniteSolutions) {
    SystemDofIds dofs(1, 6);
    dofs.setConstant(-1);
    dofs(0, 0) = 0;
    constraint::ConstraintTransformer::Options options;
    options.method = constraint::ConstraintTransformer::Method::NullSpace;
    constraint::ConstraintTransformer transformer({}, dofs, 1, options);
    SparseMatrix stiffness(1, 1);
    stiffness.insert(0, 0) = 1;
    const DynamicVector load = DynamicVector::Ones(1);
    EXPECT_TRUE(transformer.post_check_static(stiffness, load, load));
    EXPECT_FALSE(transformer.post_check_static(stiffness, load, DynamicVector::Zero(1)));
    EXPECT_FALSE(transformer.post_check_static(stiffness, load,
        DynamicVector::Constant(1, std::numeric_limits<Precision>::quiet_NaN())));
}

TEST_F(Pretension, RejectsFailedPostCheckBeforeCapturingGap) {
    replace(deck, "VALUE=1000.0", "VALUE=1e200");
    try {
        run();
        FAIL() << "Expected the numerical post-check to reject the step";
    } catch (const std::exception& error) {
        EXPECT_NE(std::string(error.what()).find("numerical post-check failed"), std::string::npos);
    }
    const auto& section = *parser.model()._data->pretension_sections.front();
    EXPECT_FALSE(section.has_solved_gap);
    EXPECT_EQ(section.last_solved_gap, 0);
    EXPECT_THROW(parser.model().lock_pretension_section(section.name), std::exception);
    EXPECT_TRUE(read_frames(read_file(directory / "model.res")).empty());
}
TEST_F(Pretension, FailedStepPreservesPreviouslyAcceptedGap) {
    replace(deck, "*PRETENSION, SECTION=PT1, MODE=HOLD",
        "*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=1e200");
    try {
        run();
        FAIL() << "Expected the second step to fail its numerical post-check";
    } catch (const std::exception& error) {
        EXPECT_NE(std::string(error.what()).find("numerical post-check failed"), std::string::npos);
    }
    const auto& section = *parser.model()._data->pretension_sections.front();
    const double expected_gap = 1000.0 * 3.0 / (210000.0 * 4.0);
    EXPECT_TRUE(section.has_solved_gap);
    EXPECT_NEAR(section.last_solved_gap, expected_gap, 1e-8);
    ASSERT_NO_THROW(parser.model().lock_pretension_section(section.name));
    EXPECT_NEAR(section.locked_gap, expected_gap, 1e-8);
    EXPECT_EQ(read_frames(read_file(directory / "model.res")).size(), 1u);
}
} // namespace