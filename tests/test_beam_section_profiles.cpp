#include "../src/io/dsl/deck_parser.h"
#include "../src/io/dsl/file.h"
#include "../src/io/reader/parser.h"
#include "../src/section/profile_geometry.h"
#include "../src/section/section_beam.h"
#include "../src/model/model.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

using namespace fem;

TEST(BeamSectionProfiles, RectAndCircularProperties) {
    const auto rect = make_beam_geometry_profile("RECT_TEST", "RECT", {20, 40});
    EXPECT_DOUBLE_EQ(rect->area_, 800);
    EXPECT_NEAR(rect->inertia_y_, 20.0 * std::pow(40.0, 3) / 12, 0.1);
    EXPECT_NEAR(rect->inertia_z_, 40.0 * std::pow(20.0, 3) / 12, 0.1);

    const auto circle = make_beam_geometry_profile("CIRC_TEST", "CIRC", {5});
    EXPECT_NEAR(circle->area_, std::acos(-1.0) * 25, 1e-4);
    EXPECT_NEAR(circle->torsion_inertia_, 2 * circle->inertia_y_, 1e-4);

    const auto pipe = make_beam_geometry_profile("PIPE_TEST", "PIPE", {5, 1});
    EXPECT_LT(pipe->area_, circle->area_);
    EXPECT_NEAR(pipe->torsion_inertia_, 2 * pipe->inertia_y_, 1e-4);
}

TEST(BeamSectionProfiles, BoxAndOpenSections) {
    const auto box = make_beam_geometry_profile("BOX_TEST", "BOX", {100, 50, 5});
    EXPECT_NEAR(box->area_, 100.0 * 50 - 90.0 * 40, 1e-4);
    EXPECT_GT(box->torsion_inertia_, 0);

    const auto unequal = make_beam_geometry_profile("BOX_UNEQUAL", "BOX", {100, 50, 4, 5, 6, 7});
    EXPECT_GT(unequal->area_, 0);
    EXPECT_GT(unequal->inertia_y_, 0);
    EXPECT_GT(unequal->inertia_z_, 0);

    const auto ibeam = make_beam_geometry_profile("I_TEST", "I", {50, 100, 60, 60, 10, 10, 8});
    EXPECT_NEAR(ibeam->area_, 2 * 60.0 * 10 + 80.0 * 8, 1e-4);
    EXPECT_NEAR(ibeam->product_inertia_yz_, 0, 1e-4);

    const auto channel = make_beam_geometry_profile(
        "CHANNEL_TEST", "CHANNEL", {50, 100, 60, 60, 10, 10, 8, 4});
    EXPECT_NEAR(channel->area_, ibeam->area_, 1e-4);
    EXPECT_GT(channel->torsion_inertia_, 0);

    const auto angle = make_beam_geometry_profile("L_TEST", "L", {60, 80, 8, 6});
    EXPECT_NEAR(angle->area_, 60.0 * 8 + 72.0 * 6, 1e-4);
    EXPECT_GT(angle->inertia_y_ * angle->inertia_z_ -
              angle->product_inertia_yz_ * angle->product_inertia_yz_, 0);
}

TEST(BeamSectionReader, BothKeywordsAcceptProfilesAndGeometry) {
    const std::string path = "tests/TMP_BEAM_SECTION_VARIANTS.INP";
    {
        std::ofstream out(path);
        ASSERT_TRUE(out.is_open());
        out << "*BEAM SECTION, ELSET=BEAMS, MATERIAL=STEEL, PROFILE=P1\n";
        out << "0, 0, -1\n";
        out << "*BEAM GENERAL SECTION, ELSET=BEAMS, MATERIAL=STEEL, PROFILE=P1\n";
        out << "0, 0, -1\n";
        out << "*BEAM SECTION, ELSET=BEAMS, MATERIAL=STEEL, SECTION=RECT\n";
        out << "20, 40\n";
        out << "0, 0, -1\n";
        out << "*BEAM GENERAL SECTION, ELSET=BEAMS, MATERIAL=STEEL, SECTION=PIPE\n";
        out << "10, 2\n";
        out << "0, 0, -1\n";
        out << "*BEAM GENERAL SECTION, ELSET=BEAMS, SECTION=GENERAL\n";
        out << "100, 50, 0, 60, 25\n";
        out << "210000, 80000, 0\n";
        out << "*BEAM SECTION, ELSET=BEAMS, SECTION=I\n";
        out << "50, 100, 60, 60, 10, 10, 8\n";
        out << "0, 0, -1\n";
        out << "210000, 80000, 0\n";
    }

    io::reader::Parser parser;
    io::dsl::File file(path);
    io::dsl::DeckParser deck_parser(parser.registry());
    const auto deck = deck_parser.parse(file);

    EXPECT_EQ(deck.root().children("BEAMSECTION").size(), 3u);
    EXPECT_EQ(deck.root().children("BEAMGENERALSECTION").size(), 3u);

    std::filesystem::remove(path);
}

TEST(BeamSectionReader, RejectsUnknownGeometry) {
    const std::string path = "tests/TMP_BEAM_SECTION_INVALID.INP";
    {
        std::ofstream out(path);
        ASSERT_TRUE(out.is_open());
        out << "*BEAM SECTION, ELSET=BEAMS, MATERIAL=STEEL, SECTION=UNSUPPORTED\n";
        out << "1, 1\n";
    }
    io::reader::Parser parser;
    io::dsl::File file(path);
    io::dsl::DeckParser deck_parser(parser.registry());
    EXPECT_THROW(deck_parser.parse(file), std::runtime_error);
    std::filesystem::remove(path);
}

TEST(BeamSectionReader, BothKeywordsConstructBeamSections) {
    const std::string input = "tests/TMP_BEAM_SECTIONS_COMPLETE.INP";
    const std::string output = "tests/TMP_BEAM_SECTIONS_COMPLETE.RES";
    {
        std::ofstream out(input);
        ASSERT_TRUE(out.is_open());
        out << "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.3\n";
        out << "*PROFILE, NAME=P1\n100, 800, 600, 400\n";
        out << "*NODE\n1, 0, 0, 0\n2, 10, 0, 0\n3, 20, 0, 0\n4, 30, 0, 0\n5, 40, 0, 0\n";
        out << "*ELEMENT, TYPE=B31\n1, 1, 2\n2, 2, 3\n3, 3, 4\n4, 4, 5\n";
        for (int id = 1; id <= 4; ++id) {
            out << "*ELSET, ELSET=SET" << id << "\n" << id << "\n";
        }
        out << "*BEAM SECTION, ELSET=SET1, MATERIAL=STEEL, PROFILE=P1\n0, 0, -1\n";
        out << "*BEAM SECTION, ELSET=SET2, MATERIAL=STEEL, SECTION=RECT\n10, 20\n0, 0, -1\n";
        out << "*BEAM GENERAL SECTION, ELSET=SET3, SECTION=GENERAL\n80, 700, 0, 900, 300\n210000, 80000, 0\n";
        out << "*BEAM GENERAL SECTION, ELSET=SET4, MATERIAL=STEEL, SECTION=PIPE\n10, 2\n0, 0, -1\n";
    }

    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input, output));
    const auto part = parser.model()._data->parts.get();
    ASSERT_NE(part, nullptr);
    ASSERT_EQ(part->sections.size(), 4u);
    for (const auto& section : part->sections) {
        auto* beam = section->as<BeamSection>();
        ASSERT_NE(beam, nullptr);
        ASSERT_NE(beam->profile_, nullptr);
        EXPECT_GT(beam->profile_->area_, 0);
        ASSERT_NE(beam->material_, nullptr);
    }

    std::filesystem::remove(input);
    std::filesystem::remove(output);
}
