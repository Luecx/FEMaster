/**
 * @file register_beam_section.cpp
 * @brief Common parsing of native and Abaqus beam section descriptions.
 *
 * Both keywords create the same generalized linear BeamSection. PROFILE
 * references an existing profile; SECTION=GENERAL reads inertias directly;
 * library SECTION types derive a Profile from their geometric dimensions.
 */

#include "register_functions.h"

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../material/isotropic_elasticity.h"
#include "../../../material/material.h"
#include "../../../model/model.h"
#include "../../../section/profile_geometry.h"
#include "../../../section/section_beam.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace fem::io::reader::commands {
namespace dsl = fem::io::dsl;

namespace {

template<std::size_t N>
dsl::Segment dimensions(std::shared_ptr<std::vector<Precision>> values) {
    auto segment = dsl::Segment::make()
        .range(dsl::LineRange{}.min(1).max(1))
        .pattern(dsl::Pattern::make().fixed<Precision, N>().name("DIMENSIONS"));
    if constexpr (N == 1) {
        segment.bind([values](Precision value) { *values = {value}; });
    } else {
        segment.bind([values](const std::array<Precision, N>& data) {
            values->assign(data.begin(), data.end());
        });
    }
    return segment;
}

dsl::Segment direction_line(std::shared_ptr<Vec3> direction, bool optional) {
    return dsl::Segment::make()
        .range(dsl::LineRange{}.min(optional ? 0 : 1).max(1))
        .pattern(dsl::Pattern::make().fixed<Precision, 3>().name("N1"))
        .bind([direction](const std::array<Precision, 3>& n1) {
            *direction = Vec3{n1[0], n1[1], n1[2]};
        });
}

dsl::Segment elastic_line(std::shared_ptr<std::array<Precision, 3>> elastic) {
    return dsl::Segment::make()
        .range(dsl::LineRange{}.min(1).max(1))
        .pattern(dsl::Pattern::make()
            .one<Precision>().name("E")
            .one<Precision>().name("G")
            .one<Precision>().name("ALPHA")
                .on_missing(Precision(0)).on_empty(Precision(0)))
        .bind([elastic](Precision E, Precision G, Precision alpha) {
            *elastic = {E, G, alpha};
        });
}

template<std::size_t N>
void add_geometry_variants(dsl::Command& command,
                           const dsl::Condition& geometry_condition,
                           std::shared_ptr<std::vector<Precision>> values,
                           std::shared_ptr<Vec3> direction,
                           std::shared_ptr<std::array<Precision, 3>> elastic) {
    const auto with_material = dsl::Condition::all_of({
        geometry_condition, dsl::Condition::key_present("MATERIAL")
    });
    const auto without_material = dsl::Condition::all_of({
        geometry_condition, dsl::Condition::negate(dsl::Condition::key_present("MATERIAL"))
    });

    command.variant(dsl::Variant::make()
        .when(with_material)
        .segment(dimensions<N>(values))
        .segment(direction_line(direction, true)));

    // n1 and direct E/G/alpha both have three components. Use separate
    // variants so an omitted n1 line is not mistaken for elastic data.
    command.variant(dsl::Variant::make()
        .when(without_material)
        .segment(dimensions<N>(values))
        .segment(direction_line(direction, false))
        .segment(elastic_line(elastic)));

    command.variant(dsl::Variant::make()
        .when(without_material)
        .segment(dimensions<N>(values))
        .segment(elastic_line(elastic)));
}

} // namespace

void register_beam_section(dsl::Registry& registry, model::Model& model) {
    // The FEMaster-specific PROFILE syntax and Abaqus SECTION syntax share
    // their grammar and final BeamSection construction.
    const auto register_section = [&](const std::string& keyword_name) {
        registry.command(keyword_name, [&, keyword_name](dsl::Command& command) {
            command.allow_if(dsl::Condition::parent_is({"ROOT", "PART"}));
            command.doc("Beam section from PROFILE, SECTION=GENERAL, or library geometry (RECT, CIRC, PIPE, BOX, I, CHANNEL, L).");

            auto elset     = std::make_shared<std::string>();
            auto material  = std::make_shared<std::string>();
            auto profile   = std::make_shared<std::string>();
            auto shape     = std::make_shared<std::string>();
            auto values    = std::make_shared<std::vector<Precision>>();
            auto direction = std::make_shared<Vec3>(Vec3{0, 0, -1});
            auto elastic   = std::make_shared<std::array<Precision, 3>>();
            auto density   = std::make_shared<Precision>(Precision(-1));

            command.keyword(dsl::KeywordSpec::make()
                .key("ELSET").required()
                .key("MATERIAL").alternative("MAT").optional()
                .key("PROFILE").optional()
                .key("SECTION").optional().allowed({"GENERAL", "RECT", "CIRC", "PIPE",
                                                    "BOX", "I", "CHANNEL", "L"})
                .key("DENSITY").optional());

            command.on_enter([=](const dsl::Keys& keys) {
                logging::error(!(keys.has("PROFILE") && keys.has("SECTION")),
                    keyword_name, ": PROFILE and SECTION are mutually exclusive");

                *elset = keys.raw("ELSET");
                *material = keys.raw("MATERIAL");
                *profile = keys.raw("PROFILE");
                *shape = keys.raw("SECTION").empty() ? "GENERAL" : keys.raw("SECTION");
                values->clear();
                elastic->fill(Precision(0));
                *direction = Vec3{0, 0, -1};
                *density = keys.has("DENSITY") ? keys.get<Precision>("DENSITY") : Precision(-1);
                logging::error(!keys.has("DENSITY") || *density >= Precision(0),
                    keyword_name, ": DENSITY must be non-negative");
                logging::error(profile->empty() || !material->empty(),
                    keyword_name, ": PROFILE requires MATERIAL");
            });

            command.on_exit([=, &model](const dsl::Keys&) {
                const auto part = model._data->parts.get();
                logging::error(part != nullptr, keyword_name, ": no active part");
                logging::error(part->elem_sets.has(*elset),
                    keyword_name, ": undefined ELSET ", *elset);
                logging::error(direction->allFinite() && direction->norm() > Precision(0),
                    keyword_name, ": n1 must be finite and non-zero");

                material::MaterialPtr assigned;
                if (!material->empty()) {
                    logging::error(model._data->materials.has(*material),
                        keyword_name, ": unknown MATERIAL ", *material);
                    assigned = model._data->materials.get(*material);
                } else {
                    const Precision E = (*elastic)[0], G = (*elastic)[1];
                    logging::error(std::isfinite(E) && std::isfinite(G) && E > 0 && G > 0,
                        keyword_name, ": direct properties require positive E and G");
                    const Precision nu = E / (Precision(2) * G) - Precision(1);
                    logging::error(nu > Precision(-1) && nu < Precision(0.5),
                        keyword_name, ": E and G are incompatible with isotropic elasticity");
                    assigned = std::make_shared<material::Material>("__BEAM_" + *elset);
                    assigned->set_elasticity<material::IsotropicElasticity>(E, nu);
                    assigned->set_thermal_expansion((*elastic)[2]);
                }

                // A section-local density override cannot mutate a named material.
                if (*density >= Precision(0)) {
                    if (!material->empty()) {
                        logging::error(assigned->has_elasticity(),
                            keyword_name, ": DENSITY override requires material elasticity");
                        const auto* iso = assigned->elasticity()->as<material::IsotropicElasticity>();
                        logging::error(iso != nullptr,
                            keyword_name, ": DENSITY override requires isotropic material");
                        auto copy = std::make_shared<material::Material>("__BEAM_" + *elset);
                        copy->set_elasticity<material::IsotropicElasticity>(iso->youngs, iso->poisson);
                        if (assigned->has_thermal_expansion()) {
                            copy->set_thermal_expansion(assigned->get_thermal_expansion());
                            copy->set_thermal_zero_temperature(assigned->get_thermal_zero_temperature());
                        }
                        assigned = std::move(copy);
                    }
                    assigned->set_density(*density);
                }

                Profile::Ptr assigned_profile;
                if (!profile->empty()) {
                    logging::error(model._data->profiles.has(*profile),
                        keyword_name, ": undefined PROFILE ", *profile);
                    assigned_profile = model._data->profiles.get(*profile);
                } else if (*shape == "GENERAL") {
                    logging::error(values->size() == 5,
                        keyword_name, ": GENERAL requires A, I11, I12, I22, J");
                    const auto& v = *values;
                    logging::error(v[0] > 0 && v[1] > 0 && v[3] > 0 && v[4] > 0 &&
                                   v[1] * v[3] > v[2] * v[2],
                        keyword_name, ": invalid general section inertia values");
                    assigned_profile = std::make_shared<Profile>(
                        "__BEAM_" + *elset, v[0], v[1], v[3], v[4], v[2]);
                } else {
                    assigned_profile = make_beam_geometry_profile(
                        "__BEAM_" + *elset, *shape, *values);
                }

                auto section = std::make_shared<BeamSection>();
                section->region_    = part->elem_sets.get(*elset);
                section->material_  = std::move(assigned);
                section->profile_   = std::move(assigned_profile);
                section->direction_ = *direction;
                model.add_section(std::move(section));
            });

            // Native PROFILE: an existing profile with optional n1 vector.
            command.variant(dsl::Variant::make()
                .when(dsl::Condition::key_present("PROFILE"))
                .segment(direction_line(direction, true)));

            // Omitted SECTION is GENERAL unless PROFILE was supplied.
            const auto general_condition = dsl::Condition::all_of({
                dsl::Condition::negate(dsl::Condition::key_present("PROFILE")),
                dsl::Condition::any_of({
                    dsl::Condition::negate(dsl::Condition::key_present("SECTION")),
                    dsl::Condition::key_equals("SECTION", "GENERAL")
                })
            });
            add_geometry_variants<5>(command, general_condition, values, direction, elastic);

            const auto shape_condition = [](const std::string& kind) {
                return dsl::Condition::all_of({
                    dsl::Condition::negate(dsl::Condition::key_present("PROFILE")),
                    dsl::Condition::key_equals("SECTION", kind)
                });
            };

            add_geometry_variants<2>(command, shape_condition("RECT"), values, direction, elastic);
            add_geometry_variants<1>(command, shape_condition("CIRC"), values, direction, elastic);
            add_geometry_variants<2>(command, shape_condition("PIPE"), values, direction, elastic);
            add_geometry_variants<3>(command, shape_condition("BOX"), values, direction, elastic);
            add_geometry_variants<6>(command, shape_condition("BOX"), values, direction, elastic);
            add_geometry_variants<7>(command, shape_condition("I"), values, direction, elastic);
            add_geometry_variants<8>(command, shape_condition("CHANNEL"), values, direction, elastic);
            add_geometry_variants<4>(command, shape_condition("L"), values, direction, elastic);
        });
    };

    register_section("BEAMSECTION");
    register_section("BEAMGENERALSECTION");
}

} // namespace fem::io::reader::commands
