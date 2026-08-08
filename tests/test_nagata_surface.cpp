/**
 * @file test_nagata_surface.cpp
 * @brief Tests Nagata reconstruction, evaluation, continuity and projection.
 *
 * The suite exercises the production Nagata and BVH implementation through a
 * lightweight finite-element surface test double. It covers supported FE
 * topologies, exact planar reduction, vertex interpolation, source-coordinate
 * mapping, analytical derivatives, G1 edge continuity and closest-point
 * projection without constructing a contact constraint.
 *
 * @see fem::model::NagataSurface
 *
 * @author Finn Eggers
 * @date 08.08.2026
 */

#include "../src/model/geometry/surface/nagata_surface.h"

#include <Eigen/Geometry>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using fem::DynamicMatrix;
using fem::DynamicVector;
using fem::ID;
using fem::Index;
using fem::Mat3;
using fem::Precision;
using fem::Vec2;
using fem::Vec3;
using fem::model::Field;
using fem::model::FieldDomain;
using fem::model::NagataSurface;
using fem::model::SurfaceInterface;
using fem::model::SurfaceRegion;

constexpr Precision tight_tolerance      = Precision(2e-11);
constexpr Precision geometry_tolerance   = Precision(2e-9);
constexpr Precision derivative_tolerance = Precision(2e-5);

/**
 * @brief Configurable FE surface used only as input to the production reconstruction.
 *
 * Natural coordinates follow FEMaster S3/S4/S6/S8 ordering. Shape functions
 * implement the corresponding standard interpolation, while nodal normals can
 * be prescribed independently to produce curved Nagata patches from simple
 * positional test meshes.
 */
class TestSurface final : public SurfaceInterface {
public:
    TestSurface(std::vector<ID> node_ids,
                DynamicMatrix  natural_coordinates,
                std::vector<Vec3> nodal_normals)
        : SurfaceInterface(edge_count(static_cast<Index>(node_ids.size())),
                           static_cast<Index>(node_ids.size()),
                           node_ids.size() > 4 ? 3 : 2),
          node_ids_(std::move(node_ids)),
          natural_coordinates_(std::move(natural_coordinates)),
          nodal_normals_(std::move(nodal_normals)) {
        require(node_ids_.size() == static_cast<std::size_t>(natural_coordinates_.rows()),
            "TestSurface natural-coordinate count mismatch");
        require(node_ids_.size() == nodal_normals_.size(),
            "TestSurface normal count mismatch");
    }

    ID* nodes() override {
        return node_ids_.data();
    }

    const ID* nodes() const override {
        return node_ids_.data();
    }

    Vec3 local_to_global(const Vec2& local, const Field& node_coords) const override {
        const DynamicVector shape = shape_function(local);
        Vec3 position = Vec3::Zero();

        for (Index local_node = 0; local_node < n_nodes; ++local_node) {
            position += shape(static_cast<Eigen::Index>(local_node))
                * node_coords.row_vec3(static_cast<Index>(node_ids_[local_node]));
        }

        return position;
    }

    Vec2 global_to_local(const Vec3&, const Field&, bool) const override {
        return Vec2::Zero();
    }

    Vec3 normal(const Field&, const Vec2& local) const override {
        const DynamicVector shape = shape_function(local);
        Vec3 result = Vec3::Zero();

        for (Index local_node = 0; local_node < n_nodes; ++local_node) {
            result += shape(static_cast<Eigen::Index>(local_node))
                * nodal_normals_[static_cast<std::size_t>(local_node)];
        }

        return result.normalized();
    }

    bool in_bounds(const Vec2& local) const override {
        if (n_nodes == 3 || n_nodes == 6) {
            return local(0) >= Precision(0)
                && local(1) >= Precision(0)
                && local.sum() <= Precision(1);
        }

        return std::abs(local(0)) <= Precision(1)
            && std::abs(local(1)) <= Precision(1);
    }

    Precision area(const Field&) const override {
        return Precision(0);
    }

    Polygon local_domain_polygon() const override {
        if (n_nodes == 3 || n_nodes == 6) {
            return Polygon{
                Vec2(Precision(0), Precision(0)),
                Vec2(Precision(1), Precision(0)),
                Vec2(Precision(0), Precision(1))
            };
        }

        return Polygon{
            Vec2(Precision(-1), Precision(-1)),
            Vec2(Precision( 1), Precision(-1)),
            Vec2(Precision( 1), Precision( 1)),
            Vec2(Precision(-1), Precision( 1))
        };
    }

    DynamicVector shape_function(const Vec2& local) const override {
        const Precision r = local(0);
        const Precision s = local(1);
        DynamicVector shape(static_cast<Eigen::Index>(n_nodes));

        if (n_nodes == 3) {
            shape << Precision(1) - r - s, r, s;
        } else if (n_nodes == 4) {
            shape <<
                Precision(0.25) * (Precision(1) - r) * (Precision(1) - s),
                Precision(0.25) * (Precision(1) + r) * (Precision(1) - s),
                Precision(0.25) * (Precision(1) + r) * (Precision(1) + s),
                Precision(0.25) * (Precision(1) - r) * (Precision(1) + s);
        } else if (n_nodes == 6) {
            const Precision l0 = Precision(1) - r - s;
            shape <<
                l0 * (Precision(2)*l0 - Precision(1)),
                r  * (Precision(2)*r  - Precision(1)),
                s  * (Precision(2)*s  - Precision(1)),
                Precision(4) * l0*r,
                Precision(4) * r*s,
                Precision(4) * s*l0;
        } else {
            shape <<
                -Precision(0.25)*(Precision(1)-r)*(Precision(1)-s)*(Precision(1)+r+s),
                -Precision(0.25)*(Precision(1)+r)*(Precision(1)-s)*(Precision(1)-r+s),
                -Precision(0.25)*(Precision(1)+r)*(Precision(1)+s)*(Precision(1)-r-s),
                -Precision(0.25)*(Precision(1)-r)*(Precision(1)+s)*(Precision(1)+r-s),
                 Precision(0.5)*(Precision(1)-r*r)*(Precision(1)-s),
                 Precision(0.5)*(Precision(1)+r)*(Precision(1)-s*s),
                 Precision(0.5)*(Precision(1)-r*r)*(Precision(1)+s),
                 Precision(0.5)*(Precision(1)-r)*(Precision(1)-s*s);
        }

        return shape;
    }

    DynamicMatrix node_coords_natural() const override {
        return natural_coordinates_;
    }

    Precision integrate_scalar_field(const Field&, const ScalarField&) const override {
        return Precision(0);
    }

    Vec3 integrate_vector_field(const Field&, const VecField&) const override {
        return Vec3::Zero();
    }

    void integrate_vector_field(const Field&, Field&, const VecField&) const override {}

    Mat3 integrate_tensor_field(const Field&, const TenField&) const override {
        return Mat3::Zero();
    }

    void integrate_triangular(
        const Field&,
        const Polygon&,
        const fem::math::quadrature::Quadrature&,
        const std::function<void(const Vec2&, const Vec3&, Precision)>&) const override {}

private:
    static Index edge_count(Index node_count) {
        return node_count > 4 ? node_count / 2 : node_count;
    }

    static void require(bool condition, const char* message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    std::vector<ID>   node_ids_;
    DynamicMatrix     natural_coordinates_;
    std::vector<Vec3> nodal_normals_;
};

struct Topology {
    DynamicMatrix                    natural_coordinates;
    std::vector<std::array<Index, 3>> triangles;
};

Topology topology(Index node_count) {
    Topology result;

    if (node_count == 3) {
        result.natural_coordinates.resize(3, 2);
        result.natural_coordinates << 0, 0, 1, 0, 0, 1;
        result.triangles = {{0, 1, 2}};
    } else if (node_count == 4) {
        result.natural_coordinates.resize(4, 2);
        result.natural_coordinates << -1, -1, 1, -1, 1, 1, -1, 1;
        result.triangles = {{0, 1, 2}, {0, 2, 3}};
    } else if (node_count == 6) {
        result.natural_coordinates.resize(6, 2);
        result.natural_coordinates << 0, 0, 1, 0, 0, 1, 0.5, 0, 0.5, 0.5, 0, 0.5;
        result.triangles = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}};
    } else if (node_count == 8) {
        result.natural_coordinates.resize(8, 2);
        result.natural_coordinates <<
            -1, -1, 1, -1, 1, 1, -1, 1, 0, -1, 1, 0, 0, 1, -1, 0;
        result.triangles = {
            {0, 4, 7}, {1, 5, 4}, {2, 6, 5},
            {3, 7, 6}, {4, 5, 6}, {4, 6, 7}
        };
    } else {
        throw std::runtime_error("Unsupported test topology");
    }

    return result;
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_near(Precision actual,
                  Precision expected,
                  Precision tolerance,
                  const std::string& message) {
    require(std::abs(actual - expected) <= tolerance,
        message + ": error=" + std::to_string(std::abs(actual - expected)));
}

void require_near(const Vec2& actual,
                  const Vec2& expected,
                  Precision tolerance,
                  const std::string& message) {
    require((actual - expected).norm() <= tolerance,
        message + ": error=" + std::to_string((actual - expected).norm()));
}

void require_near(const Vec3& actual,
                  const Vec3& expected,
                  Precision tolerance,
                  const std::string& message) {
    require((actual - expected).norm() <= tolerance,
        message + ": error=" + std::to_string((actual - expected).norm()));
}

template<class Callable>
void require_throws(Callable&& callable, const std::string& message) {
    bool threw = false;

    try {
        callable();
    } catch (const std::exception&) {
        threw = true;
    }

    require(threw, message);
}

std::vector<Vec3> constant_normals(Index count, const Vec3& normal = Vec3(0, 0, 1)) {
    return std::vector<Vec3>(static_cast<std::size_t>(count), normal.normalized());
}

Field planar_coordinates(const Topology& topology_data, bool quadrilateral) {
    const Index node_count = static_cast<Index>(topology_data.natural_coordinates.rows());
    Field coordinates("POSITION", FieldDomain::NODE, node_count, 3);

    for (Index node = 0; node < node_count; ++node) {
        const Precision r = topology_data.natural_coordinates(node, 0);
        const Precision s = topology_data.natural_coordinates(node, 1);

        coordinates(node, 0) = quadrilateral ? Precision(0.5)*(r + Precision(1)) : r;
        coordinates(node, 1) = quadrilateral ? Precision(0.5)*(s + Precision(1)) : s;
        coordinates(node, 2) = Precision(0);
    }

    return coordinates;
}

NagataSurface make_single_surface(Index node_count,
                                  const Field& coordinates,
                                  const std::vector<Vec3>& normals,
                                  SurfaceInterface::Ptr* source = nullptr) {
    const Topology topology_data = topology(node_count);
    std::vector<ID> node_ids(static_cast<std::size_t>(node_count));

    for (Index node = 0; node < node_count; ++node) {
        node_ids[static_cast<std::size_t>(node)] = static_cast<ID>(node);
    }

    auto surface = std::make_shared<TestSurface>(
        std::move(node_ids), topology_data.natural_coordinates, normals);

    std::vector<SurfaceInterface::Ptr> surfaces{surface};
    SurfaceRegion surface_set("NAGATA_TEST_SURFACE");
    surface_set.add(0);

    if (source != nullptr) {
        *source = surface;
    }

    return NagataSurface(surface_set, surfaces, coordinates);
}

void test_planar_topology(Index node_count) {
    const Topology topology_data = topology(node_count);
    const bool quadrilateral = node_count == 4 || node_count == 8;
    const Field coordinates = planar_coordinates(topology_data, quadrilateral);
    SurfaceInterface::Ptr source;
    const NagataSurface surface = make_single_surface(
        node_count, coordinates, constant_normals(node_count), &source);

    const std::array<Vec2, 4> samples{
        Vec2(Precision(1.0/3.0), Precision(1.0/3.0)),
        Vec2(Precision(0.2), Precision(0.2)),
        Vec2(Precision(0.6), Precision(0.2)),
        Vec2(Precision(0.2), Precision(0.6))
    };

    for (Index patch = 0; patch < static_cast<Index>(topology_data.triangles.size()); ++patch) {
        const auto& triangle = topology_data.triangles[static_cast<std::size_t>(patch)];

        for (const Vec2& local : samples) {
            const auto state = surface.evaluate({patch, local});
            require_near(state.position(2), Precision(0), tight_tolerance,
                "planar patch position");
            require_near(state.normal, Vec3(0, 0, 1), tight_tolerance,
                "planar patch normal");
            require_near(state.d2_xixi, Vec3::Zero(), tight_tolerance,
                "planar xixi derivative");
            require_near(state.d2_xieta, Vec3::Zero(), tight_tolerance,
                "planar xieta derivative");
            require_near(state.d2_etaeta, Vec3::Zero(), tight_tolerance,
                "planar etaeta derivative");
            require(state.surface == source.get(), "source surface mapping changed");
        }

        const std::array<Vec2, 3> patch_vertices{
            Vec2(0, 0), Vec2(1, 0), Vec2(0, 1)
        };

        for (Index local_vertex = 0; local_vertex < 3; ++local_vertex) {
            const Index source_node = triangle[static_cast<std::size_t>(local_vertex)];
            const auto state = surface.evaluate({
                patch, patch_vertices[static_cast<std::size_t>(local_vertex)]});

            require_near(state.position, coordinates.row_vec3(source_node), tight_tolerance,
                "patch vertex interpolation");
            require_near(state.normal, Vec3(0, 0, 1), tight_tolerance,
                "patch vertex normal");
            require_near(
                state.element_local,
                topology_data.natural_coordinates.row(
                    static_cast<Eigen::Index>(source_node)).transpose(),
                tight_tolerance,
                "patch vertex FE mapping");
        }

        const Vec2 local(Precision(0.2), Precision(0.3));
        const Vec3 beta(Precision(0.5), Precision(0.2), Precision(0.3));
        const Vec2 expected_element_local =
            beta(0) * topology_data.natural_coordinates.row(triangle[0]).transpose()
            + beta(1) * topology_data.natural_coordinates.row(triangle[1]).transpose()
            + beta(2) * topology_data.natural_coordinates.row(triangle[2]).transpose();

        require_near(surface.evaluate({patch, local}).element_local,
            expected_element_local, tight_tolerance, "interior FE mapping");
    }

    require_throws([&] {
        (void) surface.evaluate({
            static_cast<Index>(topology_data.triangles.size()),
            Vec2(Precision(1.0/3.0), Precision(1.0/3.0))});
    }, "unexpected Nagata patch count");
}

void test_curved_derivatives() {
    const Topology topology_data = topology(3);
    Field coordinates("POSITION", FieldDomain::NODE, 3, 3);
    coordinates(0, 0) = 0.0; coordinates(0, 1) = 0.0; coordinates(0, 2) = 0.00;
    coordinates(1, 0) = 1.1; coordinates(1, 1) = 0.0; coordinates(1, 2) = 0.08;
    coordinates(2, 0) = 0.0; coordinates(2, 1) = 1.0; coordinates(2, 2) = -0.04;

    const std::vector<Vec3> normals{
        Vec3(-0.25, -0.15, 1.0).normalized(),
        Vec3( 0.20, -0.10, 1.0).normalized(),
        Vec3(-0.15,  0.25, 1.0).normalized()
    };

    const NagataSurface surface = make_single_surface(3, coordinates, normals);
    const Vec2 local(Precision(0.31), Precision(0.27));
    const Precision first_step  = Precision(2e-6);
    const Precision second_step = Precision(2e-4);

    const auto state = surface.evaluate({0, local});
    const auto position = [&](Precision xi, Precision eta) {
        return surface.evaluate({0, Vec2(xi, eta)}).position;
    };

    const Vec3 fd_xi =
        (position(local(0) + first_step, local(1))
       - position(local(0) - first_step, local(1))) / (Precision(2)*first_step);
    const Vec3 fd_eta =
        (position(local(0), local(1) + first_step)
       - position(local(0), local(1) - first_step)) / (Precision(2)*first_step);

    const Vec3 center = position(local(0), local(1));
    const Vec3 fd_xixi =
        (position(local(0) + second_step, local(1))
       - Precision(2)*center
       + position(local(0) - second_step, local(1))) / (second_step*second_step);
    const Vec3 fd_etaeta =
        (position(local(0), local(1) + second_step)
       - Precision(2)*center
       + position(local(0), local(1) - second_step)) / (second_step*second_step);
    const Vec3 fd_xieta =
        (position(local(0) + second_step, local(1) + second_step)
       - position(local(0) + second_step, local(1) - second_step)
       - position(local(0) - second_step, local(1) + second_step)
       + position(local(0) - second_step, local(1) - second_step))
        / (Precision(4)*second_step*second_step);

    require_near(state.jacobian.col(0), fd_xi, derivative_tolerance,
        "curved first xi derivative");
    require_near(state.jacobian.col(1), fd_eta, derivative_tolerance,
        "curved first eta derivative");
    require_near(state.d2_xixi, fd_xixi, derivative_tolerance,
        "curved second xixi derivative");
    require_near(state.d2_xieta, fd_xieta, derivative_tolerance,
        "curved second xieta derivative");
    require_near(state.d2_etaeta, fd_etaeta, derivative_tolerance,
        "curved second etaeta derivative");

    const Vec3 expected_normal =
        state.jacobian.col(0).cross(state.jacobian.col(1)).normalized();
    require_near(state.normal, expected_normal, tight_tolerance,
        "normal does not match Jacobian cross product");
    require_near(state.normal.norm(), Precision(1), tight_tolerance,
        "normal is not unit length");
    require_near(state.normal.dot(state.jacobian.col(0)), Precision(0), tight_tolerance,
        "normal is not orthogonal to xi tangent");
    require_near(state.normal.dot(state.jacobian.col(1)), Precision(0), tight_tolerance,
        "normal is not orthogonal to eta tangent");
}

void test_g1_shared_edge() {
    Field coordinates("POSITION", FieldDomain::NODE, 4, 3);
    coordinates(0, 0) = -1.0; coordinates(0, 1) =  0.0; coordinates(0, 2) =  0.00;
    coordinates(1, 0) =  1.0; coordinates(1, 1) =  0.0; coordinates(1, 2) =  0.05;
    coordinates(2, 0) =  0.0; coordinates(2, 1) =  1.0; coordinates(2, 2) =  0.18;
    coordinates(3, 0) =  0.0; coordinates(3, 1) = -1.0; coordinates(3, 2) = -0.12;

    const Vec3 n0 = Vec3(-0.025, -0.12, 1.0).normalized();
    const Vec3 n1 = Vec3(-0.025,  0.10, 1.0).normalized();
    const Vec3 n2 = Vec3(-0.040, -0.20, 1.0).normalized();
    const Vec3 n3 = Vec3(-0.010,  0.18, 1.0).normalized();

    const DynamicMatrix natural = topology(3).natural_coordinates;
    std::vector<SurfaceInterface::Ptr> surfaces{
        std::make_shared<TestSurface>(
            std::vector<ID>{0, 1, 2}, natural, std::vector<Vec3>{n0, n1, n2}),
        std::make_shared<TestSurface>(
            std::vector<ID>{1, 0, 3}, natural, std::vector<Vec3>{n1, n0, n3})
    };

    SurfaceRegion surface_set("G1_EDGE");
    surface_set.add(0);
    surface_set.add(1);
    const NagataSurface surface(surface_set, surfaces, coordinates);

    for (Index sample = 1; sample < 20; ++sample) {
        const Precision t = static_cast<Precision>(sample) / Precision(20);
        const auto side_a = surface.evaluate({0, Vec2(t, 0)});
        const auto side_b = surface.evaluate({1, Vec2(Precision(1) - t, 0)});

        require_near(side_a.position, side_b.position, geometry_tolerance,
            "shared-edge position discontinuity");
        require_near(Precision(1) - side_a.normal.dot(side_b.normal),
            Precision(0), geometry_tolerance, "shared-edge normal discontinuity");
    }
}

NagataSurface make_planar_square(Field& coordinates) {
    const Topology topology_data = topology(4);
    coordinates = planar_coordinates(topology_data, true);
    return make_single_surface(4, coordinates, constant_normals(4));
}

void test_planar_projection_and_boundary() {
    Field coordinates;
    const NagataSurface square = make_planar_square(coordinates);

    const Vec3 query(0.3, 0.4, 2.7);
    const auto location = square.project(query);
    const auto state    = square.evaluate(location);

    require_near(state.position, Vec3(0.3, 0.4, 0.0), geometry_tolerance,
        "planar closest-point projection");
    require_near(state.jacobian.col(0).dot(state.position - query),
        Precision(0), geometry_tolerance, "planar xi stationarity");
    require_near(state.jacobian.col(1).dot(state.position - query),
        Precision(0), geometry_tolerance, "planar eta stationarity");

    const Topology triangle_topology = topology(3);
    const Field triangle_coordinates = planar_coordinates(triangle_topology, false);
    const NagataSurface triangle = make_single_surface(
        3, triangle_coordinates, constant_normals(3));

    const auto edge = triangle.evaluate(triangle.project(Vec3(0.5, -1.0, 1.0)));
    require_near(edge.position, Vec3(0.5, 0.0, 0.0), geometry_tolerance,
        "boundary-edge projection");

    const auto corner = triangle.evaluate(triangle.project(Vec3(-1.0, -1.0, 1.0)));
    require_near(corner.position, Vec3(0.0, 0.0, 0.0), geometry_tolerance,
        "boundary-corner projection");
}

void test_projection_walks_across_patches() {
    constexpr Index element_count = 5;
    constexpr Index node_count    = 2 * (element_count + 1);

    Field coordinates("POSITION", FieldDomain::NODE, node_count, 3);

    for (Index column = 0; column <= element_count; ++column) {
        const Index bottom = 2 * column;
        const Index top    = bottom + 1;

        coordinates(bottom, 0) = static_cast<Precision>(column);
        coordinates(bottom, 1) = Precision(0);
        coordinates(bottom, 2) = Precision(0);
        coordinates(top, 0) = static_cast<Precision>(column);
        coordinates(top, 1) = Precision(1);
        coordinates(top, 2) = Precision(0);
    }

    const DynamicMatrix natural = topology(4).natural_coordinates;
    std::vector<SurfaceInterface::Ptr> surfaces;
    SurfaceRegion surface_set("PATCH_WALK_STRIP");

    for (Index element = 0; element < element_count; ++element) {
        surfaces.push_back(std::make_shared<TestSurface>(
            std::vector<ID>{
                static_cast<ID>(2*element),
                static_cast<ID>(2*(element + 1)),
                static_cast<ID>(2*(element + 1) + 1),
                static_cast<ID>(2*element + 1)
            },
            natural,
            constant_normals(4)));
        surface_set.add(static_cast<ID>(element));
    }

    const NagataSurface surface(surface_set, surfaces, coordinates);
    NagataSurface::Location location = surface.project(Vec3(0.2, 0.35, 0.4));
    Index patch_transitions = 0;
    Index surface_transitions = 0;
    const SurfaceInterface* previous_surface = surface.evaluate(location).surface;

    for (Index step = 1; step <= 48; ++step) {
        const Precision x = Precision(0.2) + Precision(0.095) * static_cast<Precision>(step);
        const Vec3 query(x, 0.35, 0.4);
        const Index previous_patch = location.patch;

        location = surface.project(query, location);
        const auto state = surface.evaluate(location);

        if (location.patch != previous_patch) {
            ++patch_transitions;
        }
        if (state.surface != previous_surface) {
            ++surface_transitions;
            previous_surface = state.surface;
        }

        require_near(state.position, Vec3(x, 0.35, 0.0), geometry_tolerance,
            "tracked planar projection");
        require_near(state.normal, Vec3(0, 0, 1), tight_tolerance,
            "tracked planar normal");
    }

    require(patch_transitions >= 5, "tracked projection did not cross enough patches");
    require(surface_transitions >= 4, "tracked projection did not cross FE surfaces");
}

void test_invalid_inputs() {
    Field coordinates("POSITION", FieldDomain::NODE, 3, 3);
    coordinates.set_zero();
    const auto source = std::make_shared<TestSurface>(
        std::vector<ID>{0, 1, 2}, topology(3).natural_coordinates, constant_normals(3));
    const std::vector<SurfaceInterface::Ptr> surfaces{source};
    const SurfaceRegion empty_set("EMPTY");

    require_throws([&] {
        (void) NagataSurface(empty_set, surfaces, coordinates);
    }, "empty surface set was accepted");

    SurfaceRegion surface_set("ZERO_EDGE");
    surface_set.add(0);

    require_throws([&] {
        (void) NagataSurface(surface_set, surfaces, coordinates);
    }, "zero-length patch edge was accepted");
}

} // namespace

TEST(NagataSurface, PlanarS3IsExact) {
    test_planar_topology(3);
}

TEST(NagataSurface, PlanarS4IsExact) {
    test_planar_topology(4);
}

TEST(NagataSurface, PlanarS6IsExact) {
    test_planar_topology(6);
}

TEST(NagataSurface, PlanarS8IsExact) {
    test_planar_topology(8);
}

TEST(NagataSurface, EvaluateDerivativesMatchFiniteDifferences) {
    test_curved_derivatives();
}

TEST(NagataSurface, AdjacentPatchesAreG1Continuous) {
    test_g1_shared_edge();
}

TEST(NagataSurface, ProjectionAndBoundaryAreCorrect) {
    test_planar_projection_and_boundary();
}

TEST(NagataSurface, ProjectionWalksAcrossMultiplePatches) {
    test_projection_walks_across_patches();
}

TEST(NagataSurface, InvalidInputsAreRejected) {
    test_invalid_inputs();
}
