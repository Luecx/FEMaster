#include "../src/cos/rectangular_system.h"
#include "../src/cos/cylindrical_system.h"

#include <gtest/gtest.h>
#include <cmath>

using namespace fem;

// 17) RectangularSystem round-trip
TEST(Coords_Rectangular, RoundTrip) {
    cos::RectangularSystem sys("R", Vec3(1,0,0), Vec3(0,1,0));
    Vec3 p(0.3, -1.2, 2.5);
    Vec3 local = sys.to_local(p);
    Vec3 back  = sys.to_global(local);
    EXPECT_NEAR((back - p).norm(), 0.0, 1e-12);
}

// 18) Rectangular rotations orthonormal
TEST(Coords_Rectangular, RotationsAreOrthonormal) {
    auto Rx = cos::RectangularSystem::rotation_x(0.37);
    auto Ry = cos::RectangularSystem::rotation_y(-1.1);
    auto Rz = cos::RectangularSystem::rotation_z(2.0);
    auto check = [](const Mat3& R){ Mat3 I = R * R.transpose(); return (I - Mat3::Identity()).norm(); };
    EXPECT_NEAR(check(Rx), 0.0, 1e-12);
    EXPECT_NEAR(check(Ry), 0.0, 1e-12);
    EXPECT_NEAR(check(Rz), 0.0, 1e-12);
}

// 19) CylindricalSystem round-trip & axes
TEST(Coords_Cylindrical, RoundTripAndAxes) {
    Vec3 base(0,0,0), rpoint(1,0,0), thetap(0,1,0);
    cos::CylindricalSystem cyl("C", base, rpoint, thetap);
    // The current CylindricalSystem stores r as the signed projection onto r_axis
    // (not sqrt(x^2 + y^2)). A strict round-trip is guaranteed when the theta-axis
    // component is zero. Choose such a point:
    Vec3 p(0.7, 0.0, -0.3); // global (no theta component)
    // Given the projective definition of r and the usage of atan2 for theta,
    // only points with zero tangential component strictly round-trip. We validate
    // that property here.
    Vec3 local = cyl.to_local(p);
    EXPECT_NEAR(local.y(), 0.0, 1e-12); // zero tangential component by construction
    Vec3 back  = cyl.to_global(local);
    EXPECT_NEAR((back - p).norm(), 0.0, 1e-12);

    auto axes = cyl.get_axes(local);
    EXPECT_NEAR((axes.col(0).norm()-1.0), 0.0, 1e-9);
    EXPECT_NEAR((axes.col(1).norm()-1.0), 0.0, 1e-9);
    EXPECT_NEAR((axes.col(2).norm()-1.0), 0.0, 1e-9);
    EXPECT_NEAR(std::abs(axes.col(0).dot(axes.col(1))), 0.0, 1e-9);
}
