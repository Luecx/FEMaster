
#include "core/core.h"
#include "cuda/cuda.h"
#include "loadcase/linear_static.h"
#include "material/material.h"
#include "math/quadrature.h"
#include "model/c3d10.h"
#include "model/c3d15.h"
#include "model/c3d20.h"
#include "model/c3d4.h"
#include "model/c3d6.h"
#include "model/c3d8.h"
#include "model/element.h"
#include "model/model.h"
#include "reader/file.h"
#include "reader/line.h"
#include "reader/reader.h"
#include "solve/solver.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>

int main(int argc, char* argv[]) {
//    auto quad1 = fem::quadrature::Quadrature(fem::quadrature::DOMAIN_ISO_TET  , fem::quadrature::ORDER_LINEAR);
//    auto quad2 = fem::quadrature::Quadrature(fem::quadrature::DOMAIN_ISO_TET  , fem::quadrature::ORDER_QUADRATIC);
//    auto quad3 = fem::quadrature::Quadrature(fem::quadrature::DOMAIN_ISO_TET  , fem::quadrature::ORDER_CUBIC);
//    auto quad4 = fem::quadrature::Quadrature(fem::quadrature::DOMAIN_ISO_WEDGE, fem::quadrature::ORDER_SUPER_LINEAR);
//    auto quad5 = fem::quadrature::Quadrature(fem::quadrature::DOMAIN_ISO_WEDGE, fem::quadrature::ORDER_SUPER_QUADRATIC);

    // Define a cubic function wrapped in std::function
//    std::function<double(double, double, double)> f_cubic =
//        [](double x, double y, double z) {
//            return x*(x-1)*y - x*x*x + 2 * x * z -y*y;
//        };
//    std::function<double(double, double, double)> f_quad =
//        [](double x, double y, double z) {
//            return (1-x) * y + y*(3-y) + z;
//        };
//    std::function<double(double, double, double)> f_linear =
//        [](double x, double y, double z) {
//            return (3-x) + y * y + z * z * z * 3;
//        };
//
//    std::cout << quad4.integrate(f_linear) << std::endl;
//    std::cout << quad5.integrate(f_linear) << std::endl;

//    std::cout << quad1.integrate(f_cubic) << std::endl;
//    std::cout << quad2.integrate(f_cubic) << std::endl;
//    std::cout << quad3.integrate(f_cubic) << std::endl;
//
//    std::cout << quad1.integrate(f_quad) << std::endl;
//    std::cout << quad2.integrate(f_quad) << std::endl;
//    std::cout << quad3.integrate(f_quad) << std::endl;

//    std::cout << quad1.integrate(f_linear) << std::endl;
//    std::cout << quad2.integrate(f_linear) << std::endl;
//    std::cout << quad3.integrate(f_linear) << std::endl;

//    fem::model::SolidElement< 4>::test_implementation<fem::model::C3D4 >();
//    fem::model::SolidElement< 6>::test_implementation<fem::model::C3D6 >();
//    fem::model::SolidElement< 8>::test_implementation<fem::model::C3D8 >();
//    fem::model::SolidElement<10>::test_implementation<fem::model::C3D10>();
//    fem::model::SolidElement<15>::test_implementation<fem::model::C3D15>();
//    fem::model::SolidElement<20>::test_implementation<fem::model::C3D20>();

    fem::reader::Reader reader{std::string(argv[1])};
    reader.read();
//    auto c3d8 = fem::model::C3D8(0, {0,1,2,3,4,5,6,7});
//    auto local = c3d8.node_coords_local();
//    auto global = StaticMatrix<8,3>{};
//    global << 0, 0, 0,
//              1, 0, 0,
//              1, 1, 0,
//              0, 1, 0,
//              0, 0, 1,
//              1, 0, 1,
//              1, 1, 1,
//              0, 1, 1;
//
//    std::cout << local << std::endl;
//    std::cout << global << std::endl;
//
//    Precision det = 0;
//    std::cout << c3d8.strain_displacements(global, local(1,0), local(1,1), local(1,2), det) << std::endl;

    return 0;
}
