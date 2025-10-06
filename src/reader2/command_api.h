#pragma once
#include "registry.h"
#include "schema.h"
#include "context.h"
#include "types.h"
#include "keys.h"

// This header is what command files should include.
// Example usage lives in your separate commands/*.cpp files via CommandRegistrar.

//
///**
// * @file main.cpp
// * @brief Demo entry: registers commands and parses an example input.
// */
//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include "parser/command_api.h"
//#include "parser/reader.h"
//#include "parser/keyword.h" // for Keyword getters
//
//using namespace fem::reader2;
//
//// --- Command registrations that print directly ---
//
//static CommandRegistrar REG_NODE("ROOT", "NODE", [](CommandSpec& c){
//    c.lines({1,kInf},
//        Schema::make_from_columns(
//            Schema::cols<int,double,double,double>()
//                .title("Node coordinates")
//                .cols_named({"ID","X","Y","Z"})
//                .cols_meaning({"Node id","X position","Y position","Z position"})
//                .bind([](Context&, int id,double x,double y,double z,const LineMeta&){
//                    std::cout << "[NODE] id=" << id << "  (" << x << "," << y << "," << z << ")\n";
//                })
//        )
//    );
//});
//
//static CommandRegistrar REG_ELEM("ROOT", "ELEMENT", [](CommandSpec& c){
//    c.lines({1,kInf},
//        Schema::make_from_columns(
//            Schema::cols<int,int,int,int,int,int,int,int,int>()
//                .title("Element connectivity (C3D8 demo)")
//                .cols_named({"ID","N1","N2","N3","N4","N5","N6","N7","N8"})
//                .bind([](Context&,int id,int n1,int n2,int n3,int n4,int n5,int n6,int n7,int n8,const LineMeta&){
//                    std::cout << "[ELEMENT] id=" << id << "  nodes="
//                              << n1<<" "<<n2<<" "<<n3<<" "<<n4<<" "
//                              << n5<<" "<<n6<<" "<<n7<<" "<<n8 << "\n";
//                })
//        )
//    );
//});
//
//static CommandRegistrar REG_MAT("ROOT", "MATERIAL", [](CommandSpec& c){
//    c.opens_scope("MATERIAL")
//     .doc("Material Block",
//          "Defines a material and its properties. Child commands specify constitutive data.")
//     .on_enter([](Context&, const Keyword& kw){
//         auto it = kw.kv().find("NAME");
//         std::cout << "[MATERIAL] begin '"
//                   << (it != kw.kv().end() ? it->second : std::string("UNNAMED"))
//                   << "'\n";
//     })
//     .on_exit([](Context&){ std::cout << "[MATERIAL] end\n"; });
//});
//
//static CommandRegistrar REG_ELASTIC("MATERIAL","ELASTIC",[](CommandSpec& c){
//    c.doc("Linear Elastic","Young's modulus and Poisson's ratio for isotropic material.")
//     .lines({1,1},
//        Schema::make_from_columns(
//            Schema::cols<double,double>() // E, nu
//                .title("E, nu")
//                .cols_named({"E","NU"})
//                .cols_meaning({"Young's modulus","Poisson's ratio"})
//                .cols_units({"Pa","-"})
//                .bind([](Context&, double E, double nu, const LineMeta&){
//                    std::cout << "  [ELASTIC] E=" << E << "  nu=" << nu << "\n";
//                })
//        )
//    );
//});
//
//static CommandRegistrar REG_DENSITY("MATERIAL","DENSITY",[](CommandSpec& c){
//    c.doc("Mass Density","Material mass density used for inertia and loads.")
//     .lines({1,1},
//        Schema::make_from_columns(
//            Schema::cols<double>() // rho
//                .title("rho")
//                .cols_named({"RHO"})
//                .cols_meaning({"Mass density"})
//                .cols_units({"kg/m^3"})
//                .bind([](Context&, double rho, const LineMeta&){
//                    std::cout << "  [DENSITY] rho=" << rho << "\n";
//                })
//        )
//    );
//});
//
//static CommandRegistrar REG_LOADCASE("ROOT","LOADCASE",[](CommandSpec& c){
//    c.opens_scope("LOADCASE")
//     .doc("Loadcase","A group of boundary conditions and loads, solved together.")
//     .on_enter([](Context&, const Keyword& kw){
//         std::cout << "[LOADCASE] begin ("
//                   << kw.get<std::string>("TYPE", "LINEAR STATIC") << ")\n";
//     })
//     .on_exit([](Context&){ std::cout << "[LOADCASE] end\n"; });
//});
//
//static void print_support(int id, std::initializer_list<double> vals){
//    std::cout << "  [SUPPORT] node=" << id << "  ";
//    int i=0; for (auto v: vals) std::cout << "dof" << i++ << "=" << v << "  ";
//    std::cout << "\n";
//}
//
//// SUPPORT is under ROOT
//static CommandRegistrar REG_SUPPORT("ROOT","SUPPORT",[](CommandSpec& c){
//    c.doc("Nodal Supports","Fix or restrain translational/rotational DOFs per node.")
//     .lines({1,kInf}, Schema::one_of({
//        Schema::make_from_columns(
//            Schema::cols<int,double,double,double,double,double,double>() // full
//                .na_floats_as_nan(true)
//                .title("Full (6 DOF)")
//                .cols_named({"ID","UX","UY","UZ","RX","RY","RZ"})
//                .cols_meaning({"Node id","X disp","Y disp","Z disp","X rot","Y rot","Z rot"})
//                .bind([](Context&, int id,double x,double y,double z,double rx,double ry,double rz, const LineMeta&){
//                    print_support(id,{x,y,z,rx,ry,rz});
//                })
//        ),
//        Schema::make_from_columns(
//            Schema::cols<int,double,double>() // short form
//                .na_floats_as_nan(true)
//                .title("Short (2 DOF)")
//                .cols_named({"ID","UX","UY"})
//                .cols_meaning({"Node id","X disp","Y disp"})
//                .bind([](Context&, int id,double x,double y, const LineMeta&){
//                    print_support(id,{x,y,NAN,NAN,NAN,NAN});
//                })
//        )
//    }));
//});
//
//static CommandRegistrar REG_CSYS("ROOT","COORDINATESYSTEM",[](CommandSpec& c){
//    c.doc("Coordinate System","Define a named coordinate system.")
//     .keys( KeyRules{}
//        .require("TYPE").allowed("TYPE", {"RECTANGULAR","CYLINDRICAL","SPHERICAL"})
//        .require("NAME")
//        .describe("TYPE","System type (orientation implied by implementation).")
//        .describe("NAME","Identifier for later reference.")
//     )
//     .lines({1,1}, Schema::make_from_columns(
//        Schema::cols<double,double,double>() // origin x,y,z
//            .title("Origin")
//            .cols_named({"OX","OY","OZ"})
//            .cols_meaning({"Origin X","Origin Y","Origin Z"})
//            .bind([](Context&, double ox,double oy,double oz, const LineMeta&){
//                std::cout << "[CSYS] origin=("<<ox<<","<<oy<<","<<oz<<")\n";
//            })
//    ));
//});
//
//
//// --- Example input ---
//static const char* EXAMPLE = R"(
//*NODE
//1, 0, 0, 0
//2, 1, 0, 0
//3, 1, 1, 0
//*ELEMENT, TYPE=C3D8
//1, 1,2,3,1,2,3,2,3
//*MATERIAL, NAME=STEEL
//*ELASTIC
//210000, 0.3
//*DENSITY
//7800
//*LOADCASE, TYPE=LINEAR STATIC
//*COORDINATESYSTEM, TYPE=RECTANGULAR, NAME=CS0
//0,0,0
//)";
//
//int main() {
//    Registry::instance().print_tree();
//
//    std::cout << Registry::instance().document_all() << std::endl;
//
//    // // write example file
//    // std::ofstream("example.inp") << EXAMPLE;
//    //
//    // Reader reader{};
//    // try {
//    //     reader.run("example.inp");
//    // } catch (const std::exception& e) {
//    //     std::cerr << "[ERROR] " << e.what() << std::endl;
//    //     return 1;
//    // }
//    //
//    // std::cout << "\nParsing complete.\n";
//}
