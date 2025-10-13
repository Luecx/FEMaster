// ===== FILE: ./examples/example_api.cpp =====

#include <iostream>
#include <array>
#include <fstream>
#include <string>

#include "fem/dsl/registry.h"
#include "fem/dsl/condition.h"
#include "fem/dsl/engine.h"
#include "fem/dsl/file.h"

using namespace fem::dsl;

int main() {
    Registry reg;

    // *EXAMPLE under ROOT with 4 variants:
    //   V1: 3 integers (single line)
    //   V2: 3 doubles  (single line)
    //   V3: 1 string + 1 integer (single line)
    //   V4: 5 floats + 1 string (multiline allowed, floats can be on_missing-filled)
    reg.command("EXAMPLE", [&](Command& c) {
        c.allow_if( Condition::parent_is("ROOT") );
		c.variant( Variant::make()
            .doc("Two integers: i1, i2.")
            .segment( Segment::make()
                .range( LineRange{}.min(1).max(10) )
                .pattern( Pattern::make()
                    .fixed<int,2>().name("i").desc("integer inputs") )
                .bind([](std::array<int,2> i) {
                    std::cout << "[EXAMPLE/V0:int,int] "
                              << i[0] << ", " << i[1] << ", "  "\n";
                })
            )
        );
        // Variant 1 — three integers
        c.variant( Variant::make()
            .doc("Three integers: i1, i2, i3.")
            .segment( Segment::make()
                .range( LineRange{}.min(1).max(10) )
                .pattern( Pattern::make()
                    .fixed<int,3>().name("i").desc("integer inputs") )
                .bind([](std::array<int,3> i) {
                    std::cout << "[EXAMPLE/V1:int,int,int] "
                              << i[0] << ", " << i[1] << ", " << i[2] << "\n";
                })
            )
        );

        // Variant 2 — three doubles
        c.variant( Variant::make()
            .doc("Three floats: x1..x3 as doubles.")
            .segment( Segment::make()
                .range( LineRange{}.min(1).max(10) )
                .pattern( Pattern::make()
                    .fixed<double,3>().name("x").desc("floating-point inputs") )
                .bind([](std::array<double,3> x) {
                    std::cout << "[EXAMPLE/V2:double,double,double] "
                              << x[0] << ", " << x[1] << ", " << x[2] << "\n";
                })
            )
        );

        // Variant 3 — one string then one integer
        c.variant( Variant::make()
            .doc("One string label and one integer value.")
            .segment( Segment::make()
                .range( LineRange{}.min(1).max(10) )
                .pattern( Pattern::make()
                    .fixed<std::string,1>().name("label").desc("name/label")
                    .fixed<int,1>().name("n").desc("count") )
                .bind([](std::string label, int n) {
                    std::cout << "[EXAMPLE/V3:string,int] "
                              << label << ", " << n << "\n";
                })
            )
        );

         c.variant( Variant::make()
            .doc("One string label and one integer value.")
            .segment( Segment::make()
                .range( LineRange{}.min(1).max(10) )
                .pattern( Pattern::make()
                    .fixed<std::string,1>().name("label").desc("name/label")
                     )
                .bind([](std::string label) {
                    std::cout << "[EXAMPLE/V3:string] "
                              << label << "\n";
                })
            )
        );

        // Variant 4 — multiline: five floats + one string
        // Floats allow on_missing fill; string is required but can be on_missing defaulted too.
        c.variant( Variant::make()
            .doc("Five floats (multiline allowed) followed by one string.")
            .segment( Segment::make()
                .range( LineRange{}.min(1).max(10) )
                .pattern( Pattern::make()
                    .allow_multiline()
                    .fixed<float,5>().name("h").desc("measurements").on_missing(0.0f)
                    .fixed<std::string,1>().name("note").desc("comment").on_missing(std::string("default")) )
                .bind([](std::array<float,5> h, std::string note) {
                    std::cout << "[EXAMPLE/V4:floatx5,string] "
                              << h[0] << ", " << h[1] << ", " << h[2] << ", "
                              << h[3] << ", " << h[4] << " | " << note << "\n";
                })
            )
        );
    });

    // Print a compact reference so you can verify the variant docs & layouts.
    reg.print_help();

    // Build a tiny input file in-memory to exercise all four variants in order.
    // Each *EXAMPLE block is followed by data matching one specific variant.
    std::string path = "demo.inp";
    {
        std::ofstream out(path);

        // Variant 1: three integers
        out << "*EXAMPLE\n";
        out << "1,2,3\n";

        // Variant 2: three doubles
        out << "*EXAMPLE\n";
        out << "1.5,2.25,3.75\n";
        out << "1.5,2.25,3.75\n";

        // Variant 3: string + integer
        out << "*EXAMPLE\n";
        out << "steel,42\n";

        // Variant 4: multiline five floats + one string (exactly 6 tokens total)
        out << "*EXAMPLE\n";
        out << "10,20\n";        // first two floats
        out << "30,40,50\n";     // next three floats (now 5 floats total)
        out << "note-goes-here\n"; // final token: string
    }

    fem::dsl::File file(path);
    Engine engine(reg);
    engine.run(file);

    return 0;
}
