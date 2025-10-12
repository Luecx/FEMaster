// ===== FILE: ./examples/example_api.cpp =====

#include <iostream>
#include <array>
#include <fstream>

#include "fem/dsl/registry.h"
#include "fem/dsl/condition.h"
#include "fem/dsl/engine.h"
#include "fem/dsl/file.h"

using namespace fem::dsl;

int main() {
    Registry reg;

    // *ELASTIC only under *MATERIAL
    reg.command("ELASTIC", [&](Command& c) {
        c.allow_if( Condition::parent_is("MATERIAL") );
        c.variant( Variant::make()
            .segment( Segment::make()
                .range( LineRange{}.min(1).max(1) )
                .pattern( Pattern::make().fixed<double,1>().fixed<double,1>() )
                .bind([](double E, double nu) {
                    std::cout << "[ELASTIC] E=" << E << " nu=" << nu << "\n";
                })
            )
        );
    });

    // *NUMEIGENVALUES only if parent *LOADCASE with TYPE in {LINEARBUCKLING,EIGENFREQ}
    reg.command("NUMEIGENVALUES", [&](Command& c) {
        c.allow_if( Condition::all_of({
            Condition::parent_is("LOADCASE"),
            Condition::parent_key_equals("TYPE", {"LINEARBUCKLING", "EIGENFREQ"})
        }) );
        c.variant( Variant::make()
            .segment( Segment::make()
                .range( LineRange{}.min(1).max(1) )
                .pattern( Pattern::make().fixed<int,1>() )
                .bind([](int n) {
                    std::cout << "[NUMEIGENVALUES] " << n << "\n";
                })
            )
        );
    });

    // *FOO variant with two segments; first requires MODE in {A,B}
    reg.command("FOO", [&](Command& c) {
        c.allow_if( Condition::any_of({
            Condition::parent_is("ROOT"),
            Condition::parent_is("BAR")
        }) );

        c.variant( Variant::make()
            .when( Condition::key_equals("MODE", {"A", "B"}) )
            .segment( Segment::make()
                .range( LineRange{}.min(1) )
                .pattern( Pattern::make().allow_multiline()
                        .fixed<float,10>().name("h").desc("some text").on_missing(0.0f)
                        .fixed<float,1 >().name("g").desc("some other text").on_missing(123123.0f) )
                .bind([](std::array<float,10> a, float g) {
                    std::cout << "[FOO-SEG1] ints="
                              << a[0] << "," << a[1] << "," << a[2] << "," << a[3] <<
                                    "," << a[4] << "," << a[5] << "," << a[6] << "," << a[7] <<
                                    "," << a[8] << "," << a[9] << "," << g << "\n";
                })
            )
        );
    });
	reg.print_help();
    // Build a tiny input file in-memory and run
    std::string path = "demo.inp";
    {
        std::ofstream out(path);
        out << "*MATERIAL,NAME=STEEL\n";
        out << "*ELASTIC\n";
        out << "210000,0.3\n";
        out << "*LOADCASE,TYPE=EIGENFREQ\n";
        out << "*NUMEIGENVALUES\n";
        out << "12\n";
        out << "*FOO,MODE=A\n";
        out << "1,2,3,4\n";
        out << "1,2,3,4\n";
        out << "1,2\n";
        out << "1,2,3,4\n";
        out << "1,2,3,4\n";
        out << "1,2\n";
        out << "1,2,3,4\n";
    }

    fem::dsl::File file(path);
    Engine engine(reg);
    engine.run(file);

    return 0;
}
