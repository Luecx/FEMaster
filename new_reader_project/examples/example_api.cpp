#include <iostream>
#include <array>
#include <fstream>  // <-- neu: wir schreiben demo.inp
#include "fem/reader2/registry.h"
#include "fem/reader2/condition_logic.h"
#include "fem/reader2/condition_parent.h"
#include "fem/reader2/condition_keys.h"
#include "fem/reader2/engine.h"

using namespace fem::reader2;

int main(){
    Registry reg;

    // ELASTIC only under MATERIAL
    reg.command("ELASTIC", [&](Command& c){
        c.allow_if( AND({ std::make_shared<ParentIs>(std::set<std::string>{"MATERIAL"}) }) );
        c.variant( Variant::make()
                      .when(nullptr)
                      .segment( Segment::make()
                                   .range( LineRange{}.min(1).max(1) )
                                   .pattern( Pattern::make().fixed<double,1>().fixed<double,1>() )
                                   .bind([](double E, double nu){
                                       std::cout << "[ELASTIC] E="<<E<<" nu="<<nu<<"\n";
                                   })
                                   )
        );
    });

    // NUMEIGENVALUES only if parent LOADCASE with the right TYPE
    reg.command("NUMEIGENVALUES", [&](Command& c){
        c.allow_if( AND({ std::make_shared<ParentIs>(std::set<std::string>{"LOADCASE"}),
                        std::make_shared<ParentKeyEquals>("TYPE", std::set<std::string>{"LINEARBUCKLING","EIGENFREQ"}) }));
        c.variant( Variant::make()
                      .when(nullptr)
                      .segment( Segment::make()
                                   .range( LineRange{}.min(1).max(1) )
                                   .pattern( Pattern::make().fixed<int,1>() )
                                   .bind([](int n){
                                       std::cout << "[NUMEIGENVALUES] " << n << "\n";
                                   })
                                   )
        );
    });

    // FOO variant with two segments
    reg.command("FOO", [&](Command& c){
        c.allow_if( OR({ std::make_shared<ParentIs>(std::set<std::string>{"ROOT"}),
                       std::make_shared<ParentIs>(std::set<std::string>{"BAR"}) }));
        c.variant( Variant::make()
                      .when( OR({ std::make_shared<KeyEquals>("MODE", std::set<std::string>{"A","B"}) }) )
                      .segment( Segment::make()
                                   .range( LineRange{}.min(1).max(1) )
                                   .pattern( Pattern::make().fixed<int,4>() )
                                   .bind([](std::array<int,4> a){
                                       std::cout << "[FOO-SEG1] ints="<<a[0]<<","<<a[1]<<","<<a[2]<<","<<a[3]<<"\n";
                                   })
                                   )
                      .segment( Segment::make()
                                   .range( LineRange{}.min(1).max(1) )
                                   .pattern( Pattern::make().fixed<double,2>().fixed<double,2>() )
                                   .bind([](double a,double b,double c,double d){
                                       std::cout << "[FOO-SEG2] doubles="<<a<<","<<b<<","<<c<<","<<d<<"\n";
                                   })
                                   )
        );
    });

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
        out << "1.5,2.5,3.5,4.5\n";
    }

    fem::reader::File file(path);
    Engine engine(reg);
    engine.run(file);

    return 0;
}
