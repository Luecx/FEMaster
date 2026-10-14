// Times only prepare_pretension_sections: no parsing, compile or solver in timer.
#include "model/model.h"
#include "model/solid/c3d8.h"
#include <chrono>
#include <iostream>
#include <iomanip>
using namespace fem;
int main() {
    for (int n : {4, 8, 16, 32, 64}) {
        for (int repeat = 0; repeat < 6; ++repeat) {
            model::Model model;
            auto node = [n](int x, int y, int z) { return z*(n+1)*(n+1)+y*(n+1)+x; };
            for (int z=0; z<4; ++z) for(int y=0; y<=n; ++y) for(int x=0; x<=n; ++x)
                model.set_node(node(x,y,z), double(x)/n, double(y)/n, z);
            model._data->parts.get()->surface_sets.activate("CUT");
            int e=0, s=0;
            for(int z=0; z<3; ++z) for(int y=0; y<n; ++y) for(int x=0; x<n; ++x) {
                model.set_element<model::C3D8>(e,
                    node(x,y,z),node(x+1,y,z),node(x+1,y+1,z),node(x,y+1,z),
                    node(x,y,z+1),node(x+1,y,z+1),node(x+1,y+1,z+1),node(x,y+1,z+1));
                if (z==1) model.set_surface(s++,e,1);
                ++e;
            }
            model.compile();
            model.add_pretension_interface_section("PT", "CUT", "CUT");
            const auto start=std::chrono::steady_clock::now();
            model.prepare_pretension_sections();
            const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
            std::cout << "BENCH," << n << ',' << repeat << ',' << 4*(n+1)*(n+1) << ',' << e << ','
                      << s << ',' << (n+1)*(n+1) << ',' << model._data->positions->rows << ','
                      << std::setprecision(10) << ms << '\n';
        }
    }
}
