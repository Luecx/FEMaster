// Isolated unique reference-coordinate matching kernel from model_pretension.cpp.
#include "core/types_eig.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <set>
#include <map>
#include <vector>
using namespace fem;
int main() {
    std::cout << "n,repeat,matching_ms,pairs\n";
    for (int n : {4,8,16,32,64}) for (int repeat=0; repeat<6; ++repeat) {
        std::vector<Vec3> positions;
        std::set<ID> nodes_a, nodes_b;
        for (int y=0; y<=n; ++y) for (int x=0; x<=n; ++x) {
            nodes_a.insert(static_cast<ID>(positions.size()));
            positions.emplace_back(Precision(x)/n, Precision(y)/n, 1);
        }
        for (int y=0; y<=n; ++y) for (int x=0; x<=n; ++x) {
            nodes_b.insert(static_cast<ID>(positions.size()));
            positions.emplace_back(Precision(x)/n, Precision(y)/n, 1);
        }
        const auto start = std::chrono::steady_clock::now();
        std::map<ID,ID> pairs;
        std::set<ID> used_b;
        for (const ID a : nodes_a) {
            const Vec3 point = positions[a];
            ID match = -1;
            for (const ID b : nodes_b) {
                if ((point-positions[b]).norm() > Precision(1e-8)) continue;
                if (match >= 0) return 1;
                match = b;
            }
            if (match < 0 || !used_b.insert(match).second) return 1;
            pairs.emplace(a,match);
        }
        const double ms = std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-start).count();
        std::cout << n << ',' << repeat << ',' << std::setprecision(10) << ms << ',' << pairs.size() << '\n';
    }
}
