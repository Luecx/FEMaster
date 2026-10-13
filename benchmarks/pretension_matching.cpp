// Isolated matching kernel from model_pretension.cpp; excludes all topology work.
#include "core/types_eig.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <set>
#include <map>
#include <vector>
#include <limits>
using namespace fem;
int main() {
    std::cout << "n,repeat,matching_ms,pairs\n";
    for (int n : {4,8,16,32,64}) for (int repeat=0; repeat<6; ++repeat) {
        std::vector<Vec3> positions;
        std::set<ID> nodes;
        for (int y=0; y<=n; ++y) for (int x=0; x<=n; ++x) {
            nodes.insert(static_cast<ID>(positions.size()));
            positions.emplace_back(Precision(x)/n, Precision(y)/n, 1);
        }
        const auto start = std::chrono::steady_clock::now();
        std::map<ID,ID> pairs;
        auto unused = nodes;
        for (const ID a : nodes) {
            const Vec3 point = positions[a];
            ID best = -1;
            Precision distance_best = std::numeric_limits<Precision>::max();
            for (const ID b : unused) {
                const Precision distance = (point-positions[b]).norm();
                if (distance < distance_best) { distance_best=distance; best=b; }
            }
            pairs.emplace(a,best);
            unused.erase(best);
        }
        const double ms = std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-start).count();
        std::cout << n << ',' << repeat << ',' << std::setprecision(10) << ms << ',' << pairs.size() << '\n';
    }
}
