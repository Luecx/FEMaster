#pragma once

#include "../core/core.h"
#include <memory>

namespace fem {
namespace material {
struct Elasticity {
    public:
    virtual StaticMatrix<3, 3> get_2d() = 0;
    virtual StaticMatrix<6, 6> get_3d() = 0;

    virtual ~Elasticity() = default;

    template<size_t D>
    StaticMatrix<D == 3 ? 6:3,D == 3 ? 6:3> get() {
        if constexpr (D == 2){
            return get_2d();
        }
        if constexpr (D == 3){
            return get_3d();
        }
        throw std::out_of_range("indexing _material elasticity outside of valid dimensions which is 2 or 3");
    }
};

using ElasticityPtr = std::shared_ptr<Elasticity>;
}    // namespace _material
}    // namespace fem
