#pragma once

#include "../constraints/equation.h"
#include "../core/core.h"
#include "../cos/coordinate_system.h"
#include "../data/region.h"

namespace fem {

template<RegionTypes RT>
struct Support {

    typename model::Region<RT>::Ptr region = nullptr;
    Vec6 values{NAN, NAN, NAN, NAN, NAN, NAN};
    cos::CoordinateSystem::Ptr coordinate_system = nullptr;

    virtual void apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations) = 0;
};

struct NSupport : public Support<RegionTypes::NODE> {
    void apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations) override;
};

struct ESupport : public Support<RegionTypes::ELEMENT> {
    void apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations) override;
};

struct SSupport : public Support<RegionTypes::SURFACE> {
    void apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations) override;
};


} // namespace fem


