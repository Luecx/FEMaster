#pragma once

#include "../core/core.h"
#include "../data/region.h"
#include "bc.h"

namespace fem {

struct Load : public BoundaryCondition {
    using Ptr                                                      = std::shared_ptr<Load>;

    virtual void apply(model::ModelData& model_data, NodeData& bc) = 0;
};

struct CLoad : public Load {
    using Ptr = std::shared_ptr<CLoad>;

    Vec6                   values {NAN, NAN, NAN, NAN, NAN, NAN};
    model::NodeRegion::Ptr region;

    CLoad() = default;
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

struct DLoad : public Load {
    using Ptr = std::shared_ptr<DLoad>;

    Vec3                      values {NAN, NAN, NAN};
    model::SurfaceRegion::Ptr region;

    DLoad() = default;
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

struct PLoad : public Load {
    using Ptr = std::shared_ptr<PLoad>;

    Precision                 pressure {NAN};
    model::SurfaceRegion::Ptr region;

    PLoad() = default;
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

struct VLoad : public Load {
    using Ptr = std::shared_ptr<VLoad>;

    Vec3                      values {NAN, NAN, NAN};
    model::ElementRegion::Ptr region;

    VLoad() = default;
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

// struct TLoad : public Load {
//     using Ptr = std::shared_ptr<TLoad>;
//
//     model::Region<model::NodeRegion>::Ptr region;
//
//     TLoad() = default;
//     void apply(model::ModelData& model_data, NodeData& bc) override {}
// }
}    // namespace fem
