#pragma once

#include "../core/core.h"
#include "../data/region.h"
#include "../data/node_data_dict.h"
#include "bc.h"

namespace fem {

struct Load : public BoundaryCondition {
    using Ptr = std::shared_ptr<Load>;

    virtual ~Load() = default;

    virtual void apply(model::ModelData& model_data, NodeData& bc) = 0;
};

struct CLoad : public Load {
    using Ptr = std::shared_ptr<CLoad>;

    virtual ~CLoad() = default;

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
    virtual ~DLoad() = default;
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

struct PLoad : public Load {
    using Ptr = std::shared_ptr<PLoad>;

    Precision                 pressure {NAN};
    model::SurfaceRegion::Ptr region;

    PLoad() = default;
    virtual ~PLoad() = default;
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

struct VLoad : public Load {
    using Ptr = std::shared_ptr<VLoad>;

    Vec3                      values {NAN, NAN, NAN};
    model::ElementRegion::Ptr region;

    VLoad() = default;
    virtual ~VLoad() = default;
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

struct TLoad : public Load {
    using Ptr = std::shared_ptr<TLoad>;

    model::NodeField::Ptr temp_field;
    Precision             ref_temp {NAN};


    TLoad() = default;
    virtual ~TLoad() = default;
    void apply(model::ModelData& model_data, NodeData& bc) override;
};


}    // namespace fem
