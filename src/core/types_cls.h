#pragma once

#include <memory>

namespace fem::model {
    struct ModelData;
    struct Model;
    struct ElementInterface;
    struct SurfaceInterface;

    using ModelPtr = std::shared_ptr<Model>;
    using ModelDataPtr = std::shared_ptr<ModelData>;
    using ElementPtr = std::shared_ptr<ElementInterface>;
    using SurfacePtr = std::shared_ptr<SurfaceInterface>;
}

namespace fem::cos {
    struct CoordinateSystem;
    using CoordinateSystemPtr = std::shared_ptr<CoordinateSystem>;
}

namespace fem::material {
    struct Material;
    using MaterialPtr = std::shared_ptr<Material>;
}

namespace fem::loadcase {
    struct Loadcase;
    struct LinearStatic;
    struct LinearStaticTopo;
    struct LinearEigenfrequency;
}

namespace fem {
    struct Section;
    struct Profile;
    struct PointMassSection;
    struct SolidSection;
    struct BeamSection;
    struct ShellSection;
}