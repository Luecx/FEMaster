//
// Created by f_eggers on 21.11.2024.
//

#ifndef SECTION_H
#define SECTION_H

#include "../material/material.h"
#include "../data/region.h"

#include "profile.h"

namespace fem {

struct Section {
    using Ptr = std::shared_ptr<Section>;

    virtual ~Section() = default; // Make Section polymorphic

    material::Material::Ptr material = nullptr;
    model::ElementRegion::Ptr region = nullptr;

    template<typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }
    virtual void info() = 0;
};

struct SolidSection : Section {
    using Ptr = std::shared_ptr<SolidSection>;
    void info() {
        logging::info(true, "Section: ");
        logging::info(true, "   Material: ", (material ? material->name : "-"));
        logging::info(true, "   Region  : ", region->name);
    }
};

struct BeamSection : Section{
    using Ptr = std::shared_ptr<BeamSection>;
    Vec3 n1;
    Profile::Ptr profile = nullptr;

    void info() {
        logging::info(true, "BeamSection: ");
        logging::info(true, "   Material: ", (material ? material->name : "-"));
        logging::info(true, "   Region  : ", region->name);
        logging::info(true, "   Profile : ", (profile ? profile->name : "-"));
        logging::info(true, "   n1      : ", n1.transpose());
    }
};

struct ShellSection : Section {
    using Ptr = std::shared_ptr<ShellSection>;
    Precision thickness = 1.0;

    void info() {
        logging::info(true, "ShellSection: ");
        logging::info(true, "   Material: ", (material ? material->name : "-"));
        logging::info(true, "   Region  : ", (region ? region->name : "-"));
        logging::info(true, "   Thickness: ", thickness);
    }
};

struct PointMassSection : Section {
    using Ptr = std::shared_ptr<PointMassSection>;
    Precision mass = 0;
    Vec3 rotary_inertia = Vec3::Zero();
    Vec3 spring_constants = Vec3::Zero();
    Vec3 rotary_spring_constants  = Vec3::Zero();
    void info() {
        logging::info(true, "PointMassSection: ");
        logging::info(true, "   Material : ", (material ? material->name : "-"));
        logging::info(true, "   Region   : ", region->name);
        logging::info(true, "   Mass     : ", mass);
        logging::info(true, "   Inertia  : ", rotary_inertia.transpose());
        logging::info(true, "   Springs  : ", spring_constants.transpose());
        logging::info(true, "   Rotations: ", rotary_spring_constants.transpose());
    }
};

}


#endif //SECTION_H
