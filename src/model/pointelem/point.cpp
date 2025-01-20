//
// Created by f_eggers on 11.12.2024.
//

#include "point.h"
#include "../../section/section_point_mass.h"

fem::MapMatrix fem::model::Point::stiffness(Precision* buffer) {
    MapMatrix(buffer, 6, 6).setZero();
    logging::error(_section != nullptr, "Section not assigned to element ", elem_id);
    logging::error(_section->as<PointMassSection>(), "Section is not a PointMassSection");
    PointMassSection* section = _section->as<PointMassSection>();
    MapMatrix         res     = MapMatrix(buffer, 6, 6).setZero();
    res(0, 0)                 = section->spring_constants(0);
    res(1, 1)                 = section->spring_constants(1);
    res(2, 2)                 = section->spring_constants(2);
    res(3, 3)                 = section->rotary_spring_constants(0);
    res(4, 4)                 = section->rotary_spring_constants(1);
    res(5, 5)                 = section->rotary_spring_constants(2);
    return res;
}
fem::MapMatrix fem::model::Point::mass(Precision* buffer) {
    logging::error(_section != nullptr, "Section not assigned to element ", elem_id);
    logging::error(_section->as<PointMassSection>(), "Section is not a PointMassSection");
    PointMassSection* section = _section->as<PointMassSection>();
    MapMatrix         res     = MapMatrix(buffer, 6, 6).setZero();
    res(0, 0)                 = section->mass;
    res(1, 1)                 = section->mass;
    res(2, 2)                 = section->mass;
    res(3, 3)                 = section->rotary_inertia(0);
    res(4, 4)                 = section->rotary_inertia(1);
    res(5, 5)                 = section->rotary_inertia(2);
    return res;
}