/**
 * @file types_cls.h
 * @brief Forward declarations and shared-pointer aliases for key classes.
 *
 * Centralises forward declarations to reduce header dependencies between major
 * subsystems such as the model, coordinate systems, and materials.
 *
 * @see src/model/model.h
 * @see src/cos/coordinate_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include <memory>

namespace fem {
namespace model {

struct ModelData;
struct Model;
struct ElementInterface;
struct SurfaceInterface;

using ModelPtr = std::shared_ptr<Model>;
using ModelDataPtr = std::shared_ptr<ModelData>;
using ElementPtr = std::shared_ptr<ElementInterface>;
using SurfacePtr = std::shared_ptr<SurfaceInterface>;

} // namespace model

namespace cos {

struct CoordinateSystem;
using CoordinateSystemPtr = std::shared_ptr<CoordinateSystem>;

} // namespace cos

namespace material {

struct Material;
using MaterialPtr = std::shared_ptr<Material>;

} // namespace material

namespace loadcase {

struct LoadCase;
struct LinearStatic;
struct LinearStaticTopo;
struct LinearEigenfrequency;
struct LinearBuckling;

} // namespace loadcase

struct Section;
struct Profile;
struct PointMassSection;
struct SolidSection;
struct BeamSection;
struct ShellSection;

template<class T> using SPtr = std::shared_ptr<T>;
template<class T> using UPtr = std::unique_ptr<T>;
template<class T> using WPtr = std::weak_ptr<T>;

} // namespace fem
