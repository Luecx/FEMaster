/**
 * @file section_shell_composite.cpp
 * @brief Implements layer-wise through-thickness shell material integration.
 *
 * Every layer is integrated with a composite Simpson rule over its own physical
 * thickness interval. Generalized strains are rotated into the layer axes once
 * per layer; the layer resultants and tangent are rotated back into the
 * geometric shell basis and accumulated.
 *
 * A null output-state pointer denotes a state-neutral constitutive query and is
 * forwarded to every material point without pointer arithmetic.
 *
 * @see CompositeShellSection
 *
 * @author Finn Eggers
 * @date 24.11.2026
 */

#include "section_shell_composite.h"

#include "../core/logging.h"
#include "../material/strain/shell_material_strain.h"
#include "../material/stress/shell_material_stress_cauchy.h"
#include "../material/stress/shell_material_stress_pk2.h"

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace fem {

/**
 * Constructs a layered shell section.
 *
 * The total thickness is the sum of all layer thicknesses. All layers must
 * share the thermal expansion coefficient and reference temperature because
 * the shell element maps free thermal expansion as a constant membrane strain.
 * All layers must also share one density because shell mass integration uses a
 * constant density through the thickness.
 *
 * @param layers Layers ordered from the bottom to the top shell surface.
 * @param region Element region receiving the section.
 * @param orientation Optional coordinate system defining the section basis.
 * @param csys_axis Zero-based coordinate-system axis projected into the shell plane.
 */
CompositeShellSection::CompositeShellSection(
    std::vector<Layer>         layers,
    model::ElementRegion::Ptr  region,
    cos::CoordinateSystem::Ptr orientation,
    Index                      csys_axis
)
    : ShellSection(
          nullptr,
          std::move(region),
          std::accumulate(layers.begin(), layers.end(), Precision(0),
              [](Precision sum, const Layer& layer) { return sum + layer.thickness; }),
          std::move(orientation),
          csys_axis
      ),
      layers_(std::move(layers)) {
    logging::error(!layers_.empty(),
        "CompositeShellSection: at least one layer is required");

    const auto& reference = layers_.front().material;

    for (const Layer& layer : layers_) {
        logging::error(layer.material && layer.material->has_elasticity(),
            "CompositeShellSection: every layer requires a material with elasticity");
        logging::error(layer.thickness > Precision(0),
            "CompositeShellSection: layer thickness must be positive");
        logging::error(layer.points >= 3 && layer.points % 2 == 1,
            "CompositeShellSection: layer integration points must be odd and at least 3");

        // One common free expansion law for the complete section
        const auto& material = layer.material;
        logging::error(material->has_thermal_expansion() == reference->has_thermal_expansion()
                    && (!material->has_thermal_expansion()
                        || (material->get_thermal_expansion()        == reference->get_thermal_expansion()
                         && material->get_thermal_zero_temperature() == reference->get_thermal_zero_temperature())),
            "CompositeShellSection: all layers must share the same thermal expansion");

        // Mass integration currently assumes one density through the thickness
        logging::error(material->has_density() == reference->has_density()
                    && (!material->has_density()
                        || material->get_density() == reference->get_density()),
            "CompositeShellSection: all layers must share the same density");

        num_mp_ += layer.points;
    }
}

/**
 * Integrates one shell integration point through all layers.
 *
 * Within a layer with `n` points the strain is sampled at equally spaced
 * coordinates `z_k` between its bottom and top surface and weighted with the
 * composite Simpson rule. The layer axes follow from the section axes by an
 * in-plane rotation by the layer angle.
 *
 * @param position_reference Physical reference position of the shell point.
 * @param shell_basis_global Geometric shell basis in global coordinates.
 * @param strain_shell Generalized strain in the geometric shell basis.
 * @param old_material_state First through-thickness input state row.
 * @param new_material_state Optional first through-thickness output state row.
 * @param material_state_stride Scalar distance between consecutive state rows.
 * @param resultants_shell Integrated resultants in the geometric shell basis.
 * @param tangent_shell Consistent generalized tangent in the geometric shell basis.
 */
void CompositeShellSection::evaluate(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_shell,
    const Precision*              old_material_state,
    Precision*                    new_material_state,
    Index                         material_state_stride,
    ShellStressResultants&        resultants_shell,
    Mat8&                         tangent_shell
) const {
    using GeneralizedStrainComponent = ShellGeneralizedStrain::Component;
    using MaterialStrainComponent    = ShellMaterialStrain::Component;
    using MaterialStressComponent    = ShellMaterialStress::Component;
    using ResultantComponent         = ShellStressResultants::Component;

    constexpr GeneralizedStrainComponent EXX = GeneralizedStrainComponent::EpsilonXX;
    constexpr GeneralizedStrainComponent EYY = GeneralizedStrainComponent::EpsilonYY;
    constexpr GeneralizedStrainComponent GXY = GeneralizedStrainComponent::GammaXY;
    constexpr GeneralizedStrainComponent KXX = GeneralizedStrainComponent::KappaXX;
    constexpr GeneralizedStrainComponent KYY = GeneralizedStrainComponent::KappaYY;
    constexpr GeneralizedStrainComponent KXY = GeneralizedStrainComponent::KappaXY;
    constexpr GeneralizedStrainComponent GXZ = GeneralizedStrainComponent::GammaXZ;
    constexpr GeneralizedStrainComponent GYZ = GeneralizedStrainComponent::GammaYZ;

    constexpr MaterialStrainComponent MXX = MaterialStrainComponent::XX;
    constexpr MaterialStrainComponent MYY = MaterialStrainComponent::YY;
    constexpr MaterialStrainComponent MXY = MaterialStrainComponent::GammaXY;
    constexpr MaterialStrainComponent MXZ = MaterialStrainComponent::GammaXZ;
    constexpr MaterialStrainComponent MYZ = MaterialStrainComponent::GammaYZ;

    constexpr MaterialStressComponent SXX = MaterialStressComponent::XX;
    constexpr MaterialStressComponent SYY = MaterialStressComponent::YY;
    constexpr MaterialStressComponent SXY = MaterialStressComponent::XY;
    constexpr MaterialStressComponent SXZ = MaterialStressComponent::XZ;
    constexpr MaterialStressComponent SYZ = MaterialStressComponent::YZ;

    constexpr ResultantComponent NXX = ResultantComponent::NXX;
    constexpr ResultantComponent NYY = ResultantComponent::NYY;
    constexpr ResultantComponent NXY = ResultantComponent::NXY;
    constexpr ResultantComponent BXX = ResultantComponent::MXX;
    constexpr ResultantComponent BYY = ResultantComponent::MYY;
    constexpr ResultantComponent BXY = ResultantComponent::MXY;
    constexpr ResultantComponent QX  = ResultantComponent::QX;
    constexpr ResultantComponent QY  = ResultantComponent::QY;

    // Section axes in shell coordinates. Without an orientation the geometric
    // shell basis is the section basis.
    const Mat3 section_basis_global = orientation_
        ? stress_basis(position_reference, shell_basis_global)
        : shell_basis_global;

    const Mat2 section_axes_in_shell =
        shell_basis_global.template block<3, 2>(0, 0).transpose()
        * section_basis_global.template block<3, 2>(0, 0);

    const Precision shear_correction = Precision(5) / Precision(6);
    const Precision degree           = Precision(std::acos(-1.0) / 180.0);

    resultants_shell.values().setZero();
    tangent_shell.setZero();

    Precision z_bottom = Precision(-0.5) * thickness_;
    Index     mp       = 0;

    for (const Layer& layer : layers_) {
        // Rotate the section axes about the normal into the layer axes.
        const Precision c = std::cos(layer.angle * degree);
        const Precision s = std::sin(layer.angle * degree);

        Mat2 layer_rotation;
        layer_rotation << c, -s,
                          s,  c;

        const Mat2 layer_axes_in_shell     = section_axes_in_shell * layer_rotation;
        const Mat8 strain_shell_to_layer   = ShellGeneralizedStrain::transformation(layer_axes_in_shell);
        const Mat8 resultants_layer_to_shell =
            ShellStressResultants::transformation(layer_axes_in_shell.transpose());

        const ShellGeneralizedStrain strain_layer(strain_shell_to_layer * strain_shell.values());

        ShellStressResultants resultants_layer;
        Mat8                  tangent_layer;

        resultants_layer.values().setZero();
        tangent_layer.setZero();

        const Precision dz = layer.thickness / Precision(layer.points - 1);

        for (Index point = 0; point < layer.points; ++point, ++mp) {
            // Composite Simpson weights 1, 4, 2, ..., 4, 1 times dz/3
            const Precision z       = z_bottom + Precision(point) * dz;
            const bool      end     = point == 0 || point == layer.points - 1;
            const Precision simpson = end ? Precision(1) : (point % 2 == 1 ? Precision(4) : Precision(2));
            const Precision w       = simpson * dz / Precision(3);

            ShellMaterialStrain material_strain;
            material_strain[MXX] = strain_layer[EXX] + z * strain_layer[KXX];
            material_strain[MYY] = strain_layer[EYY] + z * strain_layer[KYY];
            material_strain[MXY] = strain_layer[GXY] + z * strain_layer[KXY];
            material_strain[MXZ] = strain_layer[GXZ];
            material_strain[MYZ] = strain_layer[GYZ];

            const Precision* old_state = old_material_state + mp * material_state_stride;
            Precision* new_state = new_material_state
                ? new_material_state + mp * material_state_stride
                : nullptr;

            ShellMaterialStressPK2 material_stress_pk2;
            Mat5                   material_tangent;

            layer.material->elasticity()->evaluate(
                material_strain,
                old_state,
                new_state,
                material_stress_pk2,
                &material_tangent
            );

            const ShellMaterialStress material_stress(material_stress_pk2.values());

            resultants_layer[NXX] += w * material_stress[SXX];
            resultants_layer[NYY] += w * material_stress[SYY];
            resultants_layer[NXY] += w * material_stress[SXY];

            resultants_layer[BXX] += w * z * material_stress[SXX];
            resultants_layer[BYY] += w * z * material_stress[SYY];
            resultants_layer[BXY] += w * z * material_stress[SXY];

            resultants_layer[QX] += w * shear_correction * material_stress[SXZ];
            resultants_layer[QY] += w * shear_correction * material_stress[SYZ];

            // d epsilon_material = B_z d epsilon_generalized
            StaticMatrix<5, 8> strain_map = StaticMatrix<5, 8>::Zero();
            strain_map((Index)MXX, (Index)EXX) = Precision(1);
            strain_map((Index)MXX, (Index)KXX) = z;
            strain_map((Index)MYY, (Index)EYY) = Precision(1);
            strain_map((Index)MYY, (Index)KYY) = z;
            strain_map((Index)MXY, (Index)GXY) = Precision(1);
            strain_map((Index)MXY, (Index)KXY) = z;
            strain_map((Index)MXZ, (Index)GXZ) = Precision(1);
            strain_map((Index)MYZ, (Index)GYZ) = Precision(1);

            // d R = A_z d sigma_material
            StaticMatrix<8, 5> resultant_map = StaticMatrix<8, 5>::Zero();
            resultant_map((Index)NXX, (Index)SXX) = Precision(1);
            resultant_map((Index)NYY, (Index)SYY) = Precision(1);
            resultant_map((Index)NXY, (Index)SXY) = Precision(1);
            resultant_map((Index)BXX, (Index)SXX) = z;
            resultant_map((Index)BYY, (Index)SYY) = z;
            resultant_map((Index)BXY, (Index)SXY) = z;
            resultant_map((Index)QX,  (Index)SXZ) = shear_correction;
            resultant_map((Index)QY,  (Index)SYZ) = shear_correction;

            tangent_layer.noalias() += w * resultant_map * material_tangent * strain_map;
        }

        // Rotate the layer contribution back into the geometric shell basis.
        resultants_shell.values() += resultants_layer_to_shell * resultants_layer.values();
        tangent_shell.noalias()   += resultants_layer_to_shell * tangent_layer * strain_shell_to_layer;

        z_bottom += layer.thickness;
    }
}

/**
 * Recovers physical Cauchy stress in the layer containing z.
 *
 * The material is evaluated at the Green-Lagrange base strain in the layer axes
 * using the history of the nearest stored layer point. The PK2-to-Cauchy
 * push-forward is differentiated with respect to both Delta S and Delta F, as
 * for the homogeneous integrated section.
 */
VolumeStressCauchy CompositeShellSection::recover_stress(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_base,
    const ShellGeneralizedStrain& strain_increment,
    const Precision*              old_material_state,
    Index                         material_state_stride,
    Precision                     z,
    const Mat3&                   deformation_gradient_base,
    const Mat3&                   deformation_gradient_increment
) const {
    using GeneralizedStrainComponent = ShellGeneralizedStrain::Component;
    using MaterialStrainComponent    = ShellMaterialStrain::Component;
    using MaterialStressComponent    = ShellMaterialStress::Component;

    constexpr GeneralizedStrainComponent EXX = GeneralizedStrainComponent::EpsilonXX;
    constexpr GeneralizedStrainComponent EYY = GeneralizedStrainComponent::EpsilonYY;
    constexpr GeneralizedStrainComponent GXY = GeneralizedStrainComponent::GammaXY;
    constexpr GeneralizedStrainComponent KXX = GeneralizedStrainComponent::KappaXX;
    constexpr GeneralizedStrainComponent KYY = GeneralizedStrainComponent::KappaYY;
    constexpr GeneralizedStrainComponent KXY = GeneralizedStrainComponent::KappaXY;
    constexpr GeneralizedStrainComponent GXZ = GeneralizedStrainComponent::GammaXZ;
    constexpr GeneralizedStrainComponent GYZ = GeneralizedStrainComponent::GammaYZ;

    constexpr MaterialStrainComponent MXX = MaterialStrainComponent::XX;
    constexpr MaterialStrainComponent MYY = MaterialStrainComponent::YY;
    constexpr MaterialStrainComponent MXY = MaterialStrainComponent::GammaXY;
    constexpr MaterialStrainComponent MXZ = MaterialStrainComponent::GammaXZ;
    constexpr MaterialStrainComponent MYZ = MaterialStrainComponent::GammaYZ;

    constexpr MaterialStressComponent SXX = MaterialStressComponent::XX;
    constexpr MaterialStressComponent SYY = MaterialStressComponent::YY;
    constexpr MaterialStressComponent SXY = MaterialStressComponent::XY;
    constexpr MaterialStressComponent SXZ = MaterialStressComponent::XZ;
    constexpr MaterialStressComponent SYZ = MaterialStressComponent::YZ;

    // Locate the layer containing z. Interface points belong to the lower
    // layer; points outside the section are clamped to the outer layers.
    std::size_t layer_index = layers_.size() - 1;
    Index       first_mp    = 0;
    Precision   z_bottom    = Precision(-0.5) * thickness_;

    for (std::size_t i = 0; i + 1 < layers_.size(); ++i) {
        if (z <= z_bottom + layers_[i].thickness) {
            layer_index = i;
            break;
        }
        z_bottom += layers_[i].thickness;
        first_mp += layers_[i].points;
    }

    const Layer& layer = layers_[layer_index];

    // Use the history of the nearest stored point of the layer.
    const Precision dz        = layer.thickness / Precision(layer.points - 1);
    const Precision position  = std::clamp((z - z_bottom) / dz, Precision(0), Precision(layer.points - 1));
    const Index     state_mp  = first_mp + static_cast<Index>(std::lround(position));
    const Precision* old_state = old_material_state + state_mp * material_state_stride;

    const Mat3 output_basis_global  = stress_basis(position_reference, shell_basis_global);
    const Mat3 section_basis_global = orientation_ ? output_basis_global : shell_basis_global;

    const Precision degree = Precision(std::acos(-1.0) / 180.0);
    const Precision c      = std::cos(layer.angle * degree);
    const Precision s      = std::sin(layer.angle * degree);

    Mat3 layer_rotation = Mat3::Identity();
    layer_rotation(0, 0) =  c;
    layer_rotation(0, 1) = -s;
    layer_rotation(1, 0) =  s;
    layer_rotation(1, 1) =  c;

    const Mat3 layer_basis_global = section_basis_global * layer_rotation;

    const Mat2 layer_axes_in_shell =
        shell_basis_global.template block<3, 2>(0, 0).transpose()
        * layer_basis_global.template block<3, 2>(0, 0);

    const ShellGeneralizedStrain base_layer      = strain_base.transformed(layer_axes_in_shell);
    const ShellGeneralizedStrain increment_layer = strain_increment.transformed(layer_axes_in_shell);

    ShellMaterialStrain material_strain_base;
    material_strain_base[MXX] = base_layer[EXX] + z * base_layer[KXX];
    material_strain_base[MYY] = base_layer[EYY] + z * base_layer[KYY];
    material_strain_base[MXY] = base_layer[GXY] + z * base_layer[KXY];
    material_strain_base[MXZ] = base_layer[GXZ];
    material_strain_base[MYZ] = base_layer[GYZ];

    ShellMaterialStrain material_strain_increment;
    material_strain_increment[MXX] = increment_layer[EXX] + z * increment_layer[KXX];
    material_strain_increment[MYY] = increment_layer[EYY] + z * increment_layer[KYY];
    material_strain_increment[MXY] = increment_layer[GXY] + z * increment_layer[KXY];
    material_strain_increment[MXZ] = increment_layer[GXZ];
    material_strain_increment[MYZ] = increment_layer[GYZ];

    ShellMaterialStressPK2 stress_base;
    Mat5                   material_tangent;

    layer.material->elasticity()->evaluate(material_strain_base, old_state, nullptr, stress_base, &material_tangent);

    const ShellMaterialStressPK2 stress_increment(Vec5(material_tangent * material_strain_increment.values()));

    const auto shell_stress_tensor = [=](const ShellMaterialStress& stress) {
        Mat3 tensor = Mat3::Zero();
        tensor(0, 0) = stress[SXX];
        tensor(1, 1) = stress[SYY];
        tensor(0, 1) = stress[SXY];
        tensor(1, 0) = stress[SXY];
        tensor(0, 2) = stress[SXZ];
        tensor(2, 0) = stress[SXZ];
        tensor(1, 2) = stress[SYZ];
        tensor(2, 1) = stress[SYZ];
        return tensor;
    };

    const Mat3 second_pk_base_global =
        layer_basis_global * shell_stress_tensor(stress_base) * layer_basis_global.transpose();
    const Mat3 second_pk_increment_global =
        layer_basis_global * shell_stress_tensor(stress_increment) * layer_basis_global.transpose();

    const Precision J0 = deformation_gradient_base.determinant();
    logging::error(J0 > Precision(0) && std::isfinite(J0),
        "CompositeShellSection: invalid base deformation gradient during stress linearization, J = ", J0);

    const Mat3 sigma_base =
        deformation_gradient_base * second_pk_base_global * deformation_gradient_base.transpose() / J0;

    const Mat3 sigma_increment =
        (deformation_gradient_increment * second_pk_base_global * deformation_gradient_base.transpose()
       + deformation_gradient_base * second_pk_increment_global * deformation_gradient_base.transpose()
       + deformation_gradient_base * second_pk_base_global * deformation_gradient_increment.transpose()) / J0
       - (deformation_gradient_base.inverse() * deformation_gradient_increment).trace() * sigma_base;

    const Mat3 sigma = sigma_base + sigma_increment;
    return VolumeStressCauchy(sigma).transformed(Mat3::Identity(), output_basis_global);
}

material::Material::Ptr CompositeShellSection::material_at(Index mp) const {
    Index first_mp = 0;
    for (const Layer& layer : layers_) {
        if (mp < first_mp + layer.points) {
            return layer.material;
        }
        first_mp += layer.points;
    }

    logging::error(false,
        "CompositeShellSection: material point ", mp, " exceeds the ", num_mp_, " section points");
    return nullptr;
}

Index CompositeShellSection::material_state_size() const {
    Index maximum_size = 0;
    for (const Layer& layer : layers_) {
        maximum_size = std::max(maximum_size, layer.material->elasticity()->state_size());
    }
    return maximum_size;
}

/**
 * Returns the thickness-averaged density of all layers.
 *
 * The constructor currently enforces equal layer densities, so the mean is the
 * common density. The weighted form is kept for a later extension to layer-wise
 * mass integration.
 *
 * @return Mean density, or -1 if any layer has no density.
 */
Precision CompositeShellSection::density() const {
    Precision mass = Precision(0);
    for (const Layer& layer : layers_) {
        if (!layer.material->has_density()) {
            return Precision(-1);
        }
        mass += layer.material->get_density() * layer.thickness;
    }
    return mass / thickness_;
}

void CompositeShellSection::info() {
    ShellSection::info();
    logging::info(true, "   Layers     : ", layers_.size());
    for (const Layer& layer : layers_) {
        logging::info(true, "      ", layer.material->name,
                            ", t=", layer.thickness,
                            ", angle=", layer.angle,
                            ", points=", layer.points,
                            layer.name.empty() ? "" : ", name=", layer.name);
    }
}

} // namespace fem
