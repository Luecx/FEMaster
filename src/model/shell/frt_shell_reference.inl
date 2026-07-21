/**
 * @file frt_shell_reference.cpp
 * @brief Implements reference-state initialization for finite-rotation shells.
 *
 * The routines construct nodal reference directors and evaluate the actual
 * isoparametric reference midsurface at all integration and MITC tying points.
 * Every point receives its own tangents, orthonormal basis, surface Jacobian
 * and transformed shape-function derivatives.
 *
 * No planar element projection is used. Curved S6 and S8 geometry and warped
 * S4 geometry are integrated directly through `|X_,r x X_,s|`.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"
#include "../../math/vec_util.h"

#include <cmath>

namespace fem::model {

using math::normalized;

template<Index N>
FRTShell<N>::VectorDerivatives::VectorDerivatives() {
    for (auto& value : d2) {
        value.setZero();
    }
}

template<Index N>
FRTShell<N>::EvaluationData::EvaluationData(Index num_ip, Index num_tying)
    : tying_strain_nat(static_cast<std::size_t>(num_tying), Vec8::Zero()),
      tying_B_nat     (static_cast<std::size_t>(num_tying), Mat8x6N::Zero()),
      tying_G_nat     (static_cast<std::size_t>(num_tying)),
      ip_strain       (static_cast<std::size_t>(num_ip), Vec8::Zero()),
      ip_B            (static_cast<std::size_t>(num_ip), Mat8x6N::Zero()),
      ip_G            (static_cast<std::size_t>(num_ip)),
      ip_resultants   (static_cast<std::size_t>(num_ip), Vec8::Zero()),
      ip_tangent      (static_cast<std::size_t>(num_ip), Mat8::Zero()),
      ip_weight       (static_cast<std::size_t>(num_ip), Precision(0)) {
    // Initialize every strain Hessian explicitly because std::array does not
    // guarantee zero initialization for the fixed-size Eigen matrices
    for (auto& point_G : tying_G_nat) {
        for (auto& value : point_G) {
            value.setZero();
        }
    }

    for (auto& point_G : ip_G) {
        for (auto& value : point_G) {
            value.setZero();
        }
    }
}

template<Index N>
FRTShell<N>::FRTShell(ID id, const std::array<ID, N>& nodes)
    : ShellElement<N>(id, nodes) {}

/**
 * Initializes the reference configuration used by the complete analysis step.
 *
 * The routine collects the reference nodal coordinates, constructs one
 * reference director per element node and caches the pointwise curved-surface
 * geometry at all quadrature and MITC tying points. The reference area is
 * obtained by integrating the actual isoparametric surface Jacobian.
 */
template<Index N>
void FRTShell<N>::step_begin() {
    auto data = std::make_unique<ReferenceData>();

    // Collect the undeformed midsurface geometry and nodal reference directors
    data->X  = init_reference_node_coords();
    data->d0 = init_reference_directors(data->X);

    // Cache the complete reference geometry at every quadrature point
    const auto& scheme = this->integration_scheme();
    data->ip_points.reserve(static_cast<std::size_t>(scheme.count()));

    reference_data_ = std::move(data);

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        reference_data_->ip_points.push_back(make_reference_point(point.r, point.s, point.w));
        reference_data_->area += point.w * reference_data_->ip_points.back().detJ;
    }

    // Cache the element-specific MITC tying points with zero quadrature weight
    const std::vector<Vec2> tying_coordinates = tying_point_coordinates();
    reference_data_->tying_points.reserve(tying_coordinates.size());

    for (const Vec2& local : tying_coordinates) {
        reference_data_->tying_points.push_back(
            make_reference_point(local(0), local(1), Precision(0))
        );
    }
}

template<Index N>
void FRTShell<N>::step_end() {
    reference_data_.reset();
}

template<Index N>
const typename FRTShell<N>::ReferenceData& FRTShell<N>::reference_data() const {
    logging::error(
        reference_data_ != nullptr,
        "FRTShell: reference data has not been initialized for element ",
        this->elem_id,
        ". Call step_begin() before evaluating the element."
    );

    return *reference_data_;
}

/**
 * Collects the undeformed global midsurface coordinates of the element nodes.
 *
 * @return Matrix containing one global reference position per shell node.
 */
template<Index N>
typename FRTShell<N>::MatN3 FRTShell<N>::init_reference_node_coords() const {
    logging::error(this->_model_data                      != nullptr,
                   "FRTShell: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions_reference != nullptr,
                   "FRTShell: POSITION_REFERENCE field is not set");

    const auto& positions = *this->_model_data->positions_reference;
    MatN3       X;

    for (Index node = 0; node < num_nodes; ++node) {
        X.row(node) = positions.row_vec3(static_cast<Index>(this->node_ids[node])).transpose();
    }

    return X;
}

/**
 * Constructs the nodal reference directors.
 *
 * Element-nodal normals supplied by the model take precedence. Otherwise, the
 * director at each node is obtained from the isoparametric reference tangents
 * evaluated at the natural node location. This preserves the curved geometry
 * represented by the midside nodes of S6 and S8 elements.
 *
 * @param X Global reference midsurface coordinates.
 *
 * @return Unit reference director at every shell node.
 */
template<Index N>
typename FRTShell<N>::MatN3 FRTShell<N>::init_reference_directors(const MatN3& X) const {
    MatN3 d0;

    const bool has_element_normals =
        this->_model_data != nullptr &&
        this->_model_data->shell_element_nodal_normals != nullptr;

    const MatN2 natural_nodes = node_coords_natural();

    for (Index node = 0; node < num_nodes; ++node) {
        Vec3 director = Vec3::Zero();

        if (has_element_normals) {
            const Field& normals = *this->_model_data->shell_element_nodal_normals;
            const Index  row     = static_cast<Index>(this->elem_nodal_offset) + node;

            director = normals.row_vec3(row);
            logging::error(director.allFinite() && director.norm() > Precision(0),
                           "FRTShell: invalid element-nodal shell normal for element ",
                           this->elem_id);
        } else {
            const Precision r = natural_nodes(node, 0);
            const Precision s = natural_nodes(node, 1);
            const MatN2     dshape = shape_derivative(r, s);

            Vec3 X_r = Vec3::Zero();
            Vec3 X_s = Vec3::Zero();

            for (Index i = 0; i < num_nodes; ++i) {
                const Vec3 X_i = X.row(i).transpose();
                X_r += dshape(i, 0) * X_i;
                X_s += dshape(i, 1) * X_i;
            }

            director = X_r.cross(X_s);
            logging::error(director.allFinite() && director.norm() > Precision(1e-14),
                           "FRTShell: singular reference geometry at node ", node,
                           " of element ", this->elem_id);
        }

        director.normalize();
        d0.row(node) = director.transpose();
    }

    // Keep the director orientation consistent throughout the element. The
    // first nodal director defines the positive side of the shell.
    const Vec3 orientation = d0.row(0).transpose();

    for (Index node = 1; node < num_nodes; ++node) {
        if (d0.row(node).transpose().dot(orientation) < Precision(0)) {
            d0.row(node) *= Precision(-1);
        }
    }

    return d0;
}

/**
 * Evaluates the complete reference geometry at one natural point.
 *
 * The physical surface measure is `|X_,r x X_,s|`. The first local basis
 * vector follows the first reference tangent, the third basis vector is the
 * surface normal and the second vector completes a right-handed orthonormal
 * basis. Shape derivatives are transformed into this local tangent system.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param w Quadrature weight; zero for tying and output points.
 *
 * @return Pointwise reference geometry and director interpolation.
 */
template<Index N>
typename FRTShell<N>::ReferencePoint FRTShell<N>::make_reference_point(
    Precision r,
    Precision s,
    Precision w
) const {
    const auto& ref = reference_data();

    ReferencePoint point;
    point.r         = r;
    point.s         = s;
    point.w         = w;
    point.shape     = shape_function(r, s);
    point.dshape_rs = shape_derivative(r, s);

    // Interpolate the actual isoparametric reference tangents
    for (Index node = 0; node < num_nodes; ++node) {
        const Vec3 X_i = ref.X.row(node).transpose();
        point.X_r += point.dshape_rs(node, 0) * X_i;
        point.X_s += point.dshape_rs(node, 1) * X_i;
    }

    const Vec3 normal = point.X_r.cross(point.X_s);
    point.detJ = normal.norm();

    logging::error(point.detJ > Precision(1e-14) && std::isfinite(point.detJ),
                   "FRTShell: singular reference surface Jacobian in element ",
                   this->elem_id,
                   " at (", r, ", ", s, ")");

    // Construct a pointwise orthonormal reference basis on the curved surface
    point.e1 = normalized(point.X_r);
    point.e3 = normal / point.detJ;
    point.e2 = normalized(point.e3.cross(point.e1));
    point.e1 = normalized(point.e2.cross(point.e3));

    // Map local tangent derivatives to natural derivatives
    point.A << point.X_r.dot(point.e1), point.X_r.dot(point.e2),
               point.X_s.dot(point.e1), point.X_s.dot(point.e2);

    const Precision detA = point.A.determinant();
    logging::error(std::abs(detA) > Precision(1e-14),
                   "FRTShell: singular local reference mapping in element ",
                   this->elem_id,
                   " at (", r, ", ", s, ")");

    point.invA = point.A.inverse();

    // Transform all shape gradients into the local tangent basis
    for (Index node = 0; node < num_nodes; ++node) {
        const Vec2 derivative_nat = point.dshape_rs.row(node).transpose();
        const Vec2 derivative_loc = point.invA * derivative_nat;

        point.dshape_a(node) = derivative_loc(0);
        point.dshape_b(node) = derivative_loc(1);
    }

    point.X_a = point.e1;
    point.X_b = point.e2;

    // Interpolate the reference director and all required derivatives
    for (Index node = 0; node < num_nodes; ++node) {
        const Vec3 D_i = ref.d0.row(node).transpose();

        point.D   += point.shape(node)        * D_i;
        point.D_r += point.dshape_rs(node, 0) * D_i;
        point.D_s += point.dshape_rs(node, 1) * D_i;
        point.D_a += point.dshape_a(node)     * D_i;
        point.D_b += point.dshape_b(node)     * D_i;
    }

    return point;
}

template<Index N>
Vec3 FRTShell<N>::reference_position(Precision r, Precision s) const {
    const VecN  shape = shape_function(r, s);
    const MatN3 X     = reference_data().X;

    Vec3 position = Vec3::Zero();

    for (Index node = 0; node < num_nodes; ++node) {
        position += shape(node) * X.row(node).transpose();
    }

    return position;
}

template<Index N>
Mat3 FRTShell<N>::reference_basis_global(Precision r, Precision s) const {
    const ReferencePoint* cached = cached_reference_point(r, s);
    const ReferencePoint  temporary = cached ? ReferencePoint{} : make_reference_point(r, s, Precision(0));
    const ReferencePoint& point = cached ? *cached : temporary;

    Mat3 basis;
    basis.col(0) = point.e1;
    basis.col(1) = point.e2;
    basis.col(2) = point.e3;
    return basis;
}

template<Index N>
const typename FRTShell<N>::ReferencePoint* FRTShell<N>::cached_reference_point(
    Precision r,
    Precision s
) const {
    const Precision tol = Precision(1e-14);
    const auto&     ref = reference_data();

    for (const ReferencePoint& point : ref.ip_points) {
        if (std::abs(point.r - r) <= tol && std::abs(point.s - s) <= tol) {
            return &point;
        }
    }

    for (const ReferencePoint& point : ref.tying_points) {
        if (std::abs(point.r - r) <= tol && std::abs(point.s - s) <= tol) {
            return &point;
        }
    }

    return nullptr;
}

template<Index N>
ShellSection* FRTShell<N>::shell_section() const {
    if (!this->_section) {
        logging::error(false, "FRTShell: section not set for element ", this->elem_id);
    }

    if (!this->_section->template as<ShellSection>()) {
        logging::error(false, "FRTShell: section is not a shell section for element ", this->elem_id);
    }

    return this->_section->template as<ShellSection>();
}

template<Index N>
Precision FRTShell<N>::topology_stiffness_scale() const {
    Precision scale = Precision(1);

    if (this->_model_data && this->_model_data->element_stiffness_scale) {
        const auto field = this->_model_data->element_stiffness_scale;
        logging::error(field->components == 1,
                       "Field '", field->name,
                       "': element stiffness scale requires one component");
        scale = (*field)(static_cast<Index>(this->elem_id));
    }

    return scale;
}

template<Index N>
typename FRTShell<N>::Mat8 FRTShell<N>::resultant_stiffness(Precision r, Precision s) const {
    ShellGeneralizedStrain zero_strain;
    ShellStressResultants  zero_resultants;
    Mat8                   H;

    shell_section()->evaluate(
        reference_position(r, s),
        reference_basis_global(r, s),
        zero_strain,
        false,
        zero_resultants,
        H
    );

    return topology_stiffness_scale() * H;
}

} // namespace fem::model
