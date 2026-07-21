/**
 * @file frt_shell_reference.inl
 * @brief Implements curved reference-state initialization for finite-rotation shells.
 *
 * The routines construct nodal reference directors and evaluate the actual
 * isoparametric midsurface at all integration and MITC tying points. Every
 * cached point receives its own tangents, orthonormal basis, surface Jacobian,
 * transformed shape derivatives and interpolated reference director field.
 *
 * No planar element projection is used. Curved S6/S8 geometry and warped S4
 * geometry are integrated directly through `|X_,r x X_,s|`.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"
#include "../../math/vec_util.h"

#include <cmath>

namespace fem::model {

using math::normalized;

/**
 * Constructs the common shell base from one element identifier and nodal
 * connectivity.
 *
 * @param id Element identifier.
 * @param nodes Element node identifiers in topology-specific ordering.
 */
template<Index N>
FRTShell<N>::FRTShell(ID id, const std::array<ID, N>& nodes)
    : ShellElement<N>(id, nodes) {}

/**
 * Initializes the persistent curved-reference configuration.
 *
 * The routine gathers undeformed nodal coordinates, constructs one consistent
 * reference director per node, caches the complete geometry at all numerical
 * integration and MITC tying points and integrates the physical midsurface
 * area. Thread-local evaluation memory is independent of this lifecycle and is
 * neither allocated nor released here.
 */
template<Index N>
void FRTShell<N>::step_begin() {
    auto data = std::make_unique<ReferenceData>();

    // Gather the undeformed midsurface geometry and construct the nodal
    // reference directors before any pointwise geometry is evaluated
    data->X  = init_reference_node_coords();
    data->d0 = init_reference_directors(data->X);

    // Publish the partially initialized reference data because
    // make_reference_point() accesses it through reference_data()
    reference_data_ = std::move(data);

    const auto& scheme = this->integration_scheme();
    reference_data_->ip_points.reserve(static_cast<std::size_t>(scheme.count()));

    // Cache the complete curved geometry and accumulate the physical area at
    // every numerical integration point
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        reference_data_->ip_points.push_back(
            make_reference_point(point.r, point.s, point.w)
        );
        reference_data_->area +=
            point.w * reference_data_->ip_points.back().detJ;
    }

    // Cache the topology-specific MITC tying-point geometry with zero
    // quadrature weight because tying points do not contribute area directly
    const std::vector<Vec2> tying_coordinates = tying_point_coordinates();
    reference_data_->tying_points.reserve(tying_coordinates.size());

    for (const Vec2& local : tying_coordinates) {
        reference_data_->tying_points.push_back(
            make_reference_point(local(0), local(1), Precision(0))
        );
    }


}

/**
 * Releases the persistent element reference geometry.
 *
 * The reusable `thread_local` evaluation workspaces remain owned by their
 * worker threads and intentionally survive individual analysis steps.
 */
template<Index N>
void FRTShell<N>::step_end() {
    reference_data_.reset();
}

/**
 * Returns the initialized persistent reference configuration.
 *
 * @return Constant element reference data.
 */
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

    // Gather the global reference positions in local element-node ordering
    for (Index node = 0; node < num_nodes; ++node) {
        X.row(node) = positions.row_vec3(
            static_cast<Index>(this->node_ids[node])
        ).transpose();
    }

    return X;
}

/**
 * Constructs one unit reference director at every element node.
 *
 * Element-nodal normals supplied by the model take precedence. Otherwise the
 * director is obtained from the cross product of the actual isoparametric
 * reference tangents evaluated at the natural node coordinate. This preserves
 * the curved geometry represented by midside nodes of S6 and S8 elements.
 *
 * All nodal directors are finally oriented consistently with the first node.
 *
 * @param X Global reference midsurface coordinates.
 * @return Unit reference director at every shell node.
 */
template<Index N>
typename FRTShell<N>::MatN3 FRTShell<N>::init_reference_directors(
    const MatN3& X
) const {
    MatN3 d0;

    const bool has_element_normals =
        this->_model_data != nullptr &&
        this->_model_data->shell_element_nodal_normals != nullptr;

    const MatN2 natural_nodes = node_coords_natural();

    for (Index node = 0; node < num_nodes; ++node) {
        Vec3 director = Vec3::Zero();

        if (has_element_normals) {
            // Use explicitly supplied element-nodal normals when available
            const Field& normals = *this->_model_data->shell_element_nodal_normals;
            const Index row = static_cast<Index>(this->elem_nodal_offset) + node;

            director = normals.row_vec3(row);
            logging::error(
                director.allFinite() && director.norm() > Precision(0),
                "FRTShell: invalid element-nodal shell normal for element ",
                this->elem_id
            );
        } else {
            // Recover the normal from the actual topology-specific
            // isoparametric tangents at the natural node coordinate
            const Precision r = natural_nodes(node, 0);
            const Precision s = natural_nodes(node, 1);
            const MatN2 dshape = shape_derivative(r, s);

            Vec3 X_r = Vec3::Zero();
            Vec3 X_s = Vec3::Zero();

            for (Index i = 0; i < num_nodes; ++i) {
                const Vec3 X_i = X.row(i).transpose();
                X_r += dshape(i, 0) * X_i;
                X_s += dshape(i, 1) * X_i;
            }

            director = X_r.cross(X_s);
            logging::error(
                director.allFinite() && director.norm() > Precision(1e-14),
                "FRTShell: singular reference geometry at node ",
                node,
                " of element ",
                this->elem_id
            );
        }

        director.normalize();
        d0.row(node) = director.transpose();
    }

    // Use the first nodal director to define the positive shell side and avoid
    // local sign reversals inside one higher-order element
    const Vec3 orientation = d0.row(0).transpose();

    for (Index node = 1; node < num_nodes; ++node) {
        if (d0.row(node).transpose().dot(orientation) < Precision(0)) {
            d0.row(node) *= Precision(-1);
        }
    }

    return d0;
}

/**
 * Evaluates the complete curved-reference geometry at one natural point.
 *
 * The physical surface measure is `|X_,r x X_,s|`. The first local basis vector
 * follows the first reference tangent, the third vector is the surface normal
 * and the second vector completes a right-handed orthonormal basis. Natural
 * shape derivatives are transformed into this pointwise tangent system.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param w Quadrature weight; zero for tying and arbitrary output points.
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

    logging::error(
        point.detJ > Precision(1e-14) && std::isfinite(point.detJ),
        "FRTShell: singular reference surface Jacobian in element ",
        this->elem_id,
        " at (",
        r,
        ", ",
        s,
        ")"
    );

    // Construct the pointwise orthonormal basis on the curved reference surface
    point.e1 = normalized(point.X_r);
    point.e3 = normal / point.detJ;
    point.e2 = normalized(point.e3.cross(point.e1));
    point.e1 = normalized(point.e2.cross(point.e3));

    // Express the two natural tangents in the local orthonormal tangent basis
    point.A << point.X_r.dot(point.e1), point.X_r.dot(point.e2),
               point.X_s.dot(point.e1), point.X_s.dot(point.e2);

    const Precision detA = point.A.determinant();
    logging::error(
        std::abs(detA) > Precision(1e-14),
        "FRTShell: singular local reference mapping in element ",
        this->elem_id,
        " at (",
        r,
        ", ",
        s,
        ")"
    );

    point.invA = point.A.inverse();

    // Transform every shape gradient into the orthonormal tangent coordinates
    for (Index node = 0; node < num_nodes; ++node) {
        const Vec2 derivative_nat = point.dshape_rs.row(node).transpose();
        const Vec2 derivative_loc = point.invA * derivative_nat;

        point.dshape_a(node) = derivative_loc(0);
        point.dshape_b(node) = derivative_loc(1);
    }

    // The reference derivatives with respect to the orthonormal coordinates are
    // exactly the local basis vectors by construction
    point.X_a = point.e1;
    point.X_b = point.e2;

    // Interpolate the nodal reference director field and all natural and local
    // tangent derivatives required by the shell curvature and shear measures
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

/**
 * Interpolates the physical reference midsurface position.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @return Global reference position.
 */
template<Index N>
Vec3 FRTShell<N>::reference_position(Precision r, Precision s) const {
    const VecN shape = shape_function(r, s);
    const MatN3& X   = reference_data().X;

    Vec3 position = Vec3::Zero();

    // Interpolate directly from the persistent nodal reference coordinates
    for (Index node = 0; node < num_nodes; ++node) {
        position += shape(node) * X.row(node).transpose();
    }

    return position;
}

/**
 * Returns the pointwise global orthonormal reference basis.
 *
 * Cached integration and tying points are reused. Arbitrary output coordinates
 * evaluate one temporary reference point without modifying persistent storage.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @return Matrix whose columns are `[e1,e2,e3]`.
 */
template<Index N>
Mat3 FRTShell<N>::reference_basis_global(Precision r, Precision s) const {
    const ReferencePoint* cached = cached_reference_point(r, s);
    const ReferencePoint temporary =
        cached ? ReferencePoint{} : make_reference_point(r, s, Precision(0));
    const ReferencePoint& point = cached ? *cached : temporary;

    Mat3 basis;
    basis.col(0) = point.e1;
    basis.col(1) = point.e2;
    basis.col(2) = point.e3;
    return basis;
}

/**
 * Searches the cached integration and tying points for one natural coordinate
 * pair.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @return Pointer to the cached point or `nullptr` when no match exists.
 */
template<Index N>
const typename FRTShell<N>::ReferencePoint* FRTShell<N>::cached_reference_point(
    Precision r,
    Precision s
) const {
    const Precision tolerance = Precision(1e-14);
    const auto& ref = reference_data();

    // Search numerical integration points first because most material and
    // output requests refer to these coordinates
    for (const ReferencePoint& point : ref.ip_points) {
        if (std::abs(point.r - r) <= tolerance &&
            std::abs(point.s - s) <= tolerance) {
            return &point;
        }
    }

    // Search the topology-specific MITC tying points second
    for (const ReferencePoint& point : ref.tying_points) {
        if (std::abs(point.r - r) <= tolerance &&
            std::abs(point.s - s) <= tolerance) {
            return &point;
        }
    }

    return nullptr;
}

/**
 * Returns the shell section assigned to this element.
 *
 * @return Valid shell-section pointer.
 */
template<Index N>
ShellSection* FRTShell<N>::shell_section() const {
    logging::error(this->_section != nullptr,
                   "FRTShell: section not set for element ", this->elem_id);
    logging::error(this->_section->template as<ShellSection>() != nullptr,
                   "FRTShell: section is not a shell section for element ",
                   this->elem_id);

    return this->_section->template as<ShellSection>();
}

/**
 * Returns the current optional element stiffness scaling factor.
 *
 * The value defaults to one and may be supplied by a scalar element field, for
 * example during topology optimization.
 *
 * @return Current dimensionless element stiffness scale.
 */
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

/**
 * Evaluates the zero-strain generalized shell-section tangent.
 *
 * The physical reference position and pointwise material basis are supplied to
 * the section, and the current topology stiffness scale is applied to the
 * returned tangent.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @return Scaled generalized section tangent.
 */
template<Index N>
typename FRTShell<N>::Mat8 FRTShell<N>::resultant_stiffness(
    Precision r,
    Precision s
) const {
    ShellGeneralizedStrain zero_strain(Vec8::Zero());
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
