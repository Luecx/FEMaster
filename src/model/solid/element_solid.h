
#pragma once

#include "../element/element_structural.h"
#include "../../math/interpolate.h"
#include "../geometry/surface/surface.h"
#include <functional>

namespace fem::model {

/**
 * @class SolidElement
 * @brief Base class template for a solid element in finite element analysis.
 *
 * @tparam N The number of nodes in the element.
 *
 * @details This class defines general methods that a solid element must
 * implement, including methods for calculating shape functions, Jacobians,
 * stiffness matrices, and mass matrices. Derived classes should provide
 * concrete implementations of the shape functions and integration schemes.
 */
template<Index N>
struct SolidElement : public StructuralElement {

    //-------------------------------------------------------------------------
    // Member Variables
    //-------------------------------------------------------------------------
    /**
     * @brief Stores the IDs of the nodes in the element.
     */
    std::array<ID, N> node_ids{};

protected:
    static constexpr Dim n_strain = 6; ///< Number of strain components (Voigt notation).
    static constexpr Dim D = 3;        ///< Dimensionality of the element (3D).

    friend quadrature::Quadrature;

public:
    //-------------------------------------------------------------------------
    // Constructor
    //-------------------------------------------------------------------------
    /**
     * @brief Constructs a SolidElement with the given element ID and node IDs.
     *
     * @param p_elem_id The unique ID of the element.
     * @param p_node_ids An array containing the IDs of the nodes that define the element.
     */
    SolidElement(ID p_elem_id, std::array<ID, N> p_node_ids)
        : StructuralElement(p_elem_id)
        , node_ids(p_node_ids) {}

    /**
     * @brief Function to compute the material matrix for the element.
     * Includes the material properties and the material matrix.
     * In case of topology optimization, it checks for fields which would account for rotation
     * and scaling of the material matrix.
     */
    StaticMatrix<n_strain, n_strain> material_matrix(Precision r, Precision s, Precision t);

    /**
     * @brief Computes the shape functions at the given local coordinates.
     *
     * @param r Local coordinate along the r-direction.
     * @param s Local coordinate along the s-direction.
     * @param t Local coordinate along the t-direction.
     * @return StaticMatrix<N, 1> The computed shape function values.
     */
    virtual StaticMatrix<N, 1> shape_function(Precision r, Precision s, Precision t) = 0;

    /**
     * @brief Computes the derivatives of the shape functions with respect to the local coordinates.
     *
     * @param r Local coordinate along the r-direction.
     * @param s Local coordinate along the s-direction.
     * @param t Local coordinate along the t-direction.
     * @return StaticMatrix<N, D> The computed shape function derivatives.
     */
    virtual StaticMatrix<N, D> shape_derivative(Precision r, Precision s, Precision t) = 0;

    /**
     * @brief Interpolates the given data at the given local coordinates.
     *
     * @tparam K The dimension of the data to be interpolated.
     * @param data The data to be interpolated.
     * @param r Local coordinate along the r-direction.
     * @param s Local coordinate along the s-direction.
     * @param t Local coordinate along the t-direction.
     * @return StaticVector<K> The interpolated data.
     */
    template<Dim K>
    StaticVector<K> interpolate(StaticMatrix<N, K> data, Precision r, Precision s, Precision t);

    /**
     * @brief Returns the local coordinates of the nodes in the element.
     *
     * @return StaticMatrix<N, D> The local coordinates of the nodes.
     */
    virtual StaticMatrix<N, D> node_coords_local() = 0;

    /**
     * @brief Returns the global coordinates of the nodes based on provided nodal data.
     *
     * @param node_coords The global nodal data.
     * @return StaticMatrix<N, D> The global coordinates of the nodes.
     */
    virtual StaticMatrix<N, D> node_coords_global() {
        logging::error(this->_model_data != nullptr, "no model data assigned to element ", this->elem_id);
        logging::error(this->_model_data->positions != nullptr, "positions field not set in model data");
        StaticMatrix<N, D> res {};
        const auto& positions = *this->_model_data->positions;
        for (Index i = 0; i < N; ++i) {
            const Index row = static_cast<Index>(this->node_ids[i]);
            res.row(i) = positions.row_vec3(row).transpose();
        }
        return res;
    }

    /**
     * @brief Returns the quadrature integration scheme to be used for integration over the element.
     *
     * @return const quadrature::Quadrature& The quadrature rule.
     */
    virtual const quadrature::Quadrature& integration_scheme() = 0;
    virtual const quadrature::Quadrature& integration_scheme_mass() {return this->integration_scheme();}

    //-------------------------------------------------------------------------
    // Strain-Displacement Matrix and Jacobian
    //-------------------------------------------------------------------------
    /**
     * @brief Computes the strain-displacement matrix (B-matrix) based on the shape function derivatives in global coordinates.
     *
     * @param shape_der_global The shape function derivatives in global coordinates.
     * @return StaticMatrix<n_strain, D * N> The computed strain-displacement matrix.
     */
    StaticMatrix<n_strain, D * N> strain_displacement(const StaticMatrix<N, D>& shape_der_global);

    /**
     * @brief Computes the strain-displacement matrix (B-matrix) and determinant of the Jacobian.
     *
     * @param node_coords The coordinates of the element nodes.
     * @param r Local coordinate along the r-direction.
     * @param s Local coordinate along the s-direction.
     * @param t Local coordinate along the t-direction.
     * @param det Output variable for the determinant of the Jacobian.
     * @param check_det Flag to indicate whether to check for a positive determinant.
     * @return StaticMatrix<n_strain, D * N> The computed strain-displacement matrix.
     */
    StaticMatrix<n_strain, D * N> strain_displacements(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t, Precision& det, bool check_det = true);

    /**
     * @brief Computes the Jacobian matrix at the given local coordinates.
     *
     * @param node_coords The coordinates of the element nodes.
     * @param r Local coordinate along the r-direction.
     * @param s Local coordinate along the s-direction.
     * @param t Local coordinate along the t-direction.
     * @return StaticMatrix<D, D> The Jacobian matrix.
     */
    StaticMatrix<D, D> jacobian(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t);

    //-------------------------------------------------------------------------
    // Helper Functions
    //-------------------------------------------------------------------------
    /**
     * @brief Extracts nodal data for the element from the full global nodal data.
     *
     * @tparam K The dimension of the nodal data to be extracted.
     * @param full_data The full global nodal data.
     * @param offset The offset in the data.
     * @param stride The stride for accessing the data.
     * @return StaticMatrix<N, K> The extracted nodal data for the element.
     */
    template<Dim K>
    StaticMatrix<N, K> nodal_data(const Field& full_data, Index offset = 0, Index stride = 1);

    //-------------------------------------------------------------------------
    // Element Interface Overrides
    //-------------------------------------------------------------------------
    /**
     * @brief Returns the degrees of freedom for the element.
     *
     * @return ElDofs The degrees of freedom for the element (dx, dy, dz enabled).
     */
    ElDofs dofs() override;

    /**
     * @brief Computes the stiffness matrix for the element.
     *
     * @param position The nodal positions for the element.
     * @param buffer A buffer to store the stiffness matrix.
     * @return MapMatrix The computed stiffness matrix.
     */
    MapMatrix stiffness(Precision* buffer) override;

    /**
     * @brief computes the geometric stiffness matrix for the element.
     * @param buffer
     * @param ip_stress
     * @param ip_start_idx
     * @return
     */
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override;

    /**
     * @brief Computes the mass matrix for the element.
     *
     * @param position The nodal positions for the element.
     * @param buffer A buffer to store the mass matrix.
     * @return MapMatrix The computed mass matrix.
     */
    MapMatrix mass(Precision* buffer) override;

    /**
     * @brief Returns the dimensionality of the element (typically 3D).
     *
     * @return Dim The dimensionality of the element.
     */
    Dim dimensions() override;

    /**
     * @brief Returns the number of nodes in the element.
     *
     * @return Dim The number of nodes in the element.
     */
    Dim n_nodes() override;

    /**
     * @brief Returns the array of node IDs for the element.
     *
     * @return ID* A pointer to the array of node IDs.
     */
    ID* nodes() override;

    /**
     * @brief Returns a shared pointer to a surface.
     *
     * @param surface_id The ID of the surface.
     * @return SurfacePtr A shared pointer to the surface.
     */
    virtual SurfacePtr surface(ID surface_id) override = 0;

    /**
     * @brief Returns the number of integration points in the quadrature scheme.
     *
     * @return Dim The number of integration points.
     */
    Dim n_integration_points() override;

    /**
     * @brief Computes the volume of the element.
     *
     * @param node_coords The global nodal coordinates for the element.
     * @return Precision The computed volume.
     */
    Precision volume() override;

    /**
     * @brief Applies a volume load to the element.
     *
     * @param node_coords The global nodal coordinates for the element.
     * @param node_loads The nodal loads that will be updated by the applied load.
     * @param load The external load vector.
     */
    void apply_vload(Field& node_loads, Vec3 load) override;

    /**
     * @brief Integrates a vector field f(x) over the element volume and scatters equivalent nodal loads.
     *
     * @param node_loads Reference to the global node load field (NODE x 6).
     * @param scale_by_density If true, multiply f(x) by the element material density.
     * @param field A callable f(x) with x in global coordinates returning a Vec3 body force density.
     */
    void integrate_vec_field(Field& node_loads,
                             bool scale_by_density,
                             const VecField& field) override;

    /**
     * @brief Applies a thermal load to the element.
     */
    void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) override;

    /**
     * @brief Computes the nodal stress and strain for the element.
     *
     * @param node_coords The global nodal coordinates for the element.
     * @param displacement The nodal displacement data.
     * @param stress The computed nodal stress.
     * @param strain The computed nodal strain.
     */
    void compute_stress_strain_nodal(Field& displacement, Field& stress, Field& strain) override;

    /**
     * Computes the stress at a given point in the element.
     * @param displacement
     * @param xyz
     * @return
     */
    Stresses stress(Field& displacement, std::vector<Vec3>& rst) override;
    Strains  strain(Field& displacement, std::vector<Vec3>& rst) override;

    /**
     * Computes the stress and strain at the integration points of the elements
     * @param ip_stress
     * @param ip_strain
     * @param displacement
     * @param ip_offset
     */
    void compute_stress_strain(Field& ip_stress, Field& ip_strain, Field& displacement, int ip_offset) override;

    /**
     * @brief Computes the compliance (strain energy) for the element.
     *
     * @param node_coords The global nodal coordinates for the element.
     * @param displacement The nodal displacement data.
     * @param result The computed compliance value.
     */
    void compute_compliance(Field& displacement, Field& result) override;

    /**
     * @brief Computes the derivative of the compliance w.r.t the three angles defined as material_orientation.
     * If its not defined, the derivative is zero.
     * @param displacement
     * @param result
     */
    void compute_compliance_angle_derivative(Field& displacement, Field& result) override;

    //-------------------------------------------------------------------------
    // Testing Functions
    //-------------------------------------------------------------------------
    /**
     * @brief Performs various tests to verify the implementation of a given element type.
     *
     * @tparam ElementType The specific element type to be tested.
     */
    template<class ElementType>
    static bool test_implementation(bool print=false);
};


}  // namespace fem::model


#include "element_solid_compute.ipp"
#include "element_solid_load.ipp"
#include "element_solid_.ipp"
#include "element_solid.ipp"
#include "element_solid_test.ipp"
