/******************************************************************************
 * @file SolidElement.hpp
 * @brief Defines the SolidElement template class, which serves as a base 
 * class for solid finite elements in a 3D space, with methods to calculate 
 * shape functions, strain-displacement matrices, Jacobians, and stiffness/mass 
 * matrices.
 *
 * @tparam N The number of nodes associated with the solid element.
 *
 * @details This class provides a generic interface for solid elements in FEM. 
 * It defines functions for computing shape functions and derivatives, global 
 * and local node coordinates, Jacobians, strain-displacement matrices, 
 * stiffness/mass matrices, and various element-related computations.
 *
 * @date Created on 12.06.2023
 * @author Finn Eggers
 ******************************************************************************/

#pragma once

#include "element.h"
#include "../math/interpolate.h"

namespace fem::model {

/******************************************************************************
 * @class SolidElement
 * @brief Base class template for a solid element in finite element analysis.
 * 
 * @tparam N The number of nodes in the element.
 * 
 * @details This class defines general methods that a solid element must
 * implement, including methods for calculating shape functions, Jacobians,
 * stiffness matrices, and mass matrices. Derived classes should provide
 * concrete implementations of the shape functions and integration schemes.
 ******************************************************************************/
template<Index N>
struct SolidElement : public ElementInterface {
    
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
        : ElementInterface(p_elem_id)
        , node_ids(p_node_ids) {}

    //-------------------------------------------------------------------------
    // Pure Virtual Methods
    //-------------------------------------------------------------------------
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
    virtual StaticMatrix<N, D> node_coords_global(NodeData& node_coords) {
        return this->nodal_data<D>(node_coords);
    }

    /**
     * @brief Returns the quadrature integration scheme to be used for integration over the element.
     *
     * @return const quadrature::Quadrature& The quadrature rule.
     */
    virtual const quadrature::Quadrature& integration_scheme() = 0;

    //-------------------------------------------------------------------------
    // Strain-Displacement Matrix and Jacobian
    //-------------------------------------------------------------------------
    /**
     * @brief Computes the strain-displacement matrix (B-matrix) based on the shape function derivatives in global coordinates.
     *
     * @param shape_der_global The shape function derivatives in global coordinates.
     * @return StaticMatrix<n_strain, D * N> The computed strain-displacement matrix.
     */
    virtual StaticMatrix<n_strain, D * N> strain_displacement(const StaticMatrix<N, D>& shape_der_global);

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
    virtual StaticMatrix<n_strain, D * N> strain_displacements(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t, Precision& det, bool check_det = true);

    /**
     * @brief Computes the Jacobian matrix at the given local coordinates.
     *
     * @param node_coords The coordinates of the element nodes.
     * @param r Local coordinate along the r-direction.
     * @param s Local coordinate along the s-direction.
     * @param t Local coordinate along the t-direction.
     * @return StaticMatrix<D, D> The Jacobian matrix.
     */
    virtual StaticMatrix<D, D> jacobian(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t);

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
    StaticMatrix<N, K> nodal_data(const NodeData& full_data, Index offset = 0, Index stride = 1);

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
    MapMatrix stiffness(NodeData& position, Precision* buffer) override;

    /**
     * @brief Computes the mass matrix for the element.
     *
     * @param position The nodal positions for the element.
     * @param buffer A buffer to store the mass matrix.
     * @return MapMatrix The computed mass matrix.
     */
    MapMatrix mass(NodeData& position, Precision* buffer) override;

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
    Precision volume(NodeData& node_coords) override;

    /**
     * @brief Applies a volume load to the element.
     *
     * @param node_coords The global nodal coordinates for the element.
     * @param node_loads The nodal loads that will be updated by the applied load.
     * @param load The external load vector.
     */
    void apply_vload(NodeData& node_coords, NodeData& node_loads, StaticVector<3> load) override;

    /**
     * @brief Applies a distributed load to a specific surface of the element (not yet implemented).
     *
     * @param node_coords The global nodal coordinates for the element.
     * @param node_loads The nodal loads that will be updated by the applied load.
     * @param surface The surface index on which the distributed load is applied.
     * @param load The external load vector.
     */
    void apply_dload(NodeData& node_coords, NodeData& node_loads, ID surface, StaticVector<3> load) override;

    /**
     * @brief Computes the nodal stress and strain for the element.
     *
     * @param node_coords The global nodal coordinates for the element.
     * @param displacement The nodal displacement data.
     * @param stress The computed nodal stress.
     * @param strain The computed nodal strain.
     */
    void compute_stress_strain_nodal(NodeData& node_coords, NodeData& displacement, NodeData& stress, NodeData& strain) override;

    /**
     * @brief Computes the stress and strain at integration points of the element.
     *
     * @param node_coords The global nodal coordinates for the element.
     * @param displacement The nodal displacement data.
     * @param stress The computed stress at integration points.
     * @param strain The computed strain at integration points.
     * @param xyz The computed global coordinates of the integration points.
     */
    void compute_stress_strain(NodeData& node_coords, NodeData& displacement, NodeData& stress, NodeData& strain, NodeData& xyz) override;

    /**
     * @brief Computes the compliance (strain energy) for the element.
     *
     * @param node_coords The global nodal coordinates for the element.
     * @param displacement The nodal displacement data.
     * @param result The computed compliance value.
     */
    void compute_compliance(NodeData& node_coords, NodeData& displacement, ElementData& result) override;

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


#include "element_solid.ipp"
