
#pragma once

#include "../element/element_structural.h"
#include "../../math/interpolate.h"
#include "../collection/dict.h"
#include "../geometry/surface/surface.h"

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
    SolidElement(ID p_elem_id, std::array<ID, N> p_node_ids)
        : StructuralElement(p_elem_id)
        , node_ids(p_node_ids) {}

    //-------------------------------------------------------------------------
    // Pure Virtual Methods
    //-------------------------------------------------------------------------

    virtual StaticMatrix<N, 1> shape_function(Precision r, Precision s, Precision t) = 0;
    virtual StaticMatrix<N, D> shape_derivative(Precision r, Precision s, Precision t) = 0;
    virtual StaticMatrix<N, D> node_coords_local() = 0;

    template<Dim K>
    StaticVector<K> interpolate(StaticMatrix<N, K> values, Precision r, Precision s, Precision t);

    StaticMatrix<N, D> node_coords_global(NodeData& node_coords) {
        return this->nodal_data<D>(node_coords);
    }

    virtual const quadrature::Quadrature& integration_scheme() = 0;
    virtual const quadrature::Quadrature& integration_scheme_mass() {return this->integration_scheme();}

    //-------------------------------------------------------------------------
    // Strain-Displacement Matrix and Jacobian
    //-------------------------------------------------------------------------

    StaticMatrix<n_strain, D * N> strain_displacement(const StaticMatrix<N, D>& shape_der_global);
    StaticMatrix<n_strain, D * N> strain_displacements(const StaticMatrix<N, D>& node_coords,
          Precision r, Precision s, Precision t, Precision& det, bool check_det = true);
    StaticMatrix<D, D> jacobian(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t);
    StaticMatrix<n_strain, n_strain> mat_matrix(const StaticMatrix<N, D>& node_coords,
          Precision r, Precision s, Precision t);

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

    MapMatrix stiffness(NodeData& position, Precision* buffer) override;
    MapMatrix mass(NodeData& position, Precision* buffer) override;

    Dim dimensions() override;
    Dim n_nodes() override;
    ID* nodes() override;
    virtual SurfacePtr surface(ID surface_id) override = 0;
    Dim n_integration_points() override;

    Precision volume(NodeData& node_coords) override;

    void apply_vload                (NodeData& node_coords, NodeData& node_loads, Vec3 load) override;
    void apply_tload                (NodeData& node_coords, NodeData& node_loads, NodeData& node_temp, Precision ref_temp) override;
    void compute_stress_strain_nodal(NodeData& node_coords, NodeData& displacement, NodeData& stress, NodeData& strain) override;
    void compute_stress_strain      (NodeData& node_coords, NodeData& displacement, NodeData& stress, NodeData& strain, NodeData& xyz) override;
    void compute_compliance         (NodeData& node_coords, NodeData& displacement, ElementData& result) override;

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
