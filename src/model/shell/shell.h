//
// Created by f_eggers on 11.12.2024.
//

#ifndef SHELL_H
#define SHELL_H

#include "../../material/stress.h"
#include "../../material/elasticity.h"
#include "../../core/core.h"
#include "../element/element_structural.h"
#include "../../section/section_shell.h"

#include <memory>

namespace fem::model {
template<Index N>
struct ShellElement : StructuralElement {
    std::array<ID, N> node_ids;

    ShellElement(ID p_elem_id, std::array<ID, N> p_node_ids) : StructuralElement(p_elem_id), node_ids {p_node_ids} {}
    ~ShellElement() override = default;

    ShellSection* get_section() {
        if (!this->_section) {
            logging::error(false, "Section not set for element ", this->elem_id);
        }
        if (!this->_section->template as<ShellSection>()) {
            logging::error(false, "Section is not a beam section for element ", this->elem_id);
        }
        return this->_section->template as<ShellSection>();
    }
    material::MaterialPtr get_material() {
        ShellSection* section = get_section();
        if (!section->material_) {
            logging::error(false, "Material not set for element ", this->elem_id);
        }
        return section->material_;
    }
    material::Elasticity* get_elasticity() {
        auto mat_ptr = get_material();
        if (!mat_ptr->has_elasticity()) {
            logging::error(false, "Material has no elasticity assigned");
        }
        return mat_ptr->elasticity().get();
    }
    Precision get_density(bool required = false) {
        ShellSection* section = get_section();
        if (required) {
            logging::error(section->has_density(),
                           "ShellElement: density is required for element ", this->elem_id);
        }
        return section->get_density();
    }

    // left out for childs
    SurfacePtr surface(ID surface_id) override = 0;
    Precision  volume() override = 0;
    MapMatrix  stiffness(Precision* buffer) override = 0;
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override = 0;

    MapMatrix  mass(Precision* buffer) override = 0;

    virtual const fem::quadrature::Quadrature& integration_scheme() const = 0;

    //
    StaticMatrix<N, 3> node_coords_global() {
        logging::error(this->_model_data != nullptr, "no model data assigned to element ", this->elem_id);
        logging::error(this->_model_data->positions != nullptr, "positions field not set in model data");
        const auto& positions = *this->_model_data->positions;
        StaticMatrix<N, 3> res {};
        for (Index i = 0; i < N; i++) {
            const Index row = static_cast<Index>(this->node_ids[i]);
            res.row(i) = positions.row_vec3(row).transpose();
        }
        return res;
    }

    Dim        n_integration_points() override {
        return integration_scheme().count();
    }

    ElDofs     dofs() override {
        return ElDofs{true, true, true, true, true, true};
    }
    Dim        dimensions() override {
        return 3;
    }
    Dim        n_nodes() override {
        return N;
    }
    ID*        nodes() override {
        return node_ids.data();
    }

    void compute_internal_force_nonlinear(Field& node_forces,
                                          const Field& displacement,
                                          const Field& ip_stress,
                                          int ip_offset) override {
        (void) node_forces;
        (void) displacement;
        (void) ip_stress;
        (void) ip_offset;
        logging::error(false, "ShellElement: compute_internal_force_nonlinear is not implemented yet for element ", this->elem_id);
    };

    void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) override {
        (void) node_loads;
        (void) node_temp;
        (void) ref_temp;
    };

    void compute_stress_strain(Field* strain,
                               Field* stress,
                               const Field& displacement,
                               const RowMatrix& rst,
                               int offset,
                               bool use_green_lagrange_nl) override {
        (void) strain;
        (void) stress;
        (void) displacement;
        (void) rst;
        (void) offset;
        (void) use_green_lagrange_nl;
        logging::error(false, "ShellElement: compute_stress_strain is not implemented yet for element ", this->elem_id);
    }

    void compute_stress_state(Field& stress_state,
                              const Field& displacement,
                              int offset,
                              bool use_green_lagrange_nl) override {
        RowMatrix rst = stress_strain_ip_rst();
        if (rst.rows() == 0) {
            return;
        }
        compute_stress_strain(nullptr, &stress_state, displacement, rst, offset, use_green_lagrange_nl);
    }

    bool compute_shear_flow(Field& shear_flow,
                            const Field& displacement,
                            int offset) override {
        (void) shear_flow;
        (void) displacement;
        (void) offset;
        return false;
    }

    bool compute_beam_section_forces(Field& section_forces,
                                     const Field& displacement,
                                     int offset) override {
        (void) section_forces;
        (void) displacement;
        (void) offset;
        return false;
    }

    bool compute_shell_section_forces(Field& section_forces,
                                      Field& contribution_count,
                                      const Field& displacement) override {
        (void) section_forces;
        (void) contribution_count;
        (void) displacement;
        return false;
    }

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
    StaticMatrix<N, K> nodal_data(const Field& full_data, Index offset = 0, Index stride = 1) {
        StaticMatrix<N, K> res {};
        runtime_assert(
            full_data.components >= offset + stride * (K - 1) + 1,
            "cannot extract this many elements from the data"
        );

        for (Dim m = 0; m < N; m++) {
            for (Dim j = 0; j < K; j++) {
                Index n = j * stride + offset;
                res(m, j) = full_data(static_cast<Index>(node_ids[m]), n);
            }
        }

        return res;
    }

    void compute_compliance(Field& displacement, Field& result) override {
        // Elementsteifigkeit (global gedreht) holen
        Precision buffer[6 * N * 6 * N];
        MapMatrix Ke = stiffness(buffer); // 6N × 6N

        // Element-Verschiebungsvektor (global) aufbauen: [ux,uy,uz,rx,ry,rz] je Knoten
        // nodal_data<6>(...) liefert dir genau diese 6 DOFs pro Knoten in globalen Achsen
        StaticMatrix<6, N> u_mat = StaticMatrix<6, N>(this->nodal_data<6>(displacement).transpose());
        Eigen::Map<StaticVector<6 * N>> ue(u_mat.data(), 6 * N);

        // Compliance-Beitrag des Elements:
        // Bei linearem Gleichgewicht gilt f_e = K_e u_e ⇒ u^T K u = u^T f (klassische Compliance).
        const Precision Ce = ue.dot(Ke * ue);

        // Speichern (überschreiben oder aufsummieren – je nach gewünschter Semantik)
        result(this->elem_id, 0) = Ce;
    }

    void compute_compliance_angle_derivative(Field& displacement, Field& result) override {
        (void) displacement;
        (void) result;
    };
    protected:
};
}

#endif //SHELL_H
