//
// Created by f_eggers on 11.12.2024.
//

#ifndef SHELL_H
#define SHELL_H

#include "../../material/stress.h"
#include "../../material/isotropic_elasticity.h"
#include "../../core/core.h"
#include "../element/element_structural.h"
#include "../../section/section_shell.h"

#include <memory>

namespace fem::model {
template<Index N>
struct ShellElement : StructuralElement {
    std::array<ID, N> node_ids;

    ShellElement(ID p_elem_id, std::array<ID, N> p_node_ids) : StructuralElement(p_elem_id), node_ids {p_node_ids} {}
    virtual ~ShellElement() override = default;


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
        if (!section->material) {
            logging::error(false, "Material not set for element ", this->elem_id);
        }
        return section->material;
    }
    material::IsotropicElasticity* get_elasticity() {
        auto mat_ptr = get_material();
        if (!mat_ptr->has_elasticity()) {
            logging::error(false, "Material has no elasticity assigned");
        }
        if (!mat_ptr->elasticity()->template as<material::IsotropicElasticity>()) {
            logging::error(false, "Material is not isotropic for element ", this->elem_id);
        }
        return mat_ptr->elasticity()->template as<material::IsotropicElasticity>();
    }

    // left out for childs
    virtual SurfacePtr surface(ID surface_id) override = 0;
    virtual Precision  volume() override = 0;
    virtual MapMatrix  stiffness(Precision* buffer) override = 0;
    virtual MapMatrix stiffness_geom(Precision* buffer, IPData& ip_stress, int ip_start_idx) override = 0;

    virtual MapMatrix  mass(Precision* buffer) override = 0;

    virtual const fem::quadrature::Quadrature& integration_scheme() const = 0;

    //
    StaticMatrix<N, 3> node_coords_global() {
        auto node_coords_system = this->_model_data->get(POSITION);
        StaticMatrix<N, 3> res {};
        for (Index i = 0; i < N; i++) {
            for (Index j = 0; j < 3; j++) {
                res(i, j) = node_coords_system(this->node_ids[i], j);
            }
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

    // not implemented
    void compute_stress_strain_nodal(NodeData& displacement, NodeData& stress, NodeData& strain) override {
        (void) displacement;
        (void) stress;
        (void) strain;
    };
    void compute_stress_strain(IPData& ip_stress, IPData& ip_strain, NodeData& displacement, int ip_offset) override {
        (void) ip_stress;
        (void) ip_strain;
        (void) displacement;
        (void) ip_offset;
    };
    void apply_vload(NodeData& node_loads, Vec3 load) override {
        (void) node_loads;
        (void) load;
    };
    void apply_tload(NodeData& node_loads, NodeData& node_temp, Precision ref_temp) override {
        (void) node_loads;
        (void) node_temp;
        (void) ref_temp;
    };

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
    StaticMatrix<N, K> nodal_data(const NodeData& full_data, Index offset = 0, Index stride = 1) {
        StaticMatrix<N, K> res {};
        runtime_assert(
            full_data.cols() >= offset + stride * (K - 1) + 1,
            "cannot extract this many elements from the data"
        );

        for (Dim m = 0; m < N; m++) {
            for (Dim j = 0; j < K; j++) {
                Index n = j * stride + offset;
                res(m, j) = full_data(node_ids[m], n);
            }
        }

        return res;
    }

    void compute_compliance(NodeData& displacement, ElementData& result) override {
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

    void compute_compliance_angle_derivative(NodeData& displacement, ElementData& result) override {
        (void) displacement;
        (void) result;
    };
    Stresses stress(NodeData& displacement, std::vector<Vec3>& rst) override {
        return {};
    };
    Strains  strain(NodeData& displacement, std::vector<Vec3>& rst) override {
        return {};
    };

    protected:
    Precision topo_stiffness_scale() const {
        if (this->_model_data && this->_model_data->elem_data.has(TOPO_STIFFNESS)) {
            return this->_model_data->elem_data.get(TOPO_STIFFNESS)(this->elem_id);
        }
        return Precision(1);
    }
};
}


#endif //SHELL_H
