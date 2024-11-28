////
//// Created by f_eggers on 14.10.2024.
////
//
//#ifndef BEAM_H
//#define BEAM_H
//
//#include "../../core/core.h"
//#include "../element/element_structural.h"
//
//#include <memory>
//
//namespace fem::model {
//
//template<Index N>
//struct BeamElement : public StructuralElement{
//
//        BeamElement(ID p_elem_id)
//                : StructuralElement(p_elem_id) {
//        }
//
//        virtual Precision  volume()                                                          {};
//        virtual MapMatrix  stiffness(Precision* buffer)                                       {};
//        virtual MapMatrix  mass(Precision* buffer)                                           {};
//        virtual void       compute_stress_strain_nodal(
//                                                           NodeData& displacement,
//                                                           NodeData& stress,
//                                                           NodeData& strain)                                        {};
//        virtual void       compute_stress_strain(
//                                                 NodeData& displacement,
//                                                 NodeData& stress,
//                                                 NodeData& strain,
//                                                 NodeData& xyz)                                                   {};
//        virtual void       apply_vload(NodeData& node_loads, Vec3 load) {};
//        virtual void       apply_tload(NodeData& node_loads, NodeData& node_temp, Precision ref_temp) {};
//        virtual void       compute_compliance(NodeData& displacement, ElementData& result) {};
//        virtual void       compute_compliance_angle_derivative(NodeData& displacement, ElementData& result) {};
//};
//
//}
//
//#endif //BEAM_H
