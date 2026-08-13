/**
 * @file c3d8r_assembly.ipp
 * @brief Adds the physical hourglass response to C3D8R assembly.
 *
 * @author Finn Eggers
 * @date 13.08.2026
 */

void C3D8R::assemble_local_force(Field& node_forces, const Vector24& local_force) {
    logging::error(node_forces.domain == FieldDomain::NODE,
        "C3D8R: internal force output must use NODE domain");
    logging::error(node_forces.components >= D,
        "C3D8R: internal force output requires at least three components");

    for (Index a = 0; a < N; ++a) {
        const Index node_id = static_cast<Index>(node_ids[a]);
        for (Dim d = 0; d < D; ++d) {
            node_forces(node_id, d) += local_force(D * a + d);
        }
    }
}

MapMatrix C3D8R::stiffness(Precision* buffer) {
    MapMatrix mapped{buffer, ndof, ndof};
    C3D8::stiffness(buffer);
    mapped += hourglass_stiffness();
    mapped = Precision(0.5) * (mapped + mapped.transpose());
    return mapped;
}

MapMatrix C3D8R::stiffness_tangent(Precision* buffer,
                                   Field& ip_stress_state,
                                   NodeData& nodal_forces,
                                   const Field& displacement) {
    MapMatrix mapped = SolidElement<N>::stiffness_tangent(
        buffer, ip_stress_state, nodal_forces, displacement);

    Vector24 force = Vector24::Zero();
    Matrix24 tangent = Matrix24::Zero();
    hourglass_response(force, &tangent);

    mapped += tangent;
    mapped = Precision(0.5) * (mapped + mapped.transpose());
    assemble_local_force(nodal_forces, force);
    return mapped;
}

void C3D8R::compute_internal_force_nonlinear(Field& node_forces,
                                             const Field& ip_stress) {
    C3D8::compute_internal_force_nonlinear(node_forces, ip_stress);

    Vector24 force = Vector24::Zero();
    hourglass_response(force, nullptr);
    assemble_local_force(node_forces, force);
}
