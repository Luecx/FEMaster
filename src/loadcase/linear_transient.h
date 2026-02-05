#pragma once
/**
 * @file transient.h
 * @brief Linear transient analysis (implicit Newmark-β) using affine null-space constraints.
 *
 * Pipeline (homogeneous supports assumed, u_p = 0):
 *   (0) assign sections/materials
 *   (1) build unconstrained index (node x 6 -> active dof or -1)
 *   (2) build constraint equations (supports/ties) -> C u = d
 *   (3) assemble global load matrix (node x 6)
 *   (4) assemble active K, M
 *   (5) reduce loads to active RHS f
 *   (6) build ConstraintTransformer -> T, (u_p)
 *   (7) reduce operators:  A = Tᵀ K T,  Mr = Tᵀ M T
 *   (8) Rayleigh damping in reduced space: Cr = α Mr + β A
 *   (9) reduced RHS (time-independent): fr = Tᵀ (f - K u_p)  (here u_p=0 ⇒ fr=Tᵀ f)
 *  (10) run implicit Newmark-β on (Mr,Cr,A) with constant fr
 *  (11) expand snapshots to node×DOF and write with selected cadence
 */

#include "loadcase.h"
#include "../solve/newmark.h"     // NewmarkOpts, newmark_linear
#include "../damping/rayleigh.h"   // Rayleigh
#include "../core/types_eig.h"
#include "../data/field.h"         // model::Field

#include <vector>
#include <string>
#include <optional>

namespace fem { namespace loadcase {

struct Transient : public LoadCase {
    explicit Transient(ID id,
                       reader::Writer* writer,
                       model::Model* model);

    // User inputs
    std::vector<std::string> supps;   ///< supports/ties/couplings
    std::vector<std::string> loads;   ///< load collectors (time-independent)

    // Integrator parameters (fixed step)
    double dt      = 1e-3;
    double t_start = 0.0;
    double t_end   = 1.0;
    double beta    = 0.25;  ///< Newmark β
    double gamma   = 0.5;   ///< Newmark γ

    // Damping (analysis-specific)
    std::optional<damping::Rayleigh> rayleigh;

    // Writing cadence
    int    write_every_steps = 1;     ///< write every N steps (>=1)
    double write_every_time  = 0.0;   ///< if >0, overrides steps: write every Δt_write seconds

    // Solver selection
    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;

    // Initial conditions (optional)
    // Pointer to a node FIELD providing initial velocity (6 components: ux,uy,uz,rx,ry,rz).
    model::Field::Ptr initial_velocity = nullptr;

    // Setters (fluent)
    Transient& set_time(double dt_, double t_start_, double t_end_) { dt = dt_; t_end = t_end_; t_start = t_start_; return *this; }
    Transient& set_newmark(double beta_, double gamma_) { beta = beta_; gamma = gamma_; return *this; }
    Transient& set_damping(const damping::Rayleigh& r) { rayleigh = r; return *this; }
    Transient& clear_damping() { rayleigh.reset(); return *this; }
    Transient& set_write_every(int steps) { write_every_steps = std::max(1, steps); return *this; }
    Transient& set_write_every_time(double dt_write) { write_every_time = std::max(0.0, dt_write); return *this; }

    // Diagnostics output (optional)
    std::string stiffness_file;  ///< write "_K.mtx" (active) and "_A.mtx" (reduced)
    std::string mass_file;       ///< write "_M.mtx" (active) and "_Mr.mtx" (reduced)
    std::string damping_file;    ///< write "_C.mtx" (reduced)

    // Run analysis
    virtual void run() override;
};

}} // namespace fem::loadcase
