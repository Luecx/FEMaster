// /******************************************************************************
//  * @file demo_nullspace_10k_explicit.cpp
//  * @brief 10k x 10k SPD system with ~1000 constraints; EXPLICIT T, assemble
//  *        A = Tᵀ K T and B = Tᵀ K_g T; direct solve on A; timings + nnz.
//  ******************************************************************************/
//
// #include <iostream>
// #include <random>
// #include <chrono>
// #include <Eigen/Sparse>
// #include <Eigen/Core>
//
// #include "constraints/constraint_transformer.h"
// #include "constraints/equation.h"
// #include "core/types_eig.h"
// #include "core/types_cls.h"
// #include "core/logging.h"
//
// using namespace fem;
// using namespace fem::constraint;
// using Clock = std::chrono::steady_clock;
//
// // --- knobs -------------------------------------------------------------------
// static constexpr int N_DOF        = 10000;   // total DOFs
// static constexpr int N_CONS       = 1000;    // number of constraint rows
// static constexpr int N_PER_NODE   = 6;       // (node,dof) slots; only dof 0 active
// static constexpr bool BUILD_B     = true;    // also assemble B = Tᵀ K_g T
// static constexpr unsigned SEED_K  = 2025;
// static constexpr unsigned SEED_F  = 777;
// static constexpr unsigned SEED_EQ = 12345;
//
// // --- helpers -----------------------------------------------------------------
// static SparseMatrix make_1d_laplacian(int n) {
//     SparseMatrix K(n, n);
//     std::vector<Eigen::Triplet<Precision>> trips;
//     trips.reserve(3LL * n);
//     for (int i = 0; i < n; ++i) {
//         trips.emplace_back(i, i, Precision(2.0));
//         if (i > 0)     trips.emplace_back(i, i - 1, Precision(-1.0));
//         if (i+1 < n)   trips.emplace_back(i, i + 1, Precision(-1.0));
//     }
//     K.setFromTriplets(trips.begin(), trips.end());
//     K.makeCompressed();
//     return K;
// }
//
// static SparseMatrix make_geom_like(const SparseMatrix& K) {
//     // Example geometric-like SPD-ish matrix: scaled copy of K + tiny diagonal
//     SparseMatrix Kg = K;
//     // add small diagonal to be safe
//     Eigen::VectorXd d = Eigen::VectorXd::Constant(K.rows(), 1e-6);
//     SparseMatrix D(K.rows(), K.cols());
//     std::vector<Eigen::Triplet<Precision>> trips;
//     trips.reserve(K.rows());
//     for (int i = 0; i < K.rows(); ++i) trips.emplace_back(i, i, d[i]);
//     D.setFromTriplets(trips.begin(), trips.end());
//     D.makeCompressed();
//     Kg = Kg + D;
//     return Kg;
// }
//
// int main() {
//     logging::info(true, "================== DEMO: 10k SPD with ~1000 constraints (explicit T) ==================");
//
//     const int n  = N_DOF;
//     const int nd = N_PER_NODE;
//
//     auto t0 = Clock::now();
//     SparseMatrix K = make_1d_laplacian(n);
//     auto t1 = Clock::now();
//
//     DynamicVector f(n);
//     {
//         std::mt19937 gen(SEED_F);
//         std::uniform_real_distribution<Precision> uni(-1.0, 1.0);
//         for (int i = 0; i < n; ++i) f[i] = uni(gen);
//     }
//     auto t2 = Clock::now();
//
//     // ---- Constraints: ~1/2 supports, ~1/2 ties (i != j) --------------------
//     Equations eqs;
//     eqs.reserve(N_CONS);
//     {
//         std::mt19937 gen_idx(SEED_EQ);
//         std::uniform_int_distribution<int> pick(0, n-1);
//
//         for (int k = 0; k < N_CONS; ++k) {
//             if (k % 2 == 0) {
//                 int i = pick(gen_idx);
//                 eqs.emplace_back( Equation{ { {i, 0,  1.0} }, 0.0 } );
//             } else {
//                 int i = pick(gen_idx);
//                 int j = pick(gen_idx);
//                 while (j == i) j = pick(gen_idx);
//                 eqs.emplace_back( Equation{ { {i, 0,  1.0}, {j, 0, -2.0} }, 0.0 } );
//             }
//         }
//     }
//     auto t3 = Clock::now();
//
//     // SystemDofIds: activate only dof 0 -> global index = node id; others -1
//     SystemDofIds dofmap(n, nd);
//     for (int i = 0; i < n; ++i)
//         for (int j = 0; j < nd; ++j)
//             dofmap(i, j) = (j == 0 ? i : -1);
//
//     ConstraintTransformer::BuildOptions opt;
//     opt.set.scale_columns     = false; // stays sparse, as requested
//     opt.set.scale_rows        = false;
//     opt.set.zero_row_drop_tol = 0;
//     opt.builder.rank_tol_rel  = 1e-12;
//     opt.builder.feas_tol_rel  = 1e-10;
//
//     auto t4 = Clock::now();
//     ConstraintTransformer CT(eqs, dofmap, n, opt);
//     auto t5 = Clock::now();
//
//     const auto& rep = CT.report();
//     const auto& map = CT.map();
//     logging::info(true, "[ConstraintSet] m=", rep.m, " n=", rep.n, " rank=", rep.rank,
//                   " nm=", map.n_master(), " homogeneous=", CT.homogeneous());
//
//     // ---- form T explicitly (already in map.T()) -----------------------------
//     auto t5b = Clock::now();
//     const SparseMatrix& T = map.T(); // explicit T
//     auto t5c = Clock::now();
//
//     // ---- assemble A = Tᵀ K T and b = Tᵀ(f - K u_p) -------------------------
//     auto t6 = Clock::now();
//     SparseMatrix A = map.assemble_A(K);
//     DynamicVector b = map.assemble_b(K, f);
//     auto t7 = Clock::now();
//
//     // ---- optional: build Kg and B = Tᵀ K_g T --------------------------------
//     SparseMatrix Kg, B;
//     auto t7b = Clock::now();
//     if (BUILD_B) {
//         Kg = make_geom_like(K);
//         B  = map.assemble_B(Kg);
//     }
//     auto t7c = Clock::now();
//
//     // ---- direct solve on A ---------------------------------------------------
//     auto t8 = Clock::now();
//     Eigen::SimplicialLDLT<SparseMatrix> ldlt;
//     ldlt.compute(A);
//     if (ldlt.info() != Eigen::Success) {
//         logging::error(false, "LDLT factorization failed on A (nnz=", (Index)A.nonZeros(), ").");
//         return 1;
//     }
//     DynamicVector q = ldlt.solve(b);
//     if (ldlt.info() != Eigen::Success) {
//         logging::error(false, "LDLT solve failed.");
//         return 1;
//     }
//     auto t9 = Clock::now();
//
//     // ---- recover/check -------------------------------------------------------
//     auto t10 = Clock::now();
//     DynamicVector u = CT.recover_u(q);
//     DynamicVector r = K * u - f;
//
//     double u_norm   = u.norm();
//     double r_norm   = r.norm();
//     double ttr_norm = map.apply_Tt(r).norm();
//     double cu_norm  = (CT.set().C * u - CT.set().d).norm();
//     auto t11 = Clock::now();
//
//     auto ms = [](Clock::time_point a, Clock::time_point b){
//         return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
//     };
//
//     std::cout << "---- sizes/nnz -----------------------------------------\n";
//     std::cout << "n=" << n
//               << "  rank=" << rep.rank
//               << "  nm=" << map.n_master()
//               << "  nnz(K)=" << K.nonZeros()
//               << "  nnz(T)=" << T.nonZeros()
//               << "  nnz(A)=" << A.nonZeros();
//     if (BUILD_B) std::cout << "  nnz(B)=" << B.nonZeros();
//     std::cout << "\n";
//
//     std::cout << "---- residuals -----------------------------------------\n";
//     std::cout << "[INFO] ||u||2             = " << u_norm   << "\n";
//     std::cout << "[INFO] ||K u - f||2       = " << r_norm   << "\n";
//     std::cout << ((ttr_norm > 1e-8) ? "[WARN] " : "[OK]   ")
//               << "||T^T (K u - f)||2 = " << ttr_norm << "\n";
//     std::cout << ((cu_norm > 1e-9) ? "[WARN] " : "[OK]   ")
//               << "||C u - d||2       = " << cu_norm  << "\n";
//
//     std::cout << "---- timings (ms) -------------------------------------\n";
//     std::cout << "build K:                    " << ms(t0,t1)  << "\n";
//     std::cout << "build f:                    " << ms(t1,t2)  << "\n";
//     std::cout << "build constraints (eqs):    " << ms(t2,t3)  << "\n";
//     std::cout << "QR map build (CT):          " << ms(t4,t5)  << "\n";
//     std::cout << "T materialize (get+cmp):    " << ms(t5b,t5c)<< "\n";
//     std::cout << "assemble A & b:             " << ms(t6,t7)  << "\n";
//     if (BUILD_B) {
//         std::cout << "build Kg + assemble B:      " << ms(t7b,t7c) << "\n";
//     }
//     std::cout << "direct solve (LDLT):        " << ms(t8,t9)  << "\n";
//     std::cout << "recover + checks:           " << ms(t10,t11)<< "\n";
//     std::cout << "total (all):                " << ms(t0,t11) << "\n";
//     std::cout << "-------------------------------------------------------\n";
//
//     return 0;
// }

#include "reader/reader.h"
#include "core/logging.h"

#include <argparse/argparse.hpp>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <filesystem>


#include "math/interpolate.h"

int main(int argc, char* argv[]) {
    // Set up argparse
    argparse::ArgumentParser program("FEM Solver");

    program.add_argument("input_file")
        .help("Path to the input file (.inp is optional).");

    program.add_argument("--ncpus")
        .default_value(1)
        .scan<'i', int>()
        .help("Number of CPUs to use (default: 1)");

    // Parse arguments
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    // Access parsed arguments
    std::string input_file = program.get<std::string>("input_file");
    int ncpus = program.get<int>("--ncpus");

    // Ensure the input file has ".inp" extension
    if (input_file.find(".inp") == std::string::npos) {
        input_file += ".inp";
    }

    // Check if input file exists
    if (!std::filesystem::exists(input_file)) {
        std::cerr << "Error: Input file '" << input_file << "' does not exist." << std::endl;
        return 1;
    }

    // Create the output file by replacing ".inp" with ".res"
    std::string output_file = input_file;
    size_t pos = output_file.find(".inp");
    if (pos != std::string::npos) {
        output_file.replace(pos, 4, ".res");
    }

    // Logging input and output file information
    fem::logging::info(true, "");
    fem::logging::info(true, "Input file : ", input_file);
    fem::logging::info(true, "Output file: ", output_file);
    fem::logging::info(true, "CPU(s)     : ", ncpus);
    fem::logging::info(true, "");

    // Store number of CPUs in config
    fem::global_config.max_threads = ncpus;

    // Read the input file using the reader
    fem::reader::Reader reader{input_file, output_file};
    reader.read();

    return 0;
}

