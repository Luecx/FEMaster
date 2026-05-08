#pragma once
#include <Eigen/Sparse>
#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#include <mkl.h>
#endif

template<class Scalar>
class PardisoShiftInvertOp {
    public:
    using SpMat = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, int>;
    using Vec   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    PardisoShiftInvertOp(const SpMat& A, const SpMat& B, Scalar sigma, bool assume_spd=false)
        : A_(A), B_(B), sigma_(sigma), assume_spd_(assume_spd)
    {
        // Build M = A - sigma B
        ensure_compressed(const_cast<SpMat&>(A_));
        ensure_compressed(const_cast<SpMat&>(B_));
        M_ = A_;
        if(sigma_ != Scalar(0)) {
            SpMat tmp = B_;
            tmp *= sigma_;
            M_ -= tmp;
        }

#ifdef USE_MKL
        mkl_set_num_threads(std::max(1, global_config.max_threads));
        if (assume_spd_) {
            ldlt_.analyzePattern(M_);
            ldlt_.factorize(M_);
            if (ldlt_.info() != Eigen::Success) {
                // fall back if not SPD in practice
                assume_spd_ = false;
            }
        }
        if (!assume_spd_) {
            lu_.analyzePattern(M_);
            lu_.factorize(M_);
            if (lu_.info() != Eigen::Success) throw std::runtime_error("PARDISO factorization failed");
        }
#else
        // Fallback: Eigen SparseLU (robust) / SimplicialLDLT (fast SPD)
        if (assume_spd_) {
            llt_.analyzePattern(M_); llt_.factorize(M_);
            if (llt_.info() != Eigen::Success) { assume_spd_ = false; }
        }
        if (!assume_spd_) {
            elu_.analyzePattern(M_); elu_.factorize(M_);
            if (elu_.info() != Eigen::Success) throw std::runtime_error("Eigen SparseLU factorization failed");
        }
#endif
    }

    int rows() const { return A_.rows(); }
    int cols() const { return A_.cols(); }

    // y <- (A - sigma B)^{-1} x
    void perform_op(const Scalar* x_in, Scalar* y_out) const {
        Eigen::Map<const Vec> x(x_in, rows());
        Eigen::Map<Vec>       y(y_out, rows());

#ifdef USE_MKL
        if (assume_spd_) {
            y = ldlt_.solve(x);
            if (ldlt_.info() != Eigen::Success) throw std::runtime_error("PARDISO LDLT solve failed");
        } else {
            y = lu_.solve(x);
            if (lu_.info() != Eigen::Success) throw std::runtime_error("PARDISO LU solve failed");
        }
#else
        if (assume_spd_) {
            y = llt_.solve(x);
            if (llt_.info() != Eigen::Success) throw std::runtime_error("LLT solve failed");
        } else {
            y = elu_.solve(x);
            if (elu_.info() != Eigen::Success) throw std::runtime_error("SparseLU solve failed");
        }
#endif
    }

    // Optional: fast re-factor when only sigma changes and sparsity pattern is same
    void refactor_with_shift(Scalar new_sigma) {
        if(new_sigma == sigma_) return;
        sigma_ = new_sigma;
        M_ = A_;
        if (sigma_ != Scalar(0)) {
            SpMat tmp = B_;
            tmp *= sigma_;
            M_ -= tmp;
        }
#ifdef USE_MKL
        if (assume_spd_) { ldlt_.factorize(M_); if(ldlt_.info()!=Eigen::Success) throw std::runtime_error("LDLT refactor failed"); }
        else              { lu_.factorize(M_);   if(lu_.info()!=Eigen::Success)   throw std::runtime_error("LU refactor failed"); }
#else
        if (assume_spd_) { llt_.factorize(M_); if(llt_.info()!=Eigen::Success) throw std::runtime_error("LLT refactor failed"); }
        else              { elu_.factorize(M_); if(elu_.info()!=Eigen::Success) throw std::runtime_error("SparseLU refactor failed"); }
#endif
    }

    private:
    const SpMat& A_;
    const SpMat& B_;
    Scalar       sigma_;
    bool         assume_spd_;
    SpMat        M_;

#ifdef USE_MKL
    mutable Eigen::PardisoLDLT<SpMat> ldlt_; // when SPD
    mutable Eigen::PardisoLU<SpMat>   lu_;   // robust default
#else
    mutable Eigen::SimplicialLDLT<SpMat> llt_;
    mutable Eigen::SparseLU<SpMat>       elu_;
#endif

    static void ensure_compressed(SpMat& M){ if(!M.isCompressed()) M.makeCompressed(); }
};
