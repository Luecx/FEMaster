#include "cb_substitute.h"

namespace fem { namespace constraint {

SubstituteResult substitute_fixed(const SparseMatrix& C_in,
                                  const DynamicVector& d_in,
                                  const std::vector<SimpleRow>& singles)
{
    const int n = C_in.cols();
    const int m = C_in.rows();

    // Build fixed maps from single-NNZ rows
    std::vector<char> is_fixed_col(n, 0);
    std::vector<Precision> fixed_val(n, Precision(0));
    std::vector<char> keep_row(m, 1);

    for (const auto& sr : singles) {
        const Precision ufix = sr.d / sr.a;  // (if homogeneous d==0 ⇒ ufix==0)
        is_fixed_col[sr.col] = 1;
        fixed_val[sr.col]    = ufix;
        keep_row[sr.row]     = 0;            // drop this row from QR
    }

    // Ensure d is allocated for in-place subtraction
    DynamicVector d = d_in.size() ? d_in : DynamicVector::Zero(m);

    // d ← d − Σ_j C(:,j) * fixed_val[j]   (sparse column walks)
    for (int j = 0; j < n; ++j) {
        if (!is_fixed_col[j]) continue;
        const Precision ufix = fixed_val[j];
        if (ufix == Precision(0)) continue;
        for (SparseMatrix::InnerIterator it(C_in, j); it; ++it) {
            d[it.row()] -= it.value() * ufix;
        }
    }

    // Rebuild C without fixed columns
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(std::max<std::size_t>(1, C_in.nonZeros() - singles.size()));
    for (int j = 0; j < C_in.outerSize(); ++j) {
        if (is_fixed_col[j]) continue;
        for (SparseMatrix::InnerIterator it(C_in, j); it; ++it) {
            trips.emplace_back(it.row(), j, it.value());
        }
    }

    SparseMatrix C(m, n);
    C.setFromTriplets(trips.begin(), trips.end());
    C.makeCompressed();

    return { std::move(C), std::move(d), std::move(keep_row),
             std::move(is_fixed_col), std::move(fixed_val) };
}

}} // namespace fem::constraint
