#include "cb_scan_simple.h"

namespace fem { namespace constraint {

std::vector<SimpleRow>
scan_single_nnz_rows(const SparseMatrix& C,
                     const DynamicVector& d,
                     bool homogeneous)
{
    const int m = C.rows();
    const bool has_d = (d.size() == m) && !homogeneous;

    // Count nnz per row + remember last (col, val)
    std::vector<int> nnz_row(m, 0);
    std::vector<int> last_col(m, -1);
    std::vector<Precision> last_val(m, Precision(0));

    for (int j = 0; j < C.outerSize(); ++j) {
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            const int i = it.row();
            nnz_row[i]     += 1;
            last_col[i]     = j;
            last_val[i]     = it.value();
        }
    }

    std::vector<SimpleRow> singles;
    singles.reserve(m / 8 + 1);
    for (int i = 0; i < m; ++i) {
        if (nnz_row[i] == 1 && last_col[i] >= 0 && last_val[i] != Precision(0)) {
            singles.push_back(SimpleRow{ i, last_col[i], last_val[i], has_d ? d[i] : Precision(0) });
        }
    }
    return singles;
}

}} // namespace fem::constraint
