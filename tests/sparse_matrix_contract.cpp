#include "fem/sparse_matrix.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    using BellOptimizer::detail::SparseMatrix;
    using BellOptimizer::detail::SparseTriplet;
    const std::vector<SparseTriplet> entries = {
        {0, 0, 2.0}, {0, 1, 1.5}, {0, 1, 0.5},
        {1, 0, 2.0}, {1, 1, 3.0}, {2, 2, 4.0}, {2, 1, 0.0}
    };
    const SparseMatrix a = SparseMatrix::from_triplets(3, entries);
    if (a.n != 3 || a.values.size() != 5u) return 1;
    if (std::abs(a.get(0, 1) - 2.0) > 1e-15 || a.get(2, 1) != 0.0) return 2;
    const std::vector<double> x{1.0, 2.0, -1.0};
    const auto y = a.multiply(x);
    const std::vector<double> reference{6.0, 8.0, -4.0};
    for (std::size_t i = 0; i < y.size(); ++i)
        if (std::abs(y[i] - reference[i]) > 1e-14) return 3;
    bool dimension_error = false;
    try { (void)a.multiply(std::vector<double>{1.0}); }
    catch (const std::invalid_argument&) { dimension_error = true; }
    if (!dimension_error) return 4;
    bool index_error = false;
    try { (void)a.get(3, 0); }
    catch (const std::out_of_range&) { index_error = true; }
    if (!index_error) return 5;
    std::cout << "sparse_matrix_contract=PASS\n"
              << "dimension=" << a.n << "\nnonzeros=" << a.values.size()
              << "\nstorage_bytes=" << a.storage_bytes() << '\n';
    return 0;
}
