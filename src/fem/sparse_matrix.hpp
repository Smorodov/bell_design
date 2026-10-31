#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace BellOptimizer::detail {

struct SparseTriplet {
    int row = 0;
    int column = 0;
    double value = 0.0;
};

// Dependency-free CSR matrix. Duplicate triplets are summed deterministically;
// exact zeros are omitted and no numerical threshold is applied.
struct SparseMatrix {
    int n = 0;
    std::vector<std::size_t> row_offsets;
    std::vector<int> columns;
    std::vector<double> values;

    static SparseMatrix from_triplets(int dimension, std::vector<SparseTriplet> entries) {
        if (dimension < 0) throw std::invalid_argument("negative sparse matrix dimension");
        for (const auto& entry : entries) {
            if (entry.row < 0 || entry.row >= dimension || entry.column < 0 ||
                entry.column >= dimension)
                throw std::out_of_range("sparse triplet index out of range");
        }
        std::sort(entries.begin(), entries.end(), [](const SparseTriplet& a, const SparseTriplet& b) {
            return a.row < b.row || (a.row == b.row && a.column < b.column);
        });
        SparseMatrix result;
        result.n = dimension;
        result.row_offsets.assign(static_cast<std::size_t>(dimension) + 1u, 0u);
        std::size_t i = 0;
        while (i < entries.size()) {
            const int row = entries[i].row;
            const int column = entries[i].column;
            double sum = 0.0;
            do { sum += entries[i].value; ++i; }
            while (i < entries.size() && entries[i].row == row && entries[i].column == column);
            if (sum != 0.0) {
                result.columns.push_back(column);
                result.values.push_back(sum);
                ++result.row_offsets[static_cast<std::size_t>(row) + 1u];
            }
        }
        for (std::size_t row = 1; row < result.row_offsets.size(); ++row)
            result.row_offsets[row] += result.row_offsets[row - 1u];
        return result;
    }

    double get(int row, int column) const {
        if (row < 0 || row >= n || column < 0 || column >= n)
            throw std::out_of_range("sparse matrix index out of range");
        const std::size_t begin = row_offsets[static_cast<std::size_t>(row)];
        const std::size_t end = row_offsets[static_cast<std::size_t>(row) + 1u];
        const auto first = columns.begin() + static_cast<std::ptrdiff_t>(begin);
        const auto last = columns.begin() + static_cast<std::ptrdiff_t>(end);
        const auto it = std::lower_bound(first, last, column);
        return it == last || *it != column ? 0.0 : values[static_cast<std::size_t>(it - columns.begin())];
    }

    std::vector<double> multiply(const std::vector<double>& x) const {
        if (x.size() != static_cast<std::size_t>(n))
            throw std::invalid_argument("sparse matrix/vector dimension mismatch");
        std::vector<double> y(static_cast<std::size_t>(n), 0.0);
        for (int row = 0; row < n; ++row) {
            double sum = 0.0;
            for (std::size_t k = row_offsets[static_cast<std::size_t>(row)];
                 k < row_offsets[static_cast<std::size_t>(row) + 1u]; ++k)
                sum += values[k] * x[static_cast<std::size_t>(columns[k])];
            y[static_cast<std::size_t>(row)] = sum;
        }
        return y;
    }

    std::size_t storage_bytes() const {
        return row_offsets.size() * sizeof(std::size_t) + columns.size() * sizeof(int) +
               values.size() * sizeof(double);
    }


};


struct SparseOrderingMetrics {
    std::size_t max_bandwidth = 0;
    std::size_t profile = 0;
};

inline SparseOrderingMetrics structural_metrics(const SparseMatrix& A) {
    SparseOrderingMetrics metrics;
    for (int row = 0; row < A.n; ++row) {
        std::size_t first = static_cast<std::size_t>(row);
        bool offdiag = false;
        for (std::size_t p = A.row_offsets[static_cast<std::size_t>(row)];
             p < A.row_offsets[static_cast<std::size_t>(row) + 1u]; ++p) {
            const std::size_t col = static_cast<std::size_t>(A.columns[p]);
            if (col == static_cast<std::size_t>(row)) continue;
            offdiag = true;
            metrics.max_bandwidth = std::max(metrics.max_bandwidth,
                col > static_cast<std::size_t>(row)
                    ? col - static_cast<std::size_t>(row)
                    : static_cast<std::size_t>(row) - col);
            first = std::min(first, col);
        }
        if (offdiag) metrics.profile += static_cast<std::size_t>(row) - first;
    }
    return metrics;
}

// Reverse Cuthill-McKee ordering of the symmetric structural graph. The
// returned vector is new-index -> old-index. Only the sparsity pattern is
// inspected; coefficient magnitudes are never used.
inline SparseMatrix structural_union(const SparseMatrix& a, const SparseMatrix& b) {
    if (a.n != b.n) throw std::invalid_argument("sparse structural union dimension mismatch");
    std::vector<SparseTriplet> entries;
    entries.reserve(a.values.size() + b.values.size());
    for (int i = 0; i < a.n; ++i) {
        for (std::size_t p = a.row_offsets[static_cast<std::size_t>(i)];
             p < a.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p)
            entries.push_back({i, a.columns[p], 1.0});
        for (std::size_t p = b.row_offsets[static_cast<std::size_t>(i)];
             p < b.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p)
            entries.push_back({i, b.columns[p], 1.0});
    }
    return SparseMatrix::from_triplets(a.n, std::move(entries));
}

inline std::vector<int> reverse_cuthill_mckee_order(const SparseMatrix& A) {
    if (A.n < 0) throw std::invalid_argument("negative sparse matrix dimension");
    const int n = A.n;
    std::vector<std::vector<int>> adjacency(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        for (std::size_t p = A.row_offsets[static_cast<std::size_t>(i)];
             p < A.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p) {
            const int j = A.columns[p];
            if (j != i) adjacency[static_cast<std::size_t>(i)].push_back(j);
        }
    }
    // Make the graph explicitly undirected. This also tolerates a diagnostic
    // input whose structural pattern is not exactly symmetric.
    for (int i = 0; i < n; ++i) {
        for (int j : adjacency[static_cast<std::size_t>(i)])
            adjacency[static_cast<std::size_t>(j)].push_back(i);
    }
    for (auto& neighbors : adjacency) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    std::vector<int> degree(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) degree[static_cast<std::size_t>(i)] =
        static_cast<int>(adjacency[static_cast<std::size_t>(i)].size());

    std::vector<char> visited(static_cast<std::size_t>(n), 0);
    std::vector<int> cm_order;
    cm_order.reserve(static_cast<std::size_t>(n));
    while (static_cast<int>(cm_order.size()) < n) {
        int start = -1;
        for (int i = 0; i < n; ++i) {
            if (!visited[static_cast<std::size_t>(i)] &&
                (start < 0 || degree[static_cast<std::size_t>(i)] < degree[static_cast<std::size_t>(start)]))
                start = i;
        }
        if (start < 0) break;
        std::vector<int> frontier{start};
        visited[static_cast<std::size_t>(start)] = 1;
        while (!frontier.empty()) {
            std::sort(frontier.begin(), frontier.end(), [&](int a, int b) {
                if (degree[static_cast<std::size_t>(a)] != degree[static_cast<std::size_t>(b)])
                    return degree[static_cast<std::size_t>(a)] < degree[static_cast<std::size_t>(b)];
                return a < b;
            });
            std::vector<int> next;
            for (int u : frontier) {
                cm_order.push_back(u);
                for (int v : adjacency[static_cast<std::size_t>(u)]) {
                    if (!visited[static_cast<std::size_t>(v)]) {
                        visited[static_cast<std::size_t>(v)] = 1;
                        next.push_back(v);
                    }
                }
            }
            frontier.swap(next);
        }
    }
    std::reverse(cm_order.begin(), cm_order.end());
    return cm_order;
}

// Symmetric permutation A' = P A P^T, where order[new] = old. The numerical
// matrix is unchanged up to relabeling of degrees of freedom.
inline SparseMatrix symmetric_permute(const SparseMatrix& A, const std::vector<int>& order) {
    if (static_cast<int>(order.size()) != A.n)
        throw std::invalid_argument("sparse permutation dimension mismatch");
    std::vector<int> inverse(static_cast<std::size_t>(A.n), -1);
    for (int new_index = 0; new_index < A.n; ++new_index) {
        const int old_index = order[static_cast<std::size_t>(new_index)];
        if (old_index < 0 || old_index >= A.n || inverse[static_cast<std::size_t>(old_index)] >= 0)
            throw std::invalid_argument("invalid sparse permutation");
        inverse[static_cast<std::size_t>(old_index)] = new_index;
    }
    std::vector<SparseTriplet> entries;
    entries.reserve(A.values.size());
    for (int old_row = 0; old_row < A.n; ++old_row) {
        const int new_row = inverse[static_cast<std::size_t>(old_row)];
        for (std::size_t p = A.row_offsets[static_cast<std::size_t>(old_row)];
             p < A.row_offsets[static_cast<std::size_t>(old_row) + 1u]; ++p) {
            entries.push_back({new_row, inverse[static_cast<std::size_t>(A.columns[p])], A.values[p]});
        }
    }
    return SparseMatrix::from_triplets(A.n, std::move(entries));
}

struct SymmetricPruneReport {
    SparseMatrix matrix;
    std::size_t removed_pairs = 0;
    double removed_frobenius_norm = 0.0;
    double original_frobenius_norm = 0.0;
};

// Diagnostic-only symmetric pruning. The diagonal is never removed. An
// off-diagonal pair (i,j)/(j,i) is removed only when both entries exist,
// are numerically symmetric, and |a_ij| <= tolerance*sqrt(|a_ii*a_jj|).
// This is deliberately not used by production FEM assembly or the solver.
inline SymmetricPruneReport prune_symmetric_offdiagonals(const SparseMatrix& input,
                                                         double tolerance,
                                                         double symmetry_tolerance = 1e-12) {
    if (!(tolerance >= 0.0) || !std::isfinite(tolerance) ||
        !(symmetry_tolerance >= 0.0) || !std::isfinite(symmetry_tolerance))
        throw std::invalid_argument("invalid sparse pruning tolerance");
    SymmetricPruneReport report;
    for (double value : input.values) report.original_frobenius_norm += value * value;
    report.original_frobenius_norm = std::sqrt(report.original_frobenius_norm);
    std::vector<SparseTriplet> kept;
    kept.reserve(input.values.size());
    for (int i = 0; i < input.n; ++i) {
        const double di = input.get(i, i);
        for (std::size_t p = input.row_offsets[static_cast<std::size_t>(i)];
             p < input.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p) {
            const int j = input.columns[p];
            const double aij = input.values[p];
            if (j < i) continue;
            if (j == i) {
                kept.push_back({i, j, aij});
                continue;
            }
            const double aji = input.get(j, i);
            const double pair_scale = std::max({std::abs(aij), std::abs(aji),
                                                std::numeric_limits<double>::min()});
            const double symmetry_error = std::abs(aij - aji) / pair_scale;
            const double normalized_scale = std::sqrt(std::abs(di * input.get(j, j)));
            const bool removable = symmetry_error <= symmetry_tolerance &&
                normalized_scale > 0.0 && std::abs(aij) <= tolerance * normalized_scale &&
                std::abs(aji) <= tolerance * normalized_scale;
            if (removable) {
                ++report.removed_pairs;
                report.removed_frobenius_norm += aij * aij + aji * aji;
            } else {
                kept.push_back({i, j, aij});
                kept.push_back({j, i, aji});
            }
        }
    }
    report.removed_frobenius_norm = std::sqrt(report.removed_frobenius_norm);
    report.matrix = SparseMatrix::from_triplets(input.n, std::move(kept));
    return report;
}


} // namespace BellOptimizer::detail
