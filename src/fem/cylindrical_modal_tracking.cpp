#include "cylindrical_modal_tracking.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace BellOptimizer::FEM {
namespace {

double dot_metric(const std::vector<double>& x, const std::vector<double>& y,
                  const detail::SparseMatrix& M) {
    if (x.size() != y.size() || x.size() != static_cast<std::size_t>(M.n))
        throw std::invalid_argument("modal tracking metric dimension mismatch");
    double sum = 0.0;
    for (int i = 0; i < M.n; ++i) {
        double row = 0.0;
        for (std::size_t p = M.row_offsets[static_cast<std::size_t>(i)];
             p < M.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p)
            row += M.values[p] * y[static_cast<std::size_t>(M.columns[p])];
        sum += x[static_cast<std::size_t>(i)] * row;
    }
    return sum;
}

std::vector<double> metric_orthogonalized(const std::vector<double>& v,
                                          const std::vector<std::vector<double>>& basis,
                                          const detail::SparseMatrix& M) {
    // Modified Gram-Schmidt is used in the mass inner product.  A second
    // projection pass is intentional: modal families can be nearly
    // degenerate, so the first pass can leave an O(epsilon) component along
    // an earlier vector even when the exact residual is orthogonal.  The
    // second pass removes that accumulated loss of orthogonality without
    // introducing a problem-specific tolerance.
    std::vector<double> out = v;
    for (int pass = 0; pass < 2; ++pass) {
        for (const auto& q : basis) {
            const double qq = dot_metric(q, q, M);
            if (!(qq > 0.0) || !std::isfinite(qq)) continue;
            const double qv = dot_metric(q, out, M);
            if (!std::isfinite(qv)) return {};
            const double alpha = qv / qq;
            if (!std::isfinite(alpha)) return {};
            for (std::size_t i = 0; i < out.size(); ++i) out[i] -= alpha * q[i];
        }
    }
    const double nn = dot_metric(out, out, M);
    if (!(nn > 0.0) || !std::isfinite(nn)) return {};
    const double inv = 1.0 / std::sqrt(nn);
    if (!std::isfinite(inv)) return {};
    for (double& value : out) value *= inv;
    return out;
}

double normalized_cross(const std::vector<double>& x, const std::vector<double>& y,
                        const detail::SparseMatrix& M, bool absolute_value = true) {
    const double xx = dot_metric(x, x, M);
    const double yy = dot_metric(y, y, M);
    if (!(xx > 0.0) || !(yy > 0.0)) return 0.0;
    const double value = dot_metric(x, y, M) / std::sqrt(xx * yy);
    return absolute_value ? std::abs(value) : value;
}

double smallest_singular_value_2x2(double a, double b, double c, double d) {
    // Singular values squared are the eigenvalues of C^T C. Avoid the
    // trace^2 - 4*det^2 discriminant: when the singular values are close,
    // that expression loses the small separation by cancellation. The same
    // eigenvalue gap can be formed robustly as hypot(A-D, 2B) for the
    // symmetric 2x2 matrix C^T C. Recover the smaller singular value from
    // |det(C)| / sigma_max, which is also stable for nearly equal singular
    // values.
    const double col0 = std::hypot(a, c);
    const double col1 = std::hypot(b, d);
    const double cross = a*b + c*d;
    const double trace = col0*col0 + col1*col1;
    const double gap = std::hypot(col0*col0 - col1*col1, 2.0*cross);
    const double det = a*d - b*c;
    if (!(trace > 0.0) || !std::isfinite(trace) || !std::isfinite(gap) || !std::isfinite(det)) return 0.0;
    const double lambda_max = 0.5 * (trace + gap);
    if (!(lambda_max > 0.0) || !std::isfinite(lambda_max)) return 0.0;
    const double sigma_max = std::sqrt(lambda_max);
    const double sigma_min = std::abs(det) / sigma_max;
    if (!std::isfinite(sigma_min)) return 0.0;
    // Principal overlaps are cosines and therefore belong to [0, 1].  The
    // upper clamp only removes harmless round-off above one.
    return std::min(1.0, std::max(0.0, sigma_min));
}

} // namespace

double mass_weighted_modal_overlap(const CylindricalModalResult& a, std::size_t mode_a,
                                   const CylindricalModalResult& b, std::size_t mode_b,
                                   const detail::SparseMatrix& metric_mass) {
    if (mode_a >= a.mode_shapes.size() || mode_b >= b.mode_shapes.size())
        throw std::out_of_range("modal tracking mode index out of range");
    return normalized_cross(a.mode_shapes[mode_a], b.mode_shapes[mode_b], metric_mass);
}

double mass_weighted_subspace_similarity(
    const std::vector<std::vector<double>>& a,
    const std::vector<std::vector<double>>& b,
    const detail::SparseMatrix& metric_mass) {
    if (a.empty() || a.size() != b.size() || a.size() > 2)
        throw std::invalid_argument("current modal subspace tracker supports equal dimensions 1..2");
    for (const auto& v : a)
        if (v.size() != static_cast<std::size_t>(metric_mass.n))
            throw std::invalid_argument("modal tracking subspace vector dimension mismatch");
    for (const auto& v : b)
        if (v.size() != static_cast<std::size_t>(metric_mass.n))
            throw std::invalid_argument("modal tracking subspace vector dimension mismatch");

    if (a.size() == 1) return normalized_cross(a[0], b[0], metric_mass);

    const std::vector<double> a0 = metric_orthogonalized(a[0], {}, metric_mass);
    const std::vector<double> a1 = metric_orthogonalized(a[1], {a0}, metric_mass);
    const std::vector<double> b0 = metric_orthogonalized(b[0], {}, metric_mass);
    const std::vector<double> b1 = metric_orthogonalized(b[1], {b0}, metric_mass);
    if (a0.empty() || a1.empty() || b0.empty() || b1.empty()) return 0.0;
    const double c00 = dot_metric(a0, b0, metric_mass);
    const double c01 = dot_metric(a0, b1, metric_mass);
    const double c10 = dot_metric(a1, b0, metric_mass);
    const double c11 = dot_metric(a1, b1, metric_mass);
    return smallest_singular_value_2x2(c00, c01, c10, c11);
}

double mass_weighted_subspace_similarity(const CylindricalModalResult& a,
                                          const std::vector<std::size_t>& indices_a,
                                          const CylindricalModalResult& b,
                                          const std::vector<std::size_t>& indices_b,
                                          const detail::SparseMatrix& metric_mass) {
    if (indices_a.empty() || indices_a.size() != indices_b.size() || indices_a.size() > 2)
        throw std::invalid_argument("current modal subspace tracker supports equal dimensions 1..2");
    for (std::size_t index : indices_a)
        if (index >= a.mode_shapes.size()) throw std::out_of_range("modal tracking subspace index out of bounds");
    for (std::size_t index : indices_b)
        if (index >= b.mode_shapes.size()) throw std::out_of_range("modal tracking subspace index out of bounds");
    std::vector<std::vector<double>> va, vb;
    for (std::size_t index : indices_a) va.push_back(a.mode_shapes[index]);
    for (std::size_t index : indices_b) vb.push_back(b.mode_shapes[index]);
    return mass_weighted_subspace_similarity(va, vb, metric_mass);
}

double mass_weighted_subspace_similarity(const CylindricalModalResult& a,
                                          std::size_t first_a,
                                          const CylindricalModalResult& b,
                                          std::size_t first_b,
                                          std::size_t dimension,
                                          const detail::SparseMatrix& metric_mass) {
    if (dimension == 0 || dimension > 2)
        throw std::invalid_argument("current modal subspace tracker supports dimensions 1..2");
    std::vector<std::size_t> ia;
    std::vector<std::size_t> ib;
    ia.reserve(dimension); ib.reserve(dimension);
    for (std::size_t i = 0; i < dimension; ++i) {
        ia.push_back(first_a + i);
        ib.push_back(first_b + i);
    }
    return mass_weighted_subspace_similarity(a, ia, b, ib, metric_mass);
}

std::vector<std::vector<double>> prolongate_cylindrical_modes(
    const CylindricalModalResult& source,
    const std::vector<int>& source_harmonics,
    const std::vector<int>& target_harmonics) {
    auto validate = [](const std::vector<int>& h, const char* name) {
        if (h.empty() || h.front() != 0) throw std::invalid_argument(std::string(name) + " must start at zero");
        for (std::size_t i = 1; i < h.size(); ++i)
            if (h[i] <= h[i - 1] || h[i] < 0)
                throw std::invalid_argument(std::string(name) + " must be increasing and nonnegative");
    };
    validate(source_harmonics, "source harmonics");
    validate(target_harmonics, "target harmonics");

    const std::size_t source_basis = source_harmonics.size() +
        std::count_if(source_harmonics.begin(), source_harmonics.end(), [](int h) { return h != 0; });
    const std::size_t target_basis = target_harmonics.size() +
        std::count_if(target_harmonics.begin(), target_harmonics.end(), [](int h) { return h != 0; });
    if (source_basis == 0 || source_basis * 3 == 0) throw std::invalid_argument("invalid source basis");
    if (target_basis < source_basis) throw std::invalid_argument("target harmonic basis cannot be smaller than source basis");
    if (source.mode_shapes.empty()) return {};
    if (source.mode_shapes.front().size() % (3u * source_basis) != 0)
        throw std::invalid_argument("source modal vector is incompatible with source harmonics");

    const std::size_t nodes = source.mode_shapes.front().size() / (3u * source_basis);
    std::vector<std::size_t> source_to_target(source_basis, static_cast<std::size_t>(-1));
    std::size_t source_index = 0;
    for (int h : source_harmonics) {
        auto it = std::find(target_harmonics.begin(), target_harmonics.end(), h);
        if (it == target_harmonics.end()) throw std::invalid_argument("target harmonics do not contain source harmonic");
        const std::size_t target_index = static_cast<std::size_t>(std::distance(target_harmonics.begin(), it));
        // Each nonzero harmonic contributes cosine then sine basis entries.
        const std::size_t count = h == 0 ? 1u : 2u;
        std::size_t target_offset = 0;
        for (std::size_t k = 0; k < target_index; ++k) target_offset += target_harmonics[k] == 0 ? 1u : 2u;
        for (std::size_t b = 0; b < count; ++b) source_to_target[source_index++] = target_offset + b;
    }

    std::vector<std::vector<double>> out;
    out.reserve(source.mode_shapes.size());
    for (const auto& src : source.mode_shapes) {
        if (src.size() != nodes * 3u * source_basis)
            throw std::invalid_argument("source modal vectors have inconsistent dimensions");
        std::vector<double> dst(nodes * 3u * target_basis, 0.0);
        for (std::size_t node = 0; node < nodes; ++node) {
            for (std::size_t component = 0; component < 3u; ++component) {
                for (std::size_t sb = 0; sb < source_basis; ++sb) {
                    const std::size_t tb = source_to_target[sb];
                    dst[(node * 3u + component) * target_basis + tb] =
                        src[(node * 3u + component) * source_basis + sb];
                }
            }
        }
        out.push_back(std::move(dst));
    }
    return out;
}

} // namespace BellOptimizer::FEM
