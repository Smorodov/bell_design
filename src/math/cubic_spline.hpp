#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace BellMath {

// Natural cubic spline. API intentionally mirrors the small subset formerly
// used from GSL: build, value and derivative.
class CubicSpline {
public:
    bool build(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() < 3 || x.size() != y.size()) return false;
        for (size_t i = 1; i < x.size(); ++i) if (!(x[i] > x[i - 1])) return false;
        x_ = x; y_ = y;
        const size_t n = x.size();
        second_.assign(n, 0.0);
        std::vector<double> u(n - 1, 0.0);
        for (size_t i = 1; i + 1 < n; ++i) {
            const double h0 = x_[i] - x_[i - 1];
            const double h1 = x_[i + 1] - x_[i];
            const double sig = h0 / (h0 + h1);
            const double p = sig * second_[i - 1] + 2.0;
            second_[i] = (sig - 1.0) / p;
            const double d = (y_[i + 1] - y_[i]) / h1 - (y_[i] - y_[i - 1]) / h0;
            u[i] = (6.0 * d / (h0 + h1) - sig * u[i - 1]) / p;
        }
        for (size_t k = n - 2; k > 0; --k)
            second_[k] = second_[k] * second_[k + 1] + u[k];
        built_ = true;
        return true;
    }

    bool valid() const noexcept { return built_; }
    void clear() { x_.clear(); y_.clear(); second_.clear(); built_ = false; }

    double value(double x) const { return evaluate(x, false); }
    double derivative(double x) const { return evaluate(x, true); }

private:
    double evaluate(double x, bool deriv) const {
        if (!built_) return 0.0;
        if (x <= x_.front()) x = x_.front();
        if (x >= x_.back()) x = x_.back();
        auto it = std::upper_bound(x_.begin(), x_.end(), x);
        size_t k = (it == x_.begin()) ? 0 : static_cast<size_t>(it - x_.begin() - 1);
        if (k + 1 >= x_.size()) k = x_.size() - 2;
        const double h = x_[k + 1] - x_[k];
        const double a = (x_[k + 1] - x) / h;
        const double b = (x - x_[k]) / h;
        if (!deriv) {
            return a * y_[k] + b * y_[k + 1] +
                   ((a*a*a - a) * second_[k] + (b*b*b - b) * second_[k + 1]) * h*h / 6.0;
        }
        return (y_[k + 1] - y_[k]) / h +
               ((-3.0*a*a + 1.0) * second_[k] + (3.0*b*b - 1.0) * second_[k + 1]) * h / 6.0;
    }
    std::vector<double> x_, y_, second_;
    bool built_ = false;
};

template<class F>
inline double adaptive_simpson_impl(const F& f, double a, double b,
                                    double fa, double fb, double fm,
                                    double whole, double tol, int depth) {
    const double m = 0.5 * (a + b);
    const double lm = 0.5 * (a + m);
    const double rm = 0.5 * (m + b);
    const double flm = f(lm), frm = f(rm);
    const double left = (m - a) * (fa + 4.0 * flm + fm) / 6.0;
    const double right = (b - m) * (fm + 4.0 * frm + fb) / 6.0;
    const double delta = left + right - whole;
    if (depth <= 0 || std::abs(delta) <= 15.0 * tol)
        return left + right + delta / 15.0;
    return adaptive_simpson_impl(f, a, m, fa, fm, flm, left, tol * 0.5, depth - 1) +
           adaptive_simpson_impl(f, m, b, fm, fb, frm, right, tol * 0.5, depth - 1);
}

template<class F>
double adaptive_simpson(F f, double a, double b, double tolerance = 1e-8, int max_depth = 18) {
    if (a == b) return 0.0;
    if (b < a) return -adaptive_simpson(f, b, a, tolerance, max_depth);
    const double m = 0.5 * (a + b);
    const double fa = f(a), fb = f(b), fm = f(m);
    const double whole = (b - a) * (fa + 4.0 * fm + fb) / 6.0;
    return adaptive_simpson_impl(f, a, b, fa, fb, fm, whole, tolerance, max_depth);
}

} // namespace BellMath
