// bell_generator.cpp
#include "bell_generator.h"
#include "math_constants.hpp"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstring>
#include <limits>

namespace {

constexpr double PI = 3.1415926535897932384626433832795;

struct SegmentPoint { double r; double z; };

int orientation(const SegmentPoint& a, const SegmentPoint& b, const SegmentPoint& c) {
    const double v = (b.r-a.r)*(c.z-a.z) - (b.z-a.z)*(c.r-a.r);
    const double scale = std::max({1.0, std::abs(b.r-a.r), std::abs(b.z-a.z),
                                   std::abs(c.r-a.r), std::abs(c.z-a.z)});
    const double eps = 1e-12 * scale * scale;
    if (v > eps) return 1;
    if (v < -eps) return -1;
    return 0;
}

bool on_segment(const SegmentPoint& a, const SegmentPoint& b, const SegmentPoint& p) {
    const double eps = 1e-10;
    return p.r >= std::min(a.r,b.r)-eps && p.r <= std::max(a.r,b.r)+eps &&
           p.z >= std::min(a.z,b.z)-eps && p.z <= std::max(a.z,b.z)+eps;
}

bool segments_intersect(const SegmentPoint& a, const SegmentPoint& b,
                        const SegmentPoint& c, const SegmentPoint& d) {
    const int o1=orientation(a,b,c), o2=orientation(a,b,d);
    const int o3=orientation(c,d,a), o4=orientation(c,d,b);
    if (o1 != o2 && o3 != o4) return true;
    return (o1==0 && on_segment(a,b,c)) || (o2==0 && on_segment(a,b,d)) ||
           (o3==0 && on_segment(c,d,a)) || (o4==0 && on_segment(c,d,b));
}

} // namespace


// ============================================================
// КОНСТРУКТОР
// ============================================================
BellGenerator::BellGenerator(const BellMaterial& material, const BellMeshConfig& config)
    : material_(material), config_(config) {}

BellGenerator::~BellGenerator() = default;

// ============================================================
// УПРАВЛЕНИЕ ГЕОМЕТРИЕЙ
// ============================================================
void BellGenerator::setGeometry(const BellGeometry& geom) {
    geometry_ = geom;
    spline_.clear();
    is_mesh_built_ = false;
    mesh_data_.clear();
}

// ============================================================
// ПОСТРОЕНИЕ СПЛАЙНА
// ============================================================
bool BellGenerator::buildSpline() {
    const size_t n = geometry_.r_knots.size();
    if (n < 3) { logError("At least 3 knots required"); return false; }
    std::vector<double> z_knots(n), r_knots(n);
    for (size_t i = 0; i < n; ++i) {
        z_knots[i] = geometry_.getKnotZ(i);
        r_knots[i] = geometry_.r_knots[i] * geometry_.R_mouth;
    }
    if (std::abs(z_knots.front()) > 1e-6 || std::abs(r_knots.front()) > 1e-6) {
        logError("First knot must be (0, 0)"); return false;
    }
    if (std::abs(z_knots.back() - geometry_.H) > 1e-6 || std::abs(r_knots.back() - geometry_.R_mouth) > 1e-6) {
        logError("Last knot must be (H, R_mouth)"); return false;
    }
    for (size_t i = 1; i < n; ++i) {
        if (z_knots[i] <= z_knots[i-1] || r_knots[i] <= r_knots[i-1]) {
            logError("Spline knots must be strictly increasing"); return false;
        }
    }
    if (!spline_.build(z_knots, r_knots)) { logError("Failed to build cubic spline"); return false; }
    return true;
}

// ============================================================
// ИНТЕГРАНТ ДЛЯ ДЛИНЫ ДУГИ
// ============================================================
double BellGenerator::integrandParametric(double t) {
    if (!spline_.valid()) return 0.0;
    const double dr_dz = derivative(t * geometry_.H);
    return std::hypot(geometry_.H, dr_dz * geometry_.H);
}

double BellGenerator::computeArcLengthToT(double t_target) const {
    if (!spline_.valid() || t_target <= 0.0) return 0.0;
    const double t = std::clamp(t_target, 0.0, 1.0);
    auto f = [this](double u) { return const_cast<BellGenerator*>(this)->integrandParametric(u); };
    return BellMath::adaptive_simpson(f, 0.0, t, 1e-8, 18);
}

double BellGenerator::computeArcLength() {
    if (!spline_.valid()) {
        logError("Spline not initialized");
        return 0.0;
    }
    
    if (geometry_.H <= 0.0) {
        logError("H must be > 0");
        return 0.0;
    }
    
    double total = computeArcLengthToT(1.0);
    if (total < 1e-6) {
        logError("Arc length calculation failed - got zero");
        return 0.0;
    }
    
    mesh_data_.total_arc_length = total;
    return total;
}

// ============================================================
// ПОИСК Z ПО ДЛИНЕ ДУГИ
// ============================================================
double BellGenerator::getZfromS(double s) const {
    if (s <= 0.0) return 0.0;
    if (s >= mesh_data_.total_arc_length) return geometry_.H;
    
    double t_low = 0.0;
    double t_high = 1.0;
    
    for (int i = 0; i < 100; ++i) {
        double t_mid = (t_low + t_high) / 2.0;
        double s_mid = computeArcLengthToT(t_mid);
        
        if (s_mid < s) {
            t_low = t_mid;
        } else {
            t_high = t_mid;
        }
        
        if (t_high - t_low < 1e-10) break;
    }
    
    double t = (t_low + t_high) / 2.0;
    return t * geometry_.H;
}

// ============================================================
// РАДИУС, ПРОИЗВОДНАЯ, ТОЛЩИНА
// ============================================================
double BellGenerator::radius(double z) const {
    if (!spline_.valid()) return 0.0;  // ДОБАВИЛИ ПРОВЕРКУ
    if (z <= 0.0) return 0.0;
    if (z >= geometry_.H) return geometry_.R_mouth;
    return spline_.value(z);
}

double BellGenerator::derivative(double z) const {
    if (!spline_.valid()) return 0.0;  // ДОБАВИЛИ ПРОВЕРКУ
    if (z <= 0.0) return 0.0;
    if (z >= geometry_.H) {
        return spline_.derivative(geometry_.H);
    }
    return spline_.derivative(z);
}

double BellGenerator::getThicknessAt(double z) const {
    if (z <= 0.0) return geometry_.t_top;
    if (z >= geometry_.H) return geometry_.t_bottom;
    double t = z / geometry_.H;
    return geometry_.t_top + (geometry_.t_bottom - geometry_.t_top) * t;
}

// ============================================================
// ВНЕШНИЙ И ВНУТРЕННИЙ РАДИУСЫ
// ============================================================
double BellGenerator::getOuterRadius(double z) const {
    if (z < 0.0 || z > geometry_.H) return 0.0;
    return radius(z);
}

double BellGenerator::getInnerRadius(double z) const {
    if (z < 0.0 || z > geometry_.H) return 0.0;
    if (z <= 0.0) return 0.0;
    double r = radius(z);
    const double dr_dz = derivative(z);
    const double normal_length = std::hypot(1.0, dr_dz);
    const double thickness = getThicknessAt(z);
    // The inward unit normal to the meridional curve r(z) is
    // n = (-1, dr/dz) / sqrt(1 + (dr/dz)^2).
    return std::max(0.0, r - thickness / normal_length);
}

// ============================================================
// ПОСТРОЕНИЕ ВНЕШНЕГО ПРОФИЛЯ
// ============================================================
void BellGenerator::computeOuterProfile() {
    const int n = config_.nNodes;
    mesh_data_.outer_profile.clear();
    mesh_data_.outer_profile.reserve(n);
    inner_axis_z_ = -1.0;

    // Find where the true inward normal offset first intersects the axis.
    // Clamping several negative inner radii to zero creates stacked nodes
    // and degenerate Q4 elements, so the mesh is redistributed from just
    // beyond this intersection instead.
    const auto raw_inner_radius = [this](double z) {
        const double slope = derivative(z);
        return radius(z) - getThicknessAt(z) / std::hypot(1.0, slope);
    };
    double previous_z = 0.0;
    double previous_value = raw_inner_radius(previous_z);
    for (int sample = 1; sample <= 2048; ++sample) {
        const double z = geometry_.H * static_cast<double>(sample) / 2048.0;
        const double value = raw_inner_radius(z);
        if (previous_value <= 0.0 && value > 0.0) {
            double lo = previous_z;
            double hi = z;
            for (int iteration = 0; iteration < 60; ++iteration) {
                const double mid = 0.5 * (lo + hi);
                if (raw_inner_radius(mid) > 0.0) hi = mid;
                else lo = mid;
            }
            inner_axis_z_ = 0.5 * (lo + hi);
            break;
        }
        previous_z = z;
        previous_value = value;
    }

    if (n < 3 || inner_axis_z_ < 0.0) {
        logError("Could not locate inner-profile axis intersection");
        return;
    }

    const double s_axis = computeArcLengthToT(inner_axis_z_ / geometry_.H);
    const double nominal_step = (mesh_data_.total_arc_length - s_axis) / (n - 2);
    // Keep the first non-axis inner node safely beyond the offset intersection.
    const double s_first = std::min(mesh_data_.total_arc_length,
                                    s_axis + 0.25 * nominal_step);
    mesh_data_.outer_profile.push_back({0.0, 0.0});
    for (int i = 1; i < n; ++i) {
        const double s_target = s_first + static_cast<double>(i - 1) *
            (mesh_data_.total_arc_length - s_first) / (n - 2);
        const double z = getZfromS(s_target);
        mesh_data_.outer_profile.push_back({radius(z), z});
    }
}

// ============================================================
// ВЫЧИСЛЕНИЕ ЕДИНИЧНОЙ ВНУТРЕННЕЙ НОРМАЛИ
// ============================================================
void BellGenerator::computeNormals(std::vector<double>& n_r, std::vector<double>& n_z) const {
    int n = config_.nNodes;
    n_r.resize(n);
    n_z.resize(n);
    
    for (int i = 0; i < n; ++i) {
        double z = mesh_data_.outer_profile[i].z;
        double r = mesh_data_.outer_profile[i].r;
        
        if (i == 0 || (r < 1e-12 && z < 1e-12)) {
            n_r[i] = 0.0;
            n_z[i] = 1.0;
            continue;
        }
        
        if (r < 1e-12) {
            n_r[i] = 0.0;
            n_z[i] = 1.0;
            continue;
        }
        
        const double dr_dz = derivative(z);
        const double len = std::hypot(1.0, dr_dz);

        // Unit inward normal for the curve (r(z), z): perpendicular to
        // tangent (dr/dz, 1), with a negative radial component.
        n_r[i] = -1.0 / len;
        n_z[i] = dr_dz / len;
        
        double final_len = std::sqrt(n_r[i] * n_r[i] + n_z[i] * n_z[i]);
        if (final_len > 1e-12) {
            n_r[i] /= final_len;
            n_z[i] /= final_len;
        }
    }
}

// ============================================================
// ПОСТРОЕНИЕ ВНУТРЕННЕГО ПРОФИЛЯ
// ============================================================
void BellGenerator::computeInnerProfile(const std::vector<double>& n_r, const std::vector<double>& n_z) {
    int n = config_.nNodes;
    mesh_data_.inner_profile.clear();
    mesh_data_.inner_profile.reserve(n);
    
    // ДОБАВИЛИ ПРОВЕРКУ
    if (n_r.size() < (size_t)n || n_z.size() < (size_t)n) {
        logError("Normal vectors size mismatch");
        return;
    }
    
    for (int i = 0; i < n; ++i) {
        double z_outer = mesh_data_.outer_profile[i].z;
        double r_outer = mesh_data_.outer_profile[i].r;
        double t = getThicknessAt(z_outer);
        
        double r_inner, z_inner;
        if (i == 0) {
            if (inner_axis_z_ < 0.0) {
                logError("Inner-profile axis intersection is unavailable");
                return;
            }
            const double slope = derivative(inner_axis_z_);
            z_inner = inner_axis_z_ + getThicknessAt(inner_axis_z_) * slope /
                      std::hypot(1.0, slope);
            r_inner = 0.0;
        } else {
            r_inner = r_outer + t * n_r[i];
            z_inner = z_outer + t * n_z[i];
            if (!(r_inner > 0.0) || !std::isfinite(r_inner) || !std::isfinite(z_inner)) {
                logError("Inner offset crosses the axis again or becomes non-finite");
                return;
            }
        }

        mesh_data_.inner_profile.push_back({r_inner, z_inner});
    }
}

// ============================================================
// СБОРКА СЕТКИ
// ============================================================
void BellGenerator::assembleMesh() {
    int n = config_.nNodes;
    mesh_data_.nodes.clear();
    mesh_data_.nodes.reserve(2 * n);
    mesh_data_.elements.clear();
    mesh_data_.elements.reserve(n - 1);
    
    for (int i = 0; i < n; ++i) {
        Node node;
        node.id = i + 1;
        node.r = mesh_data_.outer_profile[i].r;
        node.z = mesh_data_.outer_profile[i].z;
        mesh_data_.nodes.push_back(node);
    }
    
    for (int i = 0; i < n; ++i) {
        Node node;
        node.id = i + 1 + n;
        node.r = mesh_data_.inner_profile[i].r;
        node.z = mesh_data_.inner_profile[i].z;
        mesh_data_.nodes.push_back(node);
    }
    
    for (int i = 0; i < n - 1; ++i) {
        Element elem;
        elem.n1 = i + 1 + n;
        elem.n2 = i + 1;
        elem.n3 = i + 2;
        elem.n4 = i + 2 + n;
        mesh_data_.elements.push_back(elem);
    }
}

// ============================================================
// ВЫЧИСЛЕНИЕ МЕТРИК
// ============================================================
void BellGenerator::computeMetrics() {
    double min_t = std::numeric_limits<double>::max();
    double max_t = 0.0;
    
    for (const auto& p : mesh_data_.outer_profile) {
        double t = getThicknessAt(p.z);
        if (t < min_t) min_t = t;
        if (t > max_t) max_t = t;
    }
    
    mesh_data_.min_thickness = min_t;
    mesh_data_.max_thickness = max_t;
}

// ============================================================
// ГЕОМЕТРИЧЕСКИЕ ПРОВЕРКИ
// ============================================================
bool BellGenerator::checkMonotonicity() const {
    const auto& r = geometry_.r_knots;
    for (size_t i = 1; i < r.size(); ++i) {
        if (r[i] <= r[i-1]) return false;
    }
    return true;
}

bool BellGenerator::checkSelfIntersection() const {
    const auto& inner = mesh_data_.inner_profile;
    if (inner.size() < 2) return true;

    for (size_t i = 1; i < inner.size(); ++i) {
        const SegmentPoint a{inner[i-1].r, inner[i-1].z};
        const SegmentPoint b{inner[i].r, inner[i].z};
        if (std::hypot(b.r-a.r, b.z-a.z) < 1e-9) return false;
    }

    // Adjacent segments share an endpoint by construction. Only non-adjacent
    // segments are candidates for a genuine profile loop.
    for (size_t i = 0; i + 1 < inner.size(); ++i) {
        const SegmentPoint a{inner[i].r, inner[i].z};
        const SegmentPoint b{inner[i+1].r, inner[i+1].z};
        for (size_t j = i + 2; j + 1 < inner.size(); ++j) {
            // The first/last pair is also non-adjacent and must be checked:
            // this is an open meridional profile, so no endpoint is shared.
            const SegmentPoint c{inner[j].r, inner[j].z};
            const SegmentPoint d{inner[j+1].r, inner[j+1].z};
            if (segments_intersect(a,b,c,d)) return false;
        }
    }
    return true;
}

bool BellGenerator::checkSlopeAngles() const {
    if (!spline_.valid() || geometry_.H <= 0.0) return false;

    const double max_slope = BellConstants::getMaxSlope();
    // A spline can overshoot between control knots, so sample the actual
    // derivative rather than using secant slopes between knots.
    constexpr int samples_per_interval = 16;
    const int intervals = static_cast<int>(geometry_.r_knots.size()) - 1;
    for (int interval = 0; interval < intervals; ++interval) {
        for (int k = 0; k <= samples_per_interval; ++k) {
            const double u = (static_cast<double>(k) / samples_per_interval);
            const double z = (interval + u) * geometry_.H / intervals;
            if (std::abs(derivative(z)) > max_slope + 1e-10) return false;
        }
    }
    return true;
}

// ============================================================
// ПРОВЕРКА КАЧЕСТВА ЭЛЕМЕНТОВ
// ============================================================
bool BellGenerator::checkElementQuality() const {
    const double min_angle_deg = 10.0;
    const double max_angle_deg = 175.0;
    const double max_aspect_ratio = 100.0;
    const double min_area = 1e-12;
    
    for (const auto& elem : mesh_data_.elements) {
        const Node& n1 = mesh_data_.nodes[elem.n1 - 1];
        const Node& n2 = mesh_data_.nodes[elem.n2 - 1];
        const Node& n3 = mesh_data_.nodes[elem.n3 - 1];
        const Node& n4 = mesh_data_.nodes[elem.n4 - 1];

        // The Q4 FEM mapping must preserve orientation at every Gauss point.
        // Polygon area alone is insufficient because a bilinear element can
        // have a valid signed area while its Jacobian changes sign internally.
        const double g = 1.0 / std::sqrt(3.0);
        for (double xi : {-g, g}) {
            for (double eta : {-g, g}) {
                const double dNdxi[4] = {-0.25*(1.0-eta), 0.25*(1.0-eta),
                                          0.25*(1.0+eta), -0.25*(1.0+eta)};
                const double dNdeta[4] = {-0.25*(1.0-xi), -0.25*(1.0+xi),
                                           0.25*(1.0+xi), 0.25*(1.0-xi)};
                const Node nodes[4] = {n1,n2,n3,n4};
                double j00=0.0,j01=0.0,j10=0.0,j11=0.0;
                for (int a=0; a<4; ++a) {
                    j00 += dNdxi[a] * nodes[a].r;
                    j01 += dNdxi[a] * nodes[a].z;
                    j10 += dNdeta[a] * nodes[a].r;
                    j11 += dNdeta[a] * nodes[a].z;
                }
                const double det_j = j00*j11 - j01*j10;
                if (!(det_j > 1e-12) || !std::isfinite(det_j)) return false;
            }
        }
        
        if (std::abs(n1.r - n2.r) < 1e-12 && std::abs(n1.z - n2.z) < 1e-12) return false;
        if (std::abs(n2.r - n3.r) < 1e-12 && std::abs(n2.z - n3.z) < 1e-12) return false;
        if (std::abs(n3.r - n4.r) < 1e-12 && std::abs(n3.z - n4.z) < 1e-12) return false;
        if (std::abs(n4.r - n1.r) < 1e-12 && std::abs(n4.z - n1.z) < 1e-12) return false;
        
        double v1x = n2.r - n1.r, v1z = n2.z - n1.z;
        double v2x = n3.r - n2.r, v2z = n3.z - n2.z;
        double v3x = n4.r - n3.r, v3z = n4.z - n3.z;
        double v4x = n1.r - n4.r, v4z = n1.z - n4.z;
        
        double l1 = std::sqrt(v1x*v1x + v1z*v1z);
        double l2 = std::sqrt(v2x*v2x + v2z*v2z);
        double l3 = std::sqrt(v3x*v3x + v3z*v3z);
        double l4 = std::sqrt(v4x*v4x + v4z*v4z);
        
        if (l1 < min_area || l2 < min_area || l3 < min_area || l4 < min_area) {
            return false;
        }
        
        double min_len = std::min({l1, l2, l3, l4});
        double max_len = std::max({l1, l2, l3, l4});
        if (max_len / min_len > max_aspect_ratio) {
            return false;
        }
        
        auto calc_angle = [](double ax, double az, double bx, double bz) -> double {
            double len_a = std::sqrt(ax*ax + az*az);
            double len_b = std::sqrt(bx*bx + bz*bz);
            if (len_a < 1e-12 || len_b < 1e-12) return 0.0;
            double cos_angle = (ax*bx + az*bz) / (len_a * len_b);
            cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
            return std::acos(cos_angle) * 180.0 / bell::kPi;
        };
        
        double a1 = calc_angle(-v1x, -v1z, v4x, v4z);
        double a2 = calc_angle(-v2x, -v2z, v1x, v1z);
        double a3 = calc_angle(-v3x, -v3z, v2x, v2z);
        double a4 = calc_angle(-v4x, -v4z, v3x, v3z);
        
        if (a1 < min_angle_deg || a2 < min_angle_deg || a3 < min_angle_deg || a4 < min_angle_deg) {
            return false;
        }
        
        if (a1 > max_angle_deg || a2 > max_angle_deg || a3 > max_angle_deg || a4 > max_angle_deg) {
            return false;
        }
        
        double area = 0.5 * ((n1.r*n2.z + n2.r*n3.z + n3.r*n4.z + n4.r*n1.z) -
                             (n2.r*n1.z + n3.r*n2.z + n4.r*n3.z + n1.r*n4.z));
        if (area < min_area) {
            return false;
        }
    }
    return true;
}

// ============================================================
// ПРЕДВАРИТЕЛЬНАЯ ТЕХНОЛОГИЧЕСКАЯ ПРОВЕРКА
// ============================================================
bool BellGenerator::checkBasicCastability(GeometryValidationResult* result) const {
    if (!result || mesh_data_.outer_profile.empty()) return true;

    const double min_t = mesh_data_.min_thickness;
    const double max_t = mesh_data_.max_thickness;
    if (!(min_t > 0.0) || !(max_t >= min_t) || !std::isfinite(min_t) || !std::isfinite(max_t)) {
        result->passed = false;
        result->reason = "Invalid wall thickness for preliminary castability check";
        return false;
    }

    result->min_thickness = min_t;
    result->max_thickness = max_t;
    result->thickness_ratio = max_t / min_t;

    if (result->thickness_ratio > BellConstants::MAX_CAST_THICKNESS_RATIO) {
        result->passed = false;
        result->reason = "Excessive wall-thickness variation for preliminary castability";
        return false;
    }

    // The current wall law is linear in z. Measure the actual spatial gradient
    // along the sampled outer profile rather than relying only on endpoint data.
    double max_gradient = 0.0;
    for (size_t i = 1; i < mesh_data_.outer_profile.size(); ++i) {
        const auto& a = mesh_data_.outer_profile[i - 1];
        const auto& b = mesh_data_.outer_profile[i];
        const double ds = std::hypot(b.r - a.r, b.z - a.z);
        if (ds <= 1e-12) continue;
        const double dt = getThicknessAt(b.z) - getThicknessAt(a.z);
        max_gradient = std::max(max_gradient, std::abs(dt) / ds);
    }
    result->max_thickness_gradient = max_gradient;
    if (max_gradient > BellConstants::MAX_CAST_THICKNESS_GRADIENT) {
        result->passed = false;
        result->reason = "Excessive wall-thickness gradient for preliminary castability";
        return false;
    }

    return true;
}

void BellGenerator::computeProfileQuality(GeometryValidationResult* result) const {
    if (!result || !spline_.valid() || geometry_.H <= 0.0) return;

    constexpr int samples = 128;
    double previous_slope = derivative(0.0);
    double max_slope_change = 0.0;
    double max_angle = 0.0;
    for (int i = 1; i <= samples; ++i) {
        const double z = geometry_.H * static_cast<double>(i) / samples;
        const double slope = derivative(z);
        const double angle = std::atan(std::abs(slope)) * 180.0 / PI;
        max_angle = std::max(max_angle, angle);
        max_slope_change = std::max(max_slope_change, std::abs(slope - previous_slope));
        previous_slope = slope;
    }
    result->max_slope_angle = max_angle;
    result->max_slope_change = max_slope_change;
}

// ============================================================
// ПОЛНАЯ ГЕОМЕТРИЧЕСКАЯ ПРОВЕРКА
// ============================================================
GeometryValidationResult BellGenerator::validateGeometry() const {
    GeometryValidationResult result;
    result.passed = true;

    if (!geometry_.isValid()) {
        result.passed = false;
        result.reason = "Invalid geometry parameters";
        return result;
    }
    
    if (!checkMonotonicity()) {
        result.passed = false;
        result.reason = "Non-monotonic radius knots";
        double penalty = 0.0;
        const auto& r = geometry_.r_knots;
        for (size_t i = 1; i < r.size(); ++i) {
            if (r[i] <= r[i-1]) {
                penalty += (r[i-1] - r[i]) * (r[i-1] - r[i]);
            }
        }
        result.penalty_value = penalty / (r.size() - 1);
        return result;
    }
    
    if (!mesh_data_.inner_profile.empty()) {
        if (!checkSelfIntersection()) {
            result.passed = false;
            result.reason = "Inner profile self-intersection (loops detected)";
            return result;
        }
    }
    
    return result;
}

// ============================================================
// ПОСТРОЕНИЕ СЕТКИ
// ============================================================
bool BellGenerator::buildMesh() {
    GeometryValidationResult dummy;
    return buildMesh(&dummy);
}

bool BellGenerator::buildMesh(GeometryValidationResult* validation_out) {
    mesh_data_.clear();
    is_mesh_built_ = false;
    
    auto validation = validateGeometry();
    if (validation_out) *validation_out = validation;
    if (!validation.passed) {
        return false;
    }
    
    if (!buildSpline()) {
        logError("Failed to build spline");
        return false;
    }
    
    double arc_len = computeArcLength();
    if (arc_len < 1e-6) {
        logError("Arc length is zero");
        return false;
    }
    
    computeOuterProfile();
    if (mesh_data_.outer_profile.size() != static_cast<size_t>(config_.nNodes)) {
        if (validation_out) {
            validation_out->passed = false;
            validation_out->reason = "Could not construct valid inner-axis intersection";
        }
        return false;
    }

    std::vector<double> n_r, n_z;
    computeNormals(n_r, n_z);
    computeInnerProfile(n_r, n_z);
    if (mesh_data_.inner_profile.size() != static_cast<size_t>(config_.nNodes)) {
        if (validation_out) {
            validation_out->passed = false;
            validation_out->reason = "Could not construct a positive-radius inner offset";
        }
        return false;
    }

    if (!checkSelfIntersection()) {
        logError("Self-intersection detected in inner profile");
        if (validation_out) {
            validation_out->passed = false;
            validation_out->reason = "Inner profile self-intersection (loops)";
        }
        return false;
    }
    
    assembleMesh();
    computeMetrics();

    GeometryValidationResult castability;
    castability.passed = true;
    if (!checkBasicCastability(&castability)) {
        if (validation_out) {
            *validation_out = castability;
        }
        return false;
    }

    if (validation_out) {
        validation_out->min_thickness = castability.min_thickness;
        validation_out->max_thickness = castability.max_thickness;
        validation_out->thickness_ratio = castability.thickness_ratio;
        validation_out->max_thickness_gradient = castability.max_thickness_gradient;
    }
    computeProfileQuality(validation_out);

    if (!checkElementQuality()) {
        logError("Element quality check failed: invalid Jacobian or quality limit exceeded");
        if (validation_out) {
            validation_out->passed = false;
            validation_out->reason = "Element quality limits exceeded";
        }
        return false;
    }
    
    is_mesh_built_ = true;
    return true;
}

// ============================================================
// СБРОС
// ============================================================
void BellGenerator::reset() {
    mesh_data_.clear();
    spline_.clear();
    inner_axis_z_ = -1.0;
    is_mesh_built_ = false;
}

// ============================================================
// ЛОГИРОВАНИЕ
// ============================================================
void BellGenerator::logWarning(const std::string&) const {
    // Подавлен — не выводим предупреждения
}

void BellGenerator::logError(const std::string&) const {
    // Подавлен — не выводим ошибки
}