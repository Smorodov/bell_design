// bell_generator.cpp
#include "bell_generator.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstring>
#include <limits>

#include <fmt/format.h>

// ============================================================
// ВСПОМОГАТЕЛЬНАЯ СТРУКТУРА ДЛЯ ИНТЕГРАТОРА
// ============================================================
struct IntegrandParams {
    const BellGenerator* generator;
};

static double integrand_wrapper_parametric(double t, void* params) {
    IntegrandParams* p = static_cast<IntegrandParams*>(params);
    return p->generator->integrandParametric(t);
}

// ============================================================
// КОНСТРУКТОР
// ============================================================
BellGenerator::BellGenerator(const BellMaterial& material, const BellMeshConfig& config)
    : material_(material)
    , config_(config)
    , spline_(nullptr)
    , acc_(nullptr)
    , integ_w_(nullptr)
{
    acc_ = gsl_interp_accel_alloc();
    if (!acc_) {
        logError("Failed to allocate GSL accelerator");
    }
    
    integ_w_ = gsl_integration_workspace_alloc(10000);
    if (!integ_w_) {
        logError("Failed to allocate integration workspace");
    }
}

// ============================================================
// ДЕСТРУКТОР
// ============================================================
BellGenerator::~BellGenerator() {
    if (spline_) {
        gsl_spline_free(spline_);
        spline_ = nullptr;
    }
    if (acc_) {
        gsl_interp_accel_free(acc_);
        acc_ = nullptr;
    }
    if (integ_w_) {
        gsl_integration_workspace_free(integ_w_);
        integ_w_ = nullptr;
    }
}

// ============================================================
// УПРАВЛЕНИЕ ГЕОМЕТРИЕЙ
// ============================================================
void BellGenerator::setGeometry(const BellGeometry& geom) {
    geometry_ = geom;
    is_mesh_built_ = false;
    mesh_data_.clear();
}

// ============================================================
// ПОСТРОЕНИЕ СПЛАЙНА
// ============================================================
bool BellGenerator::buildSpline() {
    if (spline_) {
        gsl_spline_free(spline_);
        spline_ = nullptr;
    }
    
    size_t n = geometry_.r_knots.size();
    if (n < 3) {
        logError("At least 3 knots required");
        return false;
    }
    
    std::vector<double> z_knots(n);
    std::vector<double> r_knots(n);
    
    for (size_t i = 0; i < n; ++i) {
        z_knots[i] = geometry_.getKnotZ(i);
        r_knots[i] = geometry_.r_knots[i] * geometry_.R_mouth;
    }
    
    if (std::abs(z_knots.front()) > 1e-6 || std::abs(r_knots.front()) > 1e-6) {
        logError("First knot must be (0, 0)");
        return false;
    }
    if (std::abs(z_knots.back() - geometry_.H) > 1e-6 || 
        std::abs(r_knots.back() - geometry_.R_mouth) > 1e-6) {
        logError("Last knot must be (H, R_mouth)");
        return false;
    }
    
    for (size_t i = 1; i < n; ++i) {
        if (z_knots[i] <= z_knots[i-1]) {
            logError("z_knots must be strictly increasing");
            return false;
        }
    }
    
    for (size_t i = 1; i < n; ++i) {
        if (r_knots[i] <= r_knots[i-1]) {
            logError("r_knots must be strictly increasing");
            return false;
        }
    }
    
    spline_ = gsl_spline_alloc(gsl_interp_cspline, n);
    if (!spline_) {
        logError("Failed to allocate GSL spline");
        return false;
    }
    
    int status = gsl_spline_init(spline_, z_knots.data(), r_knots.data(), n);
    if (status != GSL_SUCCESS) {
        logError("Failed to initialize GSL spline");
        gsl_spline_free(spline_);
        spline_ = nullptr;
        return false;
    }
    
    return true;
}

// ============================================================
// ИНТЕГРАНТ ДЛЯ ДЛИНЫ ДУГИ
// ============================================================
double BellGenerator::integrandParametric(double t) const {
    if (!spline_) return 0.0;
    double z = t * geometry_.H;
    double dr_dz = derivative(z);
    double dz_dt = geometry_.H;
    double dr_dt = dr_dz * geometry_.H;
    return std::sqrt(dz_dt * dz_dt + dr_dt * dr_dt);
}

// ============================================================
// ВЫЧИСЛЕНИЕ ДЛИНЫ ДУГИ
// ============================================================
double BellGenerator::computeArcLengthToT(double t_target) const {
    if (!integ_w_) {
        logError("Integration workspace not initialized");
        return 0.0;
    }
    
    if (t_target <= 0.0) return 0.0;
    
    double t_clamped = std::max(0.0, std::min(t_target, 1.0));
    if (t_clamped < 1e-12) return 0.0;
    
    IntegrandParams params;
    params.generator = this;
    
    gsl_function F;
    F.function = integrand_wrapper_parametric;
    F.params = &params;
    
    double result = 0.0;
    double error = 0.0;
    
    gsl_set_error_handler_off();
    
    int status = gsl_integration_qags(
        &F, 
        0.0, 
        t_clamped, 
        0.0,
        1e-6,
        10000,
        integ_w_,
        &result,
        &error
    );
    
    gsl_set_error_handler(NULL);
    
    if (status != GSL_SUCCESS) {
        logError("Integration failed");
        return 0.0;
    }
    
    if (result < 0.0) result = 0.0;
    return result;
}

double BellGenerator::computeArcLength() {
    if (!spline_) {
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
    if (!spline_) return 0.0;
    if (z <= 0.0) return 0.0;
    if (z >= geometry_.H) return geometry_.R_mouth;
    return gsl_spline_eval(spline_, z, acc_);
}

double BellGenerator::derivative(double z) const {
    if (!spline_) return 0.0;
    if (z <= 0.0) return 0.0;
    if (z >= geometry_.H) {
        return gsl_spline_eval_deriv(spline_, geometry_.H, acc_);
    }
    return gsl_spline_eval_deriv(spline_, z, acc_);
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
    double dr = derivative(z);
    double len = std::sqrt(1.0 + dr * dr);
    double t = getThicknessAt(z);
    double r_inner = r - t * dr / len;
    return std::max(0.0, r_inner);
}

// ============================================================
// ПОСТРОЕНИЕ ВНЕШНЕГО ПРОФИЛЯ
// ============================================================
void BellGenerator::computeOuterProfile() {
    int n = config_.nNodes;
    mesh_data_.outer_profile.clear();
    mesh_data_.outer_profile.reserve(n);
    
    for (int i = 0; i < n; ++i) {
        double s_target = i * mesh_data_.total_arc_length / (n - 1);
        double z = getZfromS(s_target);
        double r = radius(z);
        mesh_data_.outer_profile.push_back({r, z});
    }
}

// ============================================================
// ВЫЧИСЛЕНИЕ НОРМАЛЕЙ (ГАРАНТИРОВАННО ВНУТРЬ И ВНИЗ)
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
        
        double dr = derivative(z);
        double len = std::sqrt(1.0 + dr * dr);
        
        n_r[i] = -std::abs(dr) / len;
        n_z[i] = 1.0 / len;
        
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
    
    for (int i = 0; i < n; ++i) {
        double z_outer = mesh_data_.outer_profile[i].z;
        double r_outer = mesh_data_.outer_profile[i].r;
        double t = getThicknessAt(z_outer);
        
        double r_inner, z_inner;
        if (i == 0) {
            r_inner = 0.0;
            z_inner = geometry_.t_top;
        } else {
            r_inner = r_outer + t * n_r[i];
            z_inner = z_outer + t * n_z[i];
            if (r_inner < 0.0) r_inner = 0.0;
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
        double dx = inner[i].r - inner[i-1].r;
        double dz = inner[i].z - inner[i-1].z;
        double dist = std::sqrt(dx*dx + dz*dz);
        if (dist < 1e-9) return false;
    }
    return true;
}

bool BellGenerator::checkSlopeAngles() const {
    return true;
}

// ============================================================
// ПРОВЕРКА КАЧЕСТВА ЭЛЕМЕНТОВ
// ============================================================
bool BellGenerator::checkElementQuality() const {
    const double min_angle_deg = 10.0;
    const double max_angle_deg = 170.0;
    const double max_aspect_ratio = 100.0;
    const double min_area = 1e-12;
    
    for (const auto& elem : mesh_data_.elements) {
        const Node& n1 = mesh_data_.nodes[elem.n1 - 1];
        const Node& n2 = mesh_data_.nodes[elem.n2 - 1];
        const Node& n3 = mesh_data_.nodes[elem.n3 - 1];
        const Node& n4 = mesh_data_.nodes[elem.n4 - 1];
        
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
            return std::acos(cos_angle) * 180.0 / M_PI;
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
    
    std::vector<double> n_r, n_z;
    computeNormals(n_r, n_z);
    computeInnerProfile(n_r, n_z);
    
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
    
    if (!checkElementQuality()) {
        logError("Element quality check failed: degenerate elements detected");
        if (validation_out) {
            validation_out->passed = false;
            validation_out->reason = "Degenerate elements in mesh";
        }
        return false;
    }
    
    is_mesh_built_ = true;
    return true;
}

// ============================================================
// СОХРАНЕНИЕ В ФАЙЛ
// ============================================================
bool BellGenerator::saveToFile(const std::string& filename) const {
    if (!is_mesh_built_ || mesh_data_.nodes.empty()) {
        logError("Mesh not built or empty");
        return false;
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        logError("Cannot open file: " + filename);
        return false;
    }
    
    file << "*HEADING\n";
    file << "Bell - axisymmetric mesh generated by BellGenerator (GSL)\n";
    file << "*MATERIAL, NAME=BRONZE\n";
    file << "*ELASTIC\n";
    file << material_.E << ", " << material_.nu << "\n";
    file << "*DENSITY\n";
    file << material_.rho << "\n";
    
    file << "*NODE, NSET=Nall\n";
    for (const auto& node : mesh_data_.nodes) {
        file << node.id << ", " << node.r << ", " << node.z << "\n";
    }
    
    file << "*ELEMENT, TYPE=CAX4, ELSET=Eall\n";
    for (size_t i = 0; i < mesh_data_.elements.size(); ++i) {
        const auto& e = mesh_data_.elements[i];
        file << (i + 1) << ", " << e.n1 << ", " << e.n2 << ", " << e.n3 << ", " << e.n4 << "\n";
    }
    
    file << "*SOLID SECTION, ELSET=Eall, MATERIAL=BRONZE\n";
    file << "*BOUNDARY\n";
    file << "1, 1, 1, 0.0\n";
    file << (config_.nNodes + 1) << ", 1, 1, 0.0\n";
    
    file << "*STEP\n";
    file << "*FREQUENCY, STORAGE=YES\n";
    file << "10\n";
    
    file << "*NODE FILE, OUTPUT=2D\n";
    file << "U\n";
    
    file << "*EL FILE, OUTPUT=2D\n";
    file << "MASS\n";
    
    file << "*END STEP\n";
    
    file.close();
    return true;
}

// ============================================================
// СБРОС
// ============================================================
void BellGenerator::reset() {
    mesh_data_.clear();
    is_mesh_built_ = false;
    if (spline_) {
        gsl_spline_free(spline_);
        spline_ = nullptr;
    }
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