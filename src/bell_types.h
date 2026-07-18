// bell_types.h
#ifndef BELL_TYPES_H
#define BELL_TYPES_H

#include <vector>
#include <string>
#include <array>
#include <cmath>

// ============================================================
// КОНСТАНТЫ (вынесенные магические числа)
// ============================================================
namespace BellConstants {
    const double ERROR_INVALID = 1e9;
    const double ERROR_THRESHOLD_VALID = 1e8;
    
    const double MIN_FREQ_HZ = 1.0;
    const double MIN_SCALE_FACTOR = 0.05;
    const double MAX_SCALE_FACTOR = 20.0;
    
    const double MIN_H = 200.0;
    const double MAX_H = 4000.0;
    const double MIN_RMOUTH_RATIO = 0.4;
    const double MAX_RMOUTH_RATIO = 0.8;
    const double MIN_TTOP_RATIO = 0.005;
    const double MAX_TTOP_RATIO = 0.03;
    const double MIN_TBOTTOM_RATIO = 0.02;
    const double MAX_TBOTTOM_RATIO = 0.09;
    const double MIN_TTOP_ABS = 2.0;
    const double MAX_TTOP_ABS = 60.0;
    const double MIN_TBOTTOM_ABS = 10.0;
    const double MAX_TBOTTOM_ABS = 160.0;
    
    const double MAX_SLOPE_ANGLE_DEG = 45.0;
    inline double getMaxSlope() {
        const double PI = 3.14159265358979323846264338328;
        return std::tan(MAX_SLOPE_ANGLE_DEG * PI / 180.0);
    }
    
    const double DIAMETER_RATIO_MIN = 0.7;
    const double DIAMETER_RATIO_MAX = 1.3;
    
    const double REFERENCE_FREQ_HZ = 440.0;
    
    const double ADAPTIVE_RADIUS_MIN = 0.03;
    const double ADAPTIVE_RADIUS_MAX = 0.25;
    const double ADAPTIVE_RADIUS_ERROR_SCALE = 1.0;
    
    const double MUTATION_RADIUS_KNOTS = 0.15;
    const double MUTATION_RADIUS_RATIO = 0.05;
    const double MUTATION_RADIUS_H = 0.05;
    const double MUTATION_RADIUS_RMOUTH = 0.05;
    const double MUTATION_RADIUS_TTOP = 0.005;
    const double MUTATION_RADIUS_TBOTTOM = 0.01;
    
    const double KNOT_FRACTION_MIN_BASE = 0.05;
    const double KNOT_FRACTION_MAX_BASE = 0.15;
    const double KNOT_FRACTION_STEP = 0.25;
    const double KNOT_FRACTION_STEP_MAX = 0.28;
    
    const double SIGMOID_K_ANGLE = 2.0;
    const double SIGMOID_ANGLE_WEIGHT = 0.5;
    const double SIGMOID_K_PROPORTION = 2.0;
    const double SIGMOID_PROPORTION_WEIGHT = 0.5;
}

// ============================================================
// ПАРАМЕТРЫ ГЕОМЕТРИИ
// ============================================================
struct BellGeometry {
    double H = 1000.0;
    double R_mouth = 400.0;
    double t_top = 10.0;
    double t_bottom = 60.0;
    std::vector<double> r_knots;
    
    void scale(double factor) {
        H *= factor;
        R_mouth *= factor;
        t_top *= factor;
        t_bottom *= factor;
        // r_knots остаются нормированными
    }
    
    bool isValid() const {
        if (H <= 0.0) return false;
        if (R_mouth <= 0.0) return false;
        if (t_top <= 0.0) return false;
        if (t_bottom <= 0.0) return false;
        if (r_knots.size() < 3) return false;
        for (double r : r_knots) {
            if (r < 0.0 || r > 1.0) return false;
        }
        return true;
    }
    
    size_t numKnots() const { return r_knots.size(); }
    
    double getKnotZ(size_t i) const {
        if (r_knots.empty()) return 0.0;
        return i * H / (r_knots.size() - 1);
    }
};

// ============================================================
// ПАРАМЕТРЫ МАТЕРИАЛА
// ============================================================
struct BellMaterial {
    double E = 1.1e4;
    double nu = 0.34;
    double rho = 8.6e-9;
};

// ============================================================
// КОНФИГУРАЦИЯ СЕТКИ
// ============================================================
struct BellMeshConfig {
    int nNodes = 15;
};

// ============================================================
// ТОЧКИ И УЗЛЫ
// ============================================================
struct Point2D {
    double r = 0.0;
    double z = 0.0;
};

struct Node {
    int id = 0;
    double r = 0.0;
    double z = 0.0;
};

struct Element {
    int n1 = 0, n2 = 0, n3 = 0, n4 = 0;
};

// ============================================================
// РЕЗУЛЬТАТЫ ПОСТРОЕНИЯ СЕТКИ
// ============================================================
struct BellMeshData {
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<Point2D> outer_profile;
    std::vector<Point2D> inner_profile;
    double total_arc_length = 0.0;
    double min_thickness = 0.0;
    double max_thickness = 0.0;
    
    void clear() {
        nodes.clear();
        elements.clear();
        outer_profile.clear();
        inner_profile.clear();
        total_arc_length = 0.0;
        min_thickness = 0.0;
        max_thickness = 0.0;
    }
};

// ============================================================
// РЕЗУЛЬТАТЫ ГЕОМЕТРИЧЕСКОЙ ПРОВЕРКИ
// ============================================================
struct GeometryValidationResult {
    bool passed = true;
    std::string reason = "";
    double penalty_value = 0.0;
    double max_slope_angle = 0.0;
};

// ============================================================
// ЦЕЛЕВЫЕ ПАРАМЕТРЫ (ДЛЯ ПРАВОСЛАВНОГО КОЛОКОЛА)
// ============================================================
struct TargetData {
    // Целевая частота Prime (второй моды)
    double prime_freq_target = 440.0;
    
    // Целевые соотношения частот относительно Prime:
    // [Hum/Prime, Tertia/Prime, Quinta/Prime, Nominal/Prime]
    std::array<double, 4> freq_ratios = {0.5, 1.25, 1.5, 2.0};
    
    // Целевые соотношения амплитуд относительно Hum:
    // [Prime/Hum, Tertia/Hum, Quinta/Hum, Nominal/Hum]
    std::array<double, 4> amp_ratios = {0.8, 0.6, 0.4, 0.2};
    
    // Веса ошибок
    double prime_weight = 1.0;
    std::array<double, 4> freq_weights = {1.0, 1.0, 1.0, 1.0};
    std::array<double, 4> amp_weights = {1.0, 1.0, 1.0, 1.0};
    
    // Крутизна сигмоидных функций
    double k_prime = 20.0;
    double k_freq = 20.0;
    double k_amp = 5.0;
};

// ============================================================
// ПАРАМЕТРЫ ЭВОЛЮЦИИ
// ============================================================
struct EvolutionParams {
    int population_size = 20;
    int generations = 10;
    int num_elite = 5;
    double mutation_probability = 0.1;
    double mutation_strength = 0.05;
    int num_knots = 5;
};

#endif // BELL_TYPES_H