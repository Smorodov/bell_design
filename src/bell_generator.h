// bell_generator.h
#ifndef BELL_GENERATOR_H
#define BELL_GENERATOR_H

#include "bell_types.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

// ============================================================
// КЛАСС ГЕНЕРАТОРА (переиспользуемый)
// ============================================================
class BellGenerator {
public:
    // Конструктор — создаётся один раз в начале оптимизации
    BellGenerator(const BellMaterial& material, const BellMeshConfig& config);
    ~BellGenerator();

    // Запрет копирования (для многопоточности — каждый поток свой объект)
    BellGenerator(const BellGenerator&) = delete;
    BellGenerator& operator=(const BellGenerator&) = delete;

    // ============================================================
    // УПРАВЛЕНИЕ ГЕОМЕТРИЕЙ
    // ============================================================
    
    // Обновить геометрию для нового кандидата
    void setGeometry(const BellGeometry& geom);
    
    // Получить текущую геометрию
    const BellGeometry& getGeometry() const { return geometry_; }
    
    // ============================================================
    // ЭТАП 1: БЫСТРАЯ ГЕОМЕТРИЧЕСКАЯ ПРОВЕРКА (БЕЗ FEM)
    // ============================================================
    
    // Проверка без построения сетки (доли секунды)
    GeometryValidationResult validateGeometry() const;
    
    // ============================================================
    // ЭТАП 2: ПОСТРОЕНИЕ СЕТКИ
    // ============================================================
    
    // Построить сетку (после успешной валидации)
    bool buildMesh();
    
    // Построить сетку с автоматической валидацией
    bool buildMesh(GeometryValidationResult* validation_out);
    
    // ============================================================
    // ДОСТУП К ДАННЫМ В ПАМЯТИ
    // ============================================================
    
    const BellMeshData& getMeshData() const { return mesh_data_; }
    
    // Удобные методы для быстрого доступа
    const std::vector<Node>& getNodes() const { return mesh_data_.nodes; }
    const std::vector<Element>& getElements() const { return mesh_data_.elements; }
    const std::vector<Point2D>& getOuterProfile() const { return mesh_data_.outer_profile; }
    const std::vector<Point2D>& getInnerProfile() const { return mesh_data_.inner_profile; }
    
    // ============================================================
    // МЕТРИКИ
    // ============================================================
    
    double getTotalArcLength() const { return mesh_data_.total_arc_length; }
    double getMinThickness() const { return mesh_data_.min_thickness; }
    double getMaxThickness() const { return mesh_data_.max_thickness; }
    double getThicknessAt(double z) const;
    
    // Радиусы в произвольной точке
    double getOuterRadius(double z) const;
    double getInnerRadius(double z) const;
    
    // ============================================================
    // ПРОИЗВОДНАЯ (ПУБЛИЧНАЯ — ДЛЯ ДИАГНОСТИКИ И АНАЛИЗА)
    // ============================================================
    double derivative(double z) const;
    
    // ============================================================
    // МЕТОД ДЛЯ ИНТЕГРАТОРА (ПУБЛИЧНЫЙ)
    // ============================================================
    double integrandParametric(double t) const;
    
    // ============================================================
    // СОХРАНЕНИЕ (только для финального результата)
    // ============================================================
    bool saveToFile(const std::string& filename) const;
    
    // ============================================================
    // СБРОС
    // ============================================================
    void reset();  // очистить данные, но сохранить материал и конфиг

private:
    // ============================================================
    // ВНУТРЕННИЕ МЕТОДЫ
    // ============================================================
    
    // Построение сплайна по текущей геометрии
    bool buildSpline();
    
    // Вычисление длины дуги
    double computeArcLength();
    double computeArcLengthToT(double t_target) const;
    
    // Построение профилей
    void computeOuterProfile();
    void computeNormals(std::vector<double>& n_r, std::vector<double>& n_z) const;
    void computeInnerProfile(const std::vector<double>& n_r, const std::vector<double>& n_z);
    void assembleMesh();
    void computeMetrics();
    
    // ============================================================
    // ГЕОМЕТРИЧЕСКИЕ ПРОВЕРКИ
    // ============================================================
    
    bool checkMonotonicity() const;
    bool checkSelfIntersection() const;
    bool checkSlopeAngles() const;
    bool checkElementQuality() const;   // <-- НОВЫЙ МЕТОД
    
    // ============================================================
    // ВСПОМОГАТЕЛЬНЫЕ МЕТОДЫ
    // ============================================================
    
    double radius(double z) const;
    double getZfromS(double s) const;
    
    // ============================================================
    // ЛОГИРОВАНИЕ
    // ============================================================
    
    void logWarning(const std::string& msg) const;
    void logError(const std::string& msg) const;

    // ============================================================
    // ДАННЫЕ
    // ============================================================
    
    BellMaterial material_;
    BellMeshConfig config_;
    BellGeometry geometry_;
    BellMeshData mesh_data_;
    bool is_mesh_built_ = false;
    
    // GSL объекты (переиспользуются)
    gsl_spline* spline_ = nullptr;
    gsl_interp_accel* acc_ = nullptr;
    mutable gsl_integration_workspace* integ_w_ = nullptr;
};

#endif // BELL_GENERATOR_H