#include "bell_model_profile.hpp"

#include "bell_generator.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace BellOptimizer::Models {
namespace {

class ProfileSeedProvider final : public GeometrySeedProvider {
public:
    explicit ProfileSeedProvider(std::string family) : family_(std::move(family)) {}

    const char* family_id() const noexcept override { return family_.c_str(); }

    std::vector<BellGeometry> make_seeds(const BellMaterial& material,
                                         int requested_knots) const override {
        std::vector<BellGeometry> geoms;
        auto add = [&](double rmouth, double ttop, double tbottom,
                       std::initializer_list<double> knots) {
            BellGeometry g;
            g.H = 1000.0;
            g.R_mouth = rmouth;
            g.t_top = ttop;
            g.t_bottom = tbottom;
            g.r_knots.assign(knots.begin(), knots.end());
            geoms.push_back(g);
        };

        if (family_ == "random") return geoms;
        if (family_ == "european") {
            add(500.0, 10.0, 30.0, {0.0, 0.12, 0.34, 0.68, 1.0});
            add(560.0, 11.0, 34.0, {0.0, 0.14, 0.38, 0.72, 1.0});
            add(460.0, 9.0, 28.0, {0.0, 0.10, 0.30, 0.64, 1.0});
        } else if (family_ == "chinese") {
            // This is an axisymmetric equivalent profile prior only. It does
            // not claim to reproduce the non-axisymmetric zhong modes.
            add(560.0, 11.0, 32.0, {0.0, 0.11, 0.32, 0.66, 1.0});
            add(610.0, 12.0, 36.0, {0.0, 0.12, 0.36, 0.70, 1.0});
            add(520.0, 10.0, 30.0, {0.0, 0.09, 0.28, 0.62, 1.0});
        } else if (family_ == "orthodox" || family_ == "russian_orthodox") {
            add(540.0, 11.0, 38.0, {0.0, 0.10, 0.31, 0.65, 1.0});
            add(600.0, 12.0, 42.0, {0.0, 0.12, 0.35, 0.68, 1.0});
            add(500.0, 10.0, 36.0, {0.0, 0.09, 0.29, 0.63, 1.0});
        } else {
            add(500.0, 10.0, 30.0, {0.0, 0.15, 0.40, 0.70, 1.0});
            add(600.0, 12.0, 35.0, {0.0, 0.12, 0.35, 0.65, 1.0});
            add(400.0, 8.0, 25.0, {0.0, 0.18, 0.45, 0.75, 1.0});
            add(520.0, 11.0, 32.0, {0.0, 0.14, 0.38, 0.68, 1.0});
            add(450.0, 9.0, 27.0, {0.0, 0.16, 0.42, 0.72, 1.0});
        }

        if (requested_knots <= 0 || requested_knots == 5) return geoms;
        const int count = std::max(3, requested_knots);
        for (auto& g : geoms) {
            BellMeshConfig config;
            config.nNodes = 15;
            BellGenerator generator(material, config);
            generator.setGeometry(g);
            GeometryValidationResult validation;
            if (!generator.buildMesh(&validation)) continue;
            std::vector<double> knots;
            knots.reserve(static_cast<size_t>(count));
            for (int i = 0; i < count; ++i) {
                const double z = g.H * static_cast<double>(i) / static_cast<double>(count - 1);
                knots.push_back(std::clamp(generator.getOuterRadius(z) / g.R_mouth, 0.0, 1.0));
            }
            knots.front() = 0.0;
            knots.back() = 1.0;
            g.r_knots = std::move(knots);
        }
        return geoms;
    }

private:
    std::string family_;
};

double range_score(double x, double lo, double hi, double softness) {
    if (x >= lo && x <= hi) return 1.0;
    const double d = (x < lo) ? (lo - x) : (x - hi);
    return std::exp(-d / std::max(softness, 1e-9));
}

double thickness_score(const BellGeometry& g) {
    if (g.H <= 0.0 || g.t_top <= 0.0 || g.t_bottom <= 0.0) return 0.0;
    const double a = g.t_top / g.H;
    const double b = g.t_bottom / g.H;
    const double ratio = b / a;
    return 0.35 * range_score(ratio, 1.5, 8.0, 3.0)
         + 0.30 * range_score(a, 0.006, 0.030, 0.012)
         + 0.35 * range_score(b, 0.020, 0.090, 0.025);
}

double profile_score(const BellGeometry& g) {
    if (g.r_knots.size() < 5 || g.H <= 0.0 || g.R_mouth <= 0.0) return 0.0;
    double max_second = 0.0;
    for (size_t i = 1; i + 1 < g.r_knots.size(); ++i)
        max_second = std::max(max_second, std::abs(g.r_knots[i + 1] - 2.0 * g.r_knots[i] + g.r_knots[i - 1]));
    return std::exp(-3.0 * max_second);
}

} // namespace

std::unique_ptr<GeometrySeedProvider> make_seed_provider(const std::string& family) {
    return std::make_unique<ProfileSeedProvider>(family.empty() ? "neutral" : family);
}

bool is_supported_tuning_model(const std::string& tuning_model) {
    return tuning_model == "western_five_partial" ||
           tuning_model == "western_five_partial_reference" ||
           tuning_model == "russian_orthodox_sysoi_2000_experimental";
}

TraditionAssessment assess_tradition(const BellGeometry& g, const std::string& tradition) {
    TraditionAssessment a;
    a.tradition = tradition;
    const double diameter_height = (g.H > 0.0) ? 2.0 * g.R_mouth / g.H : 0.0;
    const double bottom_top = (g.t_top > 0.0) ? g.t_bottom / g.t_top : 0.0;

    double proportion = 0.0;
    if (tradition == "chinese" || tradition == "chinese_zhong") {
        proportion = range_score(diameter_height, 0.65, 1.15, 0.22);
        a.compatible_with_axisymmetric_model = false;
        a.limitation = "Traditional zhong requires elliptical/asymmetric cross-section and strike-location-dependent modes; the current model is axisymmetric. This score is only a silhouette prior.";
    } else if (tradition == "orthodox" || tradition == "russian_orthodox") {
        proportion = range_score(diameter_height, 0.85, 1.25, 0.25);
        a.limitation = "Russian Orthodox tradition is not defined by a single canonical modal tuning; this is a geometry/style assessment, not a liturgical or historical certification.";
    } else if (tradition == "european" || tradition == "western") {
        proportion = range_score(diameter_height, 0.75, 1.15, 0.22);
        a.limitation = "European assessment is a broad Western bell prior, not a claim about a particular founder, era, or regional school.";
    } else {
        proportion = range_score(diameter_height, 0.70, 1.25, 0.30);
        a.limitation = "Neutral geometry assessment.";
    }

    const double th = thickness_score(g);
    const double prof = profile_score(g);
    const double lower_mass = range_score(bottom_top, 2.0, 7.0, 2.5);
    a.proportion_score = proportion;
    a.thickness_score = 0.7 * th + 0.3 * lower_mass;
    a.profile_score = prof;
    a.geometry_score = 100.0 * (0.45 * proportion + 0.30 * a.thickness_score + 0.25 * prof);

    if (!a.compatible_with_axisymmetric_model) {
        if (a.geometry_score >= 80.0) a.verdict = "strong silhouette compatibility; acoustic tradition not representable";
        else if (a.geometry_score >= 65.0) a.verdict = "reasonable silhouette compatibility; acoustic tradition not representable";
        else a.verdict = "weak silhouette compatibility; acoustic tradition not representable";
    } else if (a.geometry_score >= 80.0) a.verdict = "strong geometric/style compatibility";
    else if (a.geometry_score >= 65.0) a.verdict = "reasonable compatibility";
    else if (a.geometry_score >= 50.0) a.verdict = "mixed compatibility; inspect profile visually";
    else a.verdict = "weak compatibility; geometry is outside the selected style prior";
    return a;
}

} // namespace BellOptimizer::Models
