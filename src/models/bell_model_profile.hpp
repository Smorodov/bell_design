#pragma once

#include "bell_types.h"

#include <memory>
#include <string>
#include <vector>

namespace BellOptimizer::Models {

// Concrete bell traditions belong to this layer. Physics/FEM code must not
// branch on these identifiers.
struct BellModelSelection {
    std::string family = "neutral";
    std::string tuning_model = "western_five_partial";
};

struct TraditionAssessment {
    std::string tradition;
    double geometry_score = 0.0;
    double proportion_score = 0.0;
    double thickness_score = 0.0;
    double profile_score = 0.0;
    bool compatible_with_axisymmetric_model = true;
    std::string verdict;
    std::string limitation;
};

class GeometrySeedProvider {
public:
    virtual ~GeometrySeedProvider() = default;
    virtual const char* family_id() const noexcept = 0;
    virtual std::vector<BellGeometry> make_seeds(const BellMaterial& material,
                                                   int requested_knots) const = 0;
};

std::unique_ptr<GeometrySeedProvider> make_seed_provider(const std::string& family);

bool is_supported_tuning_model(const std::string& tuning_model);

TraditionAssessment assess_tradition(const BellGeometry& geometry,
                                     const std::string& tradition);

} // namespace BellOptimizer::Models
