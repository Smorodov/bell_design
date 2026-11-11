#pragma once

#include "bell_types.h"
#include "simple_json.hpp"
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>

namespace bell::profiles {
using bell::json::Value;

inline const Value& required(const Value& object, const char* key) {
    const Value* value = object.find(key);
    if (!value) throw std::runtime_error(std::string("missing required field '") + key + "'");
    return *value;
}

inline double finite_number(const Value& value, const std::string& field) {
    if (value.type() != Value::Type::Number || !std::isfinite(value.number(0.0)))
        throw std::runtime_error(field + " must be a finite number");
    return value.number(0.0);
}

inline std::size_t mode_index(const Value& value, const std::string& field) {
    if (!value.is_object()) throw std::runtime_error(field + " must be an object");
    const auto& selector = required(value, "selector");
    if (selector.type() != Value::Type::String || selector.string("") != "sorted_index")
        throw std::runtime_error(field + ".selector currently supports only 'sorted_index'");
    const double raw = finite_number(required(value, "index"), field + ".index");
    if (raw < 0.0 || std::floor(raw) != raw || raw > static_cast<double>(std::numeric_limits<std::size_t>::max()))
        throw std::runtime_error(field + ".index must be a non-negative integer");
    return static_cast<std::size_t>(raw);
}

inline void validate_targets(const TargetData& targets) {
    if (!std::isfinite(targets.prime_freq_target) || targets.prime_freq_target <= 0.0)
        throw std::runtime_error("frequency_reference.target_hz must be finite and positive");
    bool active = false;
    std::set<std::size_t> frequency_modes;
    for (const auto& t : targets.frequency_targets) {
        if (t.mode_index == targets.frequency_reference_mode_index) throw std::runtime_error("frequency target must not compare the reference mode with itself");
        if (!std::isfinite(t.ratio) || t.ratio <= 0.0) throw std::runtime_error("frequency target ratio must be finite and positive");
        if (!std::isfinite(t.weight) || t.weight < 0.0) throw std::runtime_error("frequency target weight must be finite and non-negative");
        if (!frequency_modes.insert(t.mode_index).second) throw std::runtime_error("duplicate frequency target mode selector");
        active = active || t.weight > 0.0;
    }
    std::set<std::pair<std::size_t, std::size_t>> amplitude_pairs;
    for (const auto& t : targets.amplitude_targets) {
        if (t.mode_index == t.reference_mode_index) throw std::runtime_error("amplitude target must not compare a mode with itself");
        if (!std::isfinite(t.ratio) || t.ratio <= 0.0) throw std::runtime_error("amplitude target ratio must be finite and positive");
        if (!std::isfinite(t.weight) || t.weight < 0.0) throw std::runtime_error("amplitude target weight must be finite and non-negative");
        if (!amplitude_pairs.insert({t.mode_index, t.reference_mode_index}).second) throw std::runtime_error("duplicate amplitude target selector pair");
        active = active || t.weight > 0.0;
    }
    if (!active) throw std::runtime_error("profile must contain at least one active frequency or amplitude target");
}

inline TargetData load(const Value& root) {
    if (!root.is_object()) throw std::runtime_error("target profile root must be an object");
    const auto& version = required(root, "schema_version");
    const double v = finite_number(version, "schema_version");
    if (v != 2.0) throw std::runtime_error("unsupported target profile schema_version (expected 2; legacy profiles are not supported)");
    TargetData out;
    out.profile_id = bell::json::string(root, "profile_id", "");
    out.display_name = bell::json::string(root, "display_name", out.profile_id);
    if (out.profile_id.empty()) throw std::runtime_error("profile_id must not be empty");
    const auto& ref = required(root, "frequency_reference");
    out.frequency_reference_mode_index = mode_index(required(ref, "mode"), "frequency_reference.mode");
    out.prime_freq_target = finite_number(required(ref, "target_hz"), "frequency_reference.target_hz");

    const auto& frequencies = required(root, "frequency_targets");
    if (frequencies.type() != Value::Type::Array) throw std::runtime_error("frequency_targets must be an array");
    for (std::size_t i = 0; i < frequencies.array().size(); ++i) {
        const auto& item = frequencies.array()[i];
        FrequencyTarget t;
        t.mode_index = mode_index(required(item, "mode"), "frequency_targets[" + std::to_string(i) + "].mode");
        t.ratio = finite_number(required(item, "ratio"), "frequency target ratio");
        t.weight = finite_number(required(item, "weight"), "frequency target weight");
        if (const auto* label = item.find("label")) t.label = label->string("");
        out.frequency_targets.push_back(std::move(t));
    }
    const auto& amplitudes = required(root, "amplitude_targets");
    if (amplitudes.type() != Value::Type::Array) throw std::runtime_error("amplitude_targets must be an array");
    for (std::size_t i = 0; i < amplitudes.array().size(); ++i) {
        const auto& item = amplitudes.array()[i];
        AmplitudeTarget t;
        t.mode_index = mode_index(required(item, "mode"), "amplitude_targets[" + std::to_string(i) + "].mode");
        t.reference_mode_index = mode_index(required(item, "reference_mode"), "amplitude_targets[" + std::to_string(i) + "].reference_mode");
        t.ratio = finite_number(required(item, "ratio"), "amplitude target ratio");
        t.weight = finite_number(required(item, "weight"), "amplitude target weight");
        if (const auto* label = item.find("label")) t.label = label->string("");
        out.amplitude_targets.push_back(std::move(t));
    }
    validate_targets(out);
    return out;
}

} // namespace bell::profiles
