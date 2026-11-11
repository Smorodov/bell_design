#include "../src/core/target_profile_loader.hpp"
#include <cmath>
#include <iostream>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

int main() {
    try {
        bool rejected_legacy = false;
        try {
            (void)bell::profiles::load(bell::json::Parser(R"({"prime_freq_target":440,"freq_ratios":[0.5,1.25,1.5,2.0]})").parse());
        } catch (const std::exception&) { rejected_legacy = true; }
        require(rejected_legacy, "legacy profile format must not be accepted");

        const auto v2 = bell::profiles::load(bell::json::Parser(R"({
          "schema_version":2, "profile_id":"free_form_profile", "display_name":"Free form",
          "frequency_reference":{"mode":{"selector":"sorted_index","index":3},"target_hz":523.25},
          "frequency_targets":[
            {"mode":{"selector":"sorted_index","index":0},"ratio":0.4,"weight":1,"label":"low"},
            {"mode":{"selector":"sorted_index","index":5},"ratio":2.1,"weight":0.5}
          ],
          "amplitude_targets":[
            {"mode":{"selector":"sorted_index","index":4},"reference_mode":{"selector":"sorted_index","index":1},"ratio":0.25,"weight":0}
          ]
        })").parse());
        require(v2.frequency_reference_mode_index == 3, "v2 arbitrary reference mode");
        require(v2.prime_freq_target == 523.25, "v2 anchor frequency");
        require(v2.frequency_targets.size() == 2 && v2.frequency_targets[1].mode_index == 5,
                "v2 variable target count/index");
        require(v2.amplitude_targets.size() == 1 && v2.amplitude_targets[0].reference_mode_index == 1,
                "v2 explicit amplitude reference");

        for (const std::size_t count : {1u, 4u, 5u, 7u}) {
            std::ostringstream json;
            json << "{\"schema_version\":2,\"profile_id\":\"count_test\","
                 << "\"frequency_reference\":{\"mode\":{\"selector\":\"sorted_index\",\"index\":" << count
                 << "},\"target_hz\":440},\"frequency_targets\":[";
            for (std::size_t i = 0; i < count; ++i) {
                if (i) json << ",";
                json << "{\"mode\":{\"selector\":\"sorted_index\",\"index\":" << i
                     << "},\"ratio\":" << (0.5 + 0.1 * i) << ",\"weight\":1}";
            }
            json << "],\"amplitude_targets\":[]}";
            const auto profile = bell::profiles::load(bell::json::Parser(json.str()).parse());
            require(profile.frequency_targets.size() == count, "v2 arbitrary target count");
        }

        bool rejected = false;
        try {
            (void)bell::profiles::load(bell::json::Parser(R"({
              "schema_version":2,"profile_id":"bad",
              "frequency_reference":{"mode":{"selector":"mode_family","index":1},"target_hz":440},
              "frequency_targets":[{"mode":{"selector":"sorted_index","index":0},"ratio":0.5,"weight":1}],
              "amplitude_targets":[]
            })").parse());
        } catch (const std::exception&) { rejected = true; }
        require(rejected, "unsupported selector must be rejected");

        rejected = false;
        try {
            (void)bell::profiles::load(bell::json::Parser(R"({
              "schema_version":2,"profile_id":"bad",
              "frequency_reference":{"mode":{"selector":"sorted_index","index":1},"target_hz":440},
              "frequency_targets":[],"amplitude_targets":[]
            })").parse());
        } catch (const std::exception&) { rejected = true; }
        require(rejected, "empty inactive profile must be rejected");

        const std::filesystem::path root = BELL_SOURCE_DIR;
        const std::filesystem::path bells = root / "bells";
        if (std::filesystem::exists(bells)) {
            for (const auto& entry : std::filesystem::directory_iterator(bells)) {
                const auto name = entry.path().filename().string();
                if (name.rfind("targets", 0) == 0 && entry.path().extension() == ".json")
                    (void)bell::profiles::load(bell::json::parse_file(entry.path().string()));
            }
        }
        std::cout << "Target profile loader contract passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Target profile loader contract failed: " << e.what() << "\n";
        return 1;
    }
}
