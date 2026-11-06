#include "synthesize_contact.hpp"
#include "../core/bell_types.h"
#include "../core/bell_material.h"
#include "../core/bell_generator.h"
#include "../core/modal_types.hpp"
#include "../core/simple_json.hpp"
#include "../fem/cylindrical_bell_model.hpp"
#include "../fem/cylindrical_modal_evaluator.hpp"
#include "../excitation/hammer_contact.hpp"
#include "../excitation/clapper_dynamics.hpp"
#include "../excitation/modal_decay.hpp"
#include "../excitation/meridional_influence.hpp"
#include "../acoustics/directional_radiation.hpp"
#include "../audio/wav_writer.hpp"
#include "../audio/modal_pressure_synthesis.hpp"
#include "../acoustics/listener_geometry.hpp"
#include "../acoustics/finite_distance_radiation.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace BellOptimizer;
constexpr double kPi = 3.1415926535897932384626433832795;


BellGeometry loadGeometry(const std::string& filename) {
    const auto data = bell::json::parse_file(filename);
    const auto* g = data.find("geometry");
    if (!g || !g->is_object()) throw std::runtime_error("geometry object is missing");
    BellGeometry out;
    out.H = bell::json::number(*g, "H", 1000.0);
    out.R_mouth = bell::json::number(*g, "R_mouth", 500.0);
    out.t_top = bell::json::number(*g, "t_top", 10.0);
    out.t_bottom = bell::json::number(*g, "t_bottom", 30.0);
    out.ellipticity = bell::json::number(*g, "ellipticity", 1.0);
    const auto* knots = g->find("r_knots");
    if (!knots || !knots->is_array()) throw std::runtime_error("r_knots are missing");
    for (const auto& v : knots->array()) out.r_knots.push_back(v.number(0.0));
    return out;
}

struct RenderMode { double omega, alpha, q0, v0; std::complex<double> radiation; };

double fourierValue(const std::vector<double>& values, std::size_t point,
                    std::size_t pointCount, const std::vector<int>& labels, double angle) {
    if (values.size() != pointCount * labels.size()) return 0.0;
    double result = 0.0;
    for (std::size_t b = 0; b < labels.size(); ++b) {
        const int label = labels[b];
        const double coefficient = values[point * labels.size() + b];
        if (label == 0) result += coefficient;
        else if (label > 0) result += coefficient * std::cos(label * angle);
        else result += coefficient * std::sin(-label * angle);
    }
    return result;
}


}

int command_synthesize_contact(int argc, char* argv[]) {
    std::string geometryFile = "best_geometry.json", output = "bell_contact.wav", workdir;
    double strikeHeight = 0.72, strikerMass = 0.045, strikerVelocity = 1.2;
    double hertzStiffness = 2.0e7, contactStep = 2.0e-6, contactDuration = 0.025;
    double t60 = 5.0, audioDuration = 6.0;
    double strikeAzimuthDeg = 0.0, listenerPolarDeg = 30.0, listenerAzimuthDeg = 0.0;
        int sampleRate = 48000, meshNodes = 71;
    bool geometrySet = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(std::string(name) + " requires a value");
            return argv[++i];
        };
        if (a == "--help" || a == "-h") {
            std::cout << "Usage: bell synth-contact [geometry.json] [--output file.wav] [--strike-height 0..1] "
                         "[--t60 seconds] [--duration seconds] [--sample-rate hz] [--mesh-nodes n] [--strike-azimuth-deg deg] [--listener-polar-deg deg] [--listener-azimuth-deg deg] [--workdir dir]\n"
                         "Fourier-resolved localized strike and directional far-field modal radiation; output is peak-normalized, not calibrated SPL.\n";
            return 0;
        } else if (a == "--output") output = value("--output");
        else if (a == "--workdir") workdir = value("--workdir");
        else if (a == "--strike-height") strikeHeight = std::stod(value("--strike-height"));
        else if (a == "--t60") t60 = std::stod(value("--t60"));
        else if (a == "--duration") audioDuration = std::stod(value("--duration"));
        else if (a == "--sample-rate") sampleRate = std::stoi(value("--sample-rate"));
        else if (a == "--mesh-nodes") meshNodes = std::stoi(value("--mesh-nodes"));
        else if (a == "--strike-azimuth-deg") strikeAzimuthDeg = std::stod(value("--strike-azimuth-deg"));
        else if (a == "--listener-polar-deg") listenerPolarDeg = std::stod(value("--listener-polar-deg"));
        else if (a == "--listener-azimuth-deg") listenerAzimuthDeg = std::stod(value("--listener-azimuth-deg"));
        else if (a == "--striker-mass") strikerMass = std::stod(value("--striker-mass"));
        else if (a == "--striker-velocity") strikerVelocity = std::stod(value("--striker-velocity"));
        else if (a == "--hertz-stiffness") hertzStiffness = std::stod(value("--hertz-stiffness"));
        else if (a == "--contact-step") contactStep = std::stod(value("--contact-step"));
        else if (a == "--contact-duration") contactDuration = std::stod(value("--contact-duration"));
        else if (a.rfind("--", 0) == 0) { std::cerr << "Unknown option: " << a << "\n"; return 2; }
        else if (!geometrySet) { geometryFile = a; geometrySet = true; }
        else { std::cerr << "Unexpected argument: " << a << "\n"; return 2; }
    }
    try {
        if (!workdir.empty()) {
            std::error_code ec; std::filesystem::current_path(workdir, ec);
            if (ec) throw std::runtime_error("cannot change workdir: " + ec.message());
        }
        if (!std::isfinite(strikeHeight) || strikeHeight < 0.05 || strikeHeight > 0.98 ||
            !std::isfinite(t60) || t60 <= 0.0 || !std::isfinite(audioDuration) || audioDuration <= 0.1 ||
            sampleRate < 8000 || sampleRate > 192000 || meshNodes < 11 || meshNodes > 101 ||
            audioDuration * sampleRate > 10000000.0 ||
            !std::isfinite(strikeAzimuthDeg) || !std::isfinite(listenerPolarDeg) ||
            !std::isfinite(listenerAzimuthDeg) || listenerPolarDeg < 0.0 || listenerPolarDeg > 180.0)
            throw std::invalid_argument("invalid strike height, decay, duration, sample rate, or mesh size");

        const BellGeometry geometry = loadGeometry(geometryFile);
        if (geometry.ellipticity < 0.999999)
            throw std::runtime_error("synth-contact currently requires circular geometry (ellipticity=1); non-axisymmetric point strike is not represented by this solver");
        BellMaterial material;
        BellMeshConfig meshConfig; meshConfig.nNodes = meshNodes;
        BellGenerator generator(material, meshConfig);
        generator.setGeometry(geometry);
        GeometryValidationResult validation;
        if (!generator.buildMesh(&validation)) throw std::runtime_error("geometry rejected: " + validation.reason);

        const auto& profile = generator.getOuterProfile();
        if (profile.size() < 2) throw std::runtime_error("outer profile is incomplete");
        BellMeshConfig directMeshConfig; directMeshConfig.nNodes = meshNodes;
        const auto directModel = FEM::make_cylindrical_bell_model(geometry, material, directMeshConfig);
        Excitation::MeridionalInfluence influence;
        FEM::CylindricalModalEvaluationConfig femConfig;
        femConfig.harmonics = {0, 1, 2, 3};
        femConfig.angular_samples = 64;
        femConfig.requested_modes = 24;
        femConfig.max_iterations = 5000;
        femConfig.tolerance = 1.0e-6;
        femConfig.minimum_frequency_hz = 1.0;
        ModalData modalData;
        if (!FEM::evaluate_cylindrical_bell_modes(directModel, influence.make(directModel.mesh_si), femConfig, modalData))
            throw std::runtime_error("Fourier cylindrical FEM modal solve failed");
        if (modalData.modes.size() < 3) throw std::runtime_error("FEM returned fewer than three modes");

        const double targetZ = geometry.H * strikeHeight;
        std::size_t strikeIndex = 0;
        for (std::size_t i = 1; i < profile.size(); ++i)
            if (std::abs(profile[i].z - targetZ) < std::abs(profile[strikeIndex].z - targetZ)) strikeIndex = i;
        const std::size_t lo = strikeIndex == 0 ? 0 : strikeIndex - 1;
        const std::size_t hi = std::min(strikeIndex + 1, profile.size() - 1);
        const double dr = profile[hi].r - profile[lo].r;
        const double dz = profile[hi].z - profile[lo].z;
        const double tangentLength = std::hypot(dr, dz);
        if (!(tangentLength > 0.0)) throw std::runtime_error("invalid strike surface tangent");
        const double nr = dz / tangentLength, nz = -dr / tangentLength;

        const double strikeAzimuth = strikeAzimuthDeg * kPi / 180.0;
        const double listenerPolar = listenerPolarDeg * kPi / 180.0;
        const double listenerAzimuth = listenerAzimuthDeg * kPi / 180.0;
        std::vector<Excitation::ContactModalMode> contactModes;
        std::vector<std::complex<double>> radiationCoefficients;
        std::vector<double> modeOmega;
        for (const auto& m : modalData.modes) {
            if (m.outer_surface_harmonics.empty() || m.outer_surface_fourier_radial.size() != profile.size() * m.outer_surface_harmonics.size() ||
                m.outer_surface_fourier_axial.size() != profile.size() * m.outer_surface_harmonics.size()) continue;
            const double ur = fourierValue(m.outer_surface_fourier_radial, strikeIndex, profile.size(), m.outer_surface_harmonics, strikeAzimuth);
            const double uz = fourierValue(m.outer_surface_fourier_axial, strikeIndex, profile.size(), m.outer_surface_harmonics, strikeAzimuth);
            const double phiNormal = ur * nr + uz * nz;
            if (!std::isfinite(phiNormal) || std::abs(phiNormal) < 1.0e-14 || !(m.modal_mass > 0.0)) continue;
            // Direct cylindrical FEM uses SI geometry/material and mass-normalized modes.
            contactModes.push_back({m.frequency_hz, m.modal_mass, phiNormal});
            const auto directionalField = Acoustics::buildDirectionalRadiationField(m, profile, listenerPolar, 1.0);
            radiationCoefficients.push_back(directionalField.evaluate(listenerAzimuth));
            modeOmega.push_back(2.0 * kPi * m.frequency_hz);
        }
        if (contactModes.size() < 3) throw std::runtime_error("insufficient Fourier modes have strike participation");

        Excitation::HammerContactConfig contact;
        contact.striker_mass_kg = strikerMass;
        contact.initial_velocity_m_per_s = strikerVelocity;
        contact.hertz_stiffness_n_per_m_1_5 = hertzStiffness;
        contact.integration_step_seconds = contactStep;
        contact.duration_seconds = contactDuration;
        Excitation::HammerContactResult impact;
        if (!Excitation::simulateHertzHammerContact(contact, contactModes, impact))
            throw std::runtime_error("Hertz contact integration failed");

        std::size_t release = impact.contact_force_n.size();
        bool hadContact = false;
        for (std::size_t i = 0; i < impact.contact_force_n.size(); ++i) {
            if (impact.contact_force_n[i] > 1.0e-8) hadContact = true;
            else if (hadContact) { release = i; break; }
        }
        if (release >= impact.contact_force_n.size())
            throw std::runtime_error("hammer did not separate within --contact-duration; increase it or review contact parameters");

        std::vector<RenderMode> renderModes;
        renderModes.reserve(contactModes.size());
        for (std::size_t i = 0; i < contactModes.size(); ++i) {
            const double alpha = std::log(1000.0) / t60;
            const double omega = modeOmega[i];
            const double wd2 = omega * omega - alpha * alpha;
            if (!(wd2 > 0.0)) continue;
            if (std::abs(radiationCoefficients[i]) < 1.0e-18) continue;
            renderModes.push_back({omega, alpha, impact.modal_displacement[i][release],
                                   impact.modal_velocity[i][release], radiationCoefficients[i]});
        }
        if (renderModes.size() < 3) throw std::runtime_error("fewer than three underdamped acoustic modes remain");

        const std::size_t count = static_cast<std::size_t>(std::llround(audioDuration * sampleRate));
        std::vector<double> audio(count, 0.0);
        const double releaseTime = static_cast<double>(release) * contactStep;
        double peak = 0.0;
        for (std::size_t n = 0; n < count; ++n) {
            const double t = static_cast<double>(n) / sampleRate;
            const double tau = t - releaseTime;
            if (tau < 0.0) continue;
            double pressureProxy = 0.0;
            for (const auto& m : renderModes) {
                const double wd = std::sqrt(m.omega * m.omega - m.alpha * m.alpha);
                const double b = (m.v0 + m.alpha * m.q0) / wd;
                const std::complex<double> analyticQddot =
                    std::complex<double>(m.q0, -b) * std::exp(std::complex<double>(-m.alpha, wd) * tau) *
                    std::pow(std::complex<double>(-m.alpha, wd), 2);
                pressureProxy += std::real(m.radiation * analyticQddot);
            }
            audio[n] = pressureProxy;
            peak = std::max(peak, std::abs(pressureProxy));
        }
        if (!(peak > 1.0e-30) || !std::isfinite(peak)) throw std::runtime_error("acoustic radiation proxy is numerically zero");
        const double scale = 0.95 / peak;
        for (double& x : audio) x *= scale;
        if (!Audio::writeMonoPcm16(output, audio, sampleRate)) throw std::runtime_error("failed to write WAV file");

        std::filesystem::path metadataPath(output);
        metadataPath.replace_extension();
        metadataPath += "_metadata.txt";
        std::ofstream meta(metadataPath);
        if (meta) {
            meta << "BellOptimizer contact-driven modal audio\n"
                 << "Geometry: " << geometryFile << "\nMesh nodes: " << meshNodes << "\n"
                 << "Modes used: " << renderModes.size() << "\nStrike height fraction: " << strikeHeight << "\n"
                 << "Strike azimuth deg: " << strikeAzimuthDeg << "\nListener polar deg from bell axis: " << listenerPolarDeg
                 << "\nListener azimuth deg: " << listenerAzimuthDeg << "\n"
                 << "Striker mass kg: " << strikerMass << "\nStriker velocity m/s: " << strikerVelocity << "\n"
                 << "Hertz stiffness N/m^1.5: " << hertzStiffness << "\nContact duration s: " << impact.contact_duration_seconds << "\n"
                 << "Peak force N: " << *std::max_element(impact.contact_force_n.begin(), impact.contact_force_n.end()) << "\n"
                 << "Contact energy relative error: " << impact.maximum_relative_energy_error << "\n"
                 << "T60 s: " << t60 << "\nSample rate Hz: " << sampleRate << "\n"
                 << "Output normalization: peak 0.95 (not calibrated SPL)\n"
                 << "Excitation scope: full real Fourier basis m=0..3 evaluated at strike azimuth " << strikeAzimuthDeg << " deg\n"
                 << "Radiation: directional far-field surface integral at elevation " << listenerPolarDeg << " deg from bell axis, azimuth " << listenerAzimuthDeg << " deg; linear acoustic far-field approximation, not calibrated SPL\n"
                 << "Modal frequencies Hz:";
            for (const auto& m : modalData.modes) meta << " " << m.frequency_hz;
            meta << "\n";
        }
        std::cout << std::fixed << std::setprecision(4)
                  << "Generated " << output << " | modes=" << renderModes.size()
                  << " | contact=" << impact.contact_duration_seconds << " s"
                  << " | peak force=" << *std::max_element(impact.contact_force_n.begin(), impact.contact_force_n.end())
                  << " N | energy error=" << std::scientific << impact.maximum_relative_energy_error << std::fixed << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "synth-contact error: " << e.what() << "\n";
        return 1;
    }
}


int command_synthesize_clapper(int argc, char* argv[]) {
    std::string geometryFile = "best_geometry.json", output = "bell_clapper.wav", workdir;
    double clapperMass = 5.0, comLength = 0.46, tipLength = 0.52, tipRadius = 0.035;
    double inertia = 0.0, pivotZ = std::numeric_limits<double>::quiet_NaN(), pivotR = 0.0;
    double initialAngleDeg = 0.0, initialAngularVelocity = 3.2;
    double tipYoung = 200.0e9, tipPoisson = 0.29, hertzStiffness = 0.0;
    double contactStep = 2.0e-6, impactHorizon = 0.25;
    double t60 = 5.0, structuralLossFactor = std::numeric_limits<double>::quiet_NaN(), audioDuration = 6.0;
    std::string decayModelName = "fixed-t60", materialPresetName = "bell_bronze_reference";
    bool lossFactorOverride = false;
    double strikeAzimuthDeg = 0.0, listenerPolarDeg = 30.0, listenerAzimuthDeg = 0.0;
    double listenerDistance = 3.0, earSpacing = 0.18, headYawDeg = 0.0;
    bool stereoOutput = false;
    std::string stereoRadiationModel = "far-field";
    int sampleRate = 48000, meshNodes = 71;
    bool geometrySet = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(std::string(name) + " requires a value");
            return argv[++i];
        };
        if (a == "--help" || a == "-h") {
            std::cout << "Usage: bell synth-clapper [geometry.json] [--output file.wav] [--clapper-mass kg] "
                         "[--clapper-inertia kg-m2] [--com-length m] [--tip-length m] [--tip-radius m] "
                         "[--pivot-z m] [--initial-angle-deg deg] [--initial-angular-velocity rad/s] "
                         "[--contact-step s] [--impact-horizon s] [--material-preset id] [--t60 s] [--decay-model fixed-t60|structural-loss] [--loss-factor eta] [--duration s] "
                         "[--strike-azimuth-deg deg] [--listener-polar-deg deg] [--listener-azimuth-deg deg] [--stereo] [--stereo-radiation far-field|finite-distance] [--listener-distance m] [--ear-spacing m] [--head-yaw-deg deg]\n"
                         "Planar physical pendulum + internal spherical Hertz contact + Fourier FEM + far-field modal radiation.\n";
            return 0;
        } else if (a == "--output") output = value("--output");
        else if (a == "--workdir") workdir = value("--workdir");
        else if (a == "--mesh-nodes") meshNodes = std::stoi(value("--mesh-nodes"));
        else if (a == "--clapper-mass") clapperMass = std::stod(value("--clapper-mass"));
        else if (a == "--clapper-inertia") inertia = std::stod(value("--clapper-inertia"));
        else if (a == "--com-length") comLength = std::stod(value("--com-length"));
        else if (a == "--tip-length") tipLength = std::stod(value("--tip-length"));
        else if (a == "--tip-radius") tipRadius = std::stod(value("--tip-radius"));
        else if (a == "--pivot-z") pivotZ = std::stod(value("--pivot-z"));
        else if (a == "--pivot-r") pivotR = std::stod(value("--pivot-r"));
        else if (a == "--initial-angle-deg") initialAngleDeg = std::stod(value("--initial-angle-deg"));
        else if (a == "--initial-angular-velocity") initialAngularVelocity = std::stod(value("--initial-angular-velocity"));
        else if (a == "--tip-young-modulus") tipYoung = std::stod(value("--tip-young-modulus"));
        else if (a == "--tip-poisson") tipPoisson = std::stod(value("--tip-poisson"));
        else if (a == "--hertz-stiffness") hertzStiffness = std::stod(value("--hertz-stiffness"));
        else if (a == "--contact-step") contactStep = std::stod(value("--contact-step"));
        else if (a == "--impact-horizon") impactHorizon = std::stod(value("--impact-horizon"));
        else if (a == "--material-preset") materialPresetName = value("--material-preset");
        else if (a == "--t60") t60 = std::stod(value("--t60"));
        else if (a == "--decay-model") decayModelName = value("--decay-model");
        else if (a == "--loss-factor") { structuralLossFactor = std::stod(value("--loss-factor")); lossFactorOverride = true; }
        else if (a == "--duration") audioDuration = std::stod(value("--duration"));
        else if (a == "--sample-rate") sampleRate = std::stoi(value("--sample-rate"));
        else if (a == "--strike-azimuth-deg") strikeAzimuthDeg = std::stod(value("--strike-azimuth-deg"));
        else if (a == "--listener-polar-deg") listenerPolarDeg = std::stod(value("--listener-polar-deg"));
        else if (a == "--listener-azimuth-deg") listenerAzimuthDeg = std::stod(value("--listener-azimuth-deg"));
        else if (a == "--stereo") stereoOutput = true;
        else if (a == "--stereo-radiation") stereoRadiationModel = value("--stereo-radiation");
        else if (a == "--listener-distance") listenerDistance = std::stod(value("--listener-distance"));
        else if (a == "--ear-spacing") earSpacing = std::stod(value("--ear-spacing"));
        else if (a == "--head-yaw-deg") headYawDeg = std::stod(value("--head-yaw-deg"));
        else if (a.rfind("--", 0) == 0) { std::cerr << "Unknown option: " << a << "\n"; return 2; }
        else if (!geometrySet) { geometryFile = a; geometrySet = true; }
        else { std::cerr << "Unexpected argument: " << a << "\n"; return 2; }
    }
    try {
        if (!workdir.empty()) {
            std::error_code ec; std::filesystem::current_path(workdir, ec);
            if (ec) throw std::runtime_error("cannot change workdir: " + ec.message());
        }
        if (!std::isfinite(clapperMass) || clapperMass <= 0.0 || !std::isfinite(comLength) || comLength <= 0.0 ||
            !std::isfinite(tipLength) || tipLength <= 0.0 || !std::isfinite(tipRadius) || tipRadius <= 0.0 ||
            !std::isfinite(pivotR) || pivotR < 0.0 || !std::isfinite(initialAngleDeg) ||
            !std::isfinite(initialAngularVelocity) || !std::isfinite(tipYoung) || tipYoung <= 0.0 ||
            !std::isfinite(tipPoisson) || tipPoisson <= -1.0 || tipPoisson >= 0.5 ||
            !std::isfinite(contactStep) || contactStep <= 0.0 || !std::isfinite(impactHorizon) || impactHorizon <= 0.0 ||
            !std::isfinite(t60) || t60 <= 0.0 || (lossFactorOverride && !std::isfinite(structuralLossFactor)) ||
            !std::isfinite(audioDuration) || audioDuration <= 0.1 ||
            sampleRate < 8000 || sampleRate > 192000 || meshNodes < 11 || meshNodes > 101 ||
            audioDuration * sampleRate > 10000000.0 || !std::isfinite(strikeAzimuthDeg) ||
            !std::isfinite(listenerPolarDeg) || listenerPolarDeg < 0.0 || listenerPolarDeg > 180.0 ||
            !std::isfinite(listenerAzimuthDeg) || !std::isfinite(listenerDistance) || listenerDistance <= 0.25 ||
            !std::isfinite(earSpacing) || earSpacing <= 0.0 || earSpacing >= 0.5 ||
            !std::isfinite(headYawDeg) || headYawDeg < -90.0 || headYawDeg > 90.0)
            throw std::invalid_argument("invalid clapper, contact, audio, listener distance, ear spacing, or mesh parameter");
        if (stereoRadiationModel != "far-field" && stereoRadiationModel != "finite-distance")
            throw std::invalid_argument("--stereo-radiation must be far-field or finite-distance");
        if (!stereoOutput && stereoRadiationModel != "far-field")
            throw std::invalid_argument("--stereo-radiation finite-distance requires --stereo");
        Excitation::ModalDecayModel decayModel;
        if (decayModelName == "fixed-t60") decayModel = Excitation::ModalDecayModel::FixedT60;
        else if (decayModelName == "structural-loss") decayModel = Excitation::ModalDecayModel::StructuralLossFactor;
        else throw std::invalid_argument("--decay-model must be fixed-t60 or structural-loss");
        const auto* materialPreset = find_material_preset(materialPresetName);
        if (!materialPreset) {
            std::string available;
            for (const auto& id : material_preset_ids()) available += " " + id;
            throw std::invalid_argument("unknown --material-preset '" + materialPresetName + "'; available:" + available);
        }
        if (!lossFactorOverride && materialPreset->has_structural_loss_factor)
            structuralLossFactor = materialPreset->structural_loss_factor;
        if (decayModel == Excitation::ModalDecayModel::StructuralLossFactor &&
            (!std::isfinite(structuralLossFactor) || !(structuralLossFactor > 0.0 && structuralLossFactor < 0.1)))
            throw std::invalid_argument("selected material preset has no measured structural loss factor; use --loss-factor only for an explicit sensitivity study");

        const BellGeometry geometry = loadGeometry(geometryFile);
        if (geometry.ellipticity < 0.999999)
            throw std::runtime_error("synth-clapper currently requires circular geometry; elliptical impact needs a 3D clapper path");
        BellMaterial material = materialPreset->material;
        BellMeshConfig meshConfig; meshConfig.nNodes = meshNodes;
        BellGenerator generator(material, meshConfig);
        generator.setGeometry(geometry);
        GeometryValidationResult validation;
        if (!generator.buildMesh(&validation)) throw std::runtime_error("geometry rejected: " + validation.reason);
        const auto& outerProfile = generator.getOuterProfile();
        const auto& innerProfileMm = generator.getInnerProfile();
        if (innerProfileMm.size() < 3 || outerProfile.size() < 3) throw std::runtime_error("bell boundary profiles are incomplete");

        BellMeshConfig directMeshConfig; directMeshConfig.nNodes = meshNodes;
        const auto directModel = FEM::make_cylindrical_bell_model(geometry, material, directMeshConfig);
        Excitation::MeridionalInfluence influence;
        FEM::CylindricalModalEvaluationConfig femConfig;
        femConfig.harmonics = {0, 1, 2, 3};
        femConfig.angular_samples = 64;
        femConfig.requested_modes = 24;
        femConfig.max_iterations = 5000;
        femConfig.tolerance = 1.0e-6;
        femConfig.minimum_frequency_hz = 1.0;
        ModalData modalData;
        if (!FEM::evaluate_cylindrical_bell_modes(directModel, influence.make(directModel.mesh_si), femConfig, modalData))
            throw std::runtime_error("Fourier cylindrical FEM modal solve failed");

        std::vector<const ModalMode*> usedModes;
        for (const auto& m : modalData.modes) {
            const std::size_t basis = m.inner_surface_harmonics.size();
            if (basis == 0 || m.inner_surface_fourier_radial.size() != innerProfileMm.size() * basis ||
                m.inner_surface_fourier_axial.size() != innerProfileMm.size() * basis ||
                m.outer_surface_harmonics.empty()) continue;
            usedModes.push_back(&m);
        }
        if (usedModes.size() < 3) throw std::runtime_error("insufficient modes with both inner and outer boundary fields");

        std::vector<Point2D> innerProfileSi = directModel.mesh_si.inner_profile;
        if (!std::isfinite(pivotZ)) pivotZ = innerProfileSi.front().z - 0.07;
        const double effectiveModulus = 1.0 / ((1.0 - tipPoisson * tipPoisson) / tipYoung +
            (1.0 - material.nu * material.nu) / directModel.young_modulus_pa);
        if (!(hertzStiffness > 0.0)) hertzStiffness = (4.0 / 3.0) * effectiveModulus * std::sqrt(tipRadius);
        if (!std::isfinite(hertzStiffness) || hertzStiffness <= 0.0)
            throw std::invalid_argument("Hertz stiffness must be positive");
        if (!(inertia > 0.0)) inertia = clapperMass * comLength * comLength;

        Excitation::ClapperDynamicsConfig clapper;
        clapper.mass_kg = clapperMass;
        clapper.moment_of_inertia_kg_m2 = inertia;
        clapper.centre_of_mass_length_m = comLength;
        clapper.tip_length_m = tipLength;
        clapper.tip_radius_m = tipRadius;
        clapper.pivot_r_m = pivotR;
        clapper.pivot_z_m = pivotZ;
        clapper.initial_angle_rad = initialAngleDeg * kPi / 180.0;
        clapper.initial_angular_velocity_rad_s = initialAngularVelocity;
        clapper.hertz_stiffness_n_per_m_1_5 = hertzStiffness;
        clapper.strike_azimuth_rad = strikeAzimuthDeg * kPi / 180.0;
        clapper.integration_step_seconds = contactStep;
        clapper.maximum_duration_seconds = impactHorizon;
        Excitation::ClapperDynamicsResult impact;
        std::vector<ModalMode> contactModes;
        contactModes.reserve(usedModes.size());
        for (const ModalMode* m : usedModes) contactModes.push_back(*m);
        if (!Excitation::simulateInternalClapperContact(clapper, innerProfileSi, contactModes, impact))
            throw std::runtime_error("clapper did not complete a contact-and-separation event; adjust pivot, tip path, initial angular velocity, or impact horizon");
        for (std::size_t i = 0; i < impact.modal_displacement_at_separation.size(); ++i) {
            if (!std::isfinite(impact.modal_displacement_at_separation[i]) ||
                !std::isfinite(impact.modal_velocity_at_separation[i]))
                throw std::runtime_error("non-finite modal release state at mode " + std::to_string(i));
        }

        struct ModeAuditRow {
            double frequency_hz;
            double modal_mass;
            double q_release;
            double v_release;
            double alpha;
            double t60_mode;
            std::complex<double> radiation;
        };
        std::vector<RenderMode> renderModes;
        std::vector<std::complex<double>> radiationLeft, radiationRight;
        std::vector<ModeAuditRow> modeAuditRows;
        const double polar = listenerPolarDeg * kPi / 180.0;
        const double azimuth = listenerAzimuthDeg * kPi / 180.0;
        const auto listener = Acoustics::makeStereoListenerGeometry(
            listenerDistance, polar, azimuth, earSpacing,
            headYawDeg * kPi / 180.0);
        double leftDelay = listener.left.propagation_delay_s;
        double rightDelay = listener.right.propagation_delay_s;
        Acoustics::PreparedStereoRadiation preparedStereoRadiation;
        std::vector<Acoustics::StereoRadiationTransfer> preparedStereoTransfers;
        if (stereoOutput && stereoRadiationModel != "finite-distance") {
            std::vector<const ModalMode*> radiationModes(usedModes.begin(), usedModes.end());
            preparedStereoRadiation = Acoustics::prepareStereoRadiation(
                radiationModes, outerProfile, listenerDistance, polar, azimuth, earSpacing,
                0.5 * kPi, 2.0 * kPi / 180.0);
            preparedStereoTransfers = preparedStereoRadiation.query(headYawDeg * kPi / 180.0);
            if (preparedStereoTransfers.size() != usedModes.size())
                throw std::runtime_error("prepared stereo radiation mode count mismatch");
            if (!preparedStereoTransfers.empty()) {
                leftDelay = preparedStereoTransfers.front().left_delay_s;
                rightDelay = preparedStereoTransfers.front().right_delay_s;
            }
        }
        for (std::size_t i = 0; i < usedModes.size(); ++i) {
            Excitation::ModalDecayParameters decay;
            if (!Excitation::calculateModalDecay(usedModes[i]->frequency_hz, decayModel, t60,
                                                  structuralLossFactor, decay))
                throw std::runtime_error("invalid modal decay parameters for mode " + std::to_string(i));
            const double alpha = decay.alpha_per_second;
            const double omega = 2.0 * kPi * usedModes[i]->frequency_hz;
            std::complex<double> centreRadiation{}, radL{}, radR{};
            if (stereoOutput && stereoRadiationModel == "finite-distance") {
                // The structural/contact/modal history remains shared. Only the
                // receiver transfer is evaluated independently at the two exact
                // ear points using the finite-distance Green-kernel proxy.
                radL = Acoustics::evaluateFiniteDistanceRadiation(
                    *usedModes[i], outerProfile, listener.left, 720, 343.0, 16);
                radR = Acoustics::evaluateFiniteDistanceRadiation(
                    *usedModes[i], outerProfile, listener.right, 720, 343.0, 16);
                centreRadiation = 0.5 * (radL + radR);
            } else if (!stereoOutput) {
                // Mono needs only its declared listener direction; retain the
                // single-polar-field path instead of paying for a stereo grid.
                const auto field = Acoustics::buildDirectionalRadiationField(
                    *usedModes[i], outerProfile, polar, listenerDistance);
                centreRadiation = field.evaluate(azimuth);
                radL = radR = centreRadiation;
            } else {
                // The prepared immutable polar family is built once before this
                // loop and can serve repeated listener-yaw queries without FEM,
                // contact, modal release, or surface-integration work.
                const auto& transfer = preparedStereoTransfers.at(i);
                radL = transfer.left;
                radR = transfer.right;
                centreRadiation = 0.5 * (radL + radR);
            }
            if (omega * omega <= alpha * alpha ||
                (std::abs(radL) < 1.0e-18 && std::abs(radR) < 1.0e-18)) continue;
            const double qRelease = impact.modal_displacement_at_separation[i];
            const double vRelease = impact.modal_velocity_at_separation[i];
            renderModes.push_back({omega, alpha, qRelease, vRelease, stereoOutput ? 0.5 * (radL + radR) : centreRadiation});
            const double leftRangeScale = (!stereoOutput || stereoRadiationModel == "finite-distance") ? 1.0 :
                preparedStereoTransfers[i].left_range_scale;
            const double rightRangeScale = (!stereoOutput || stereoRadiationModel == "finite-distance") ? 1.0 :
                preparedStereoTransfers[i].right_range_scale;
            radiationLeft.push_back(radL * leftRangeScale);
            radiationRight.push_back(radR * rightRangeScale);
            modeAuditRows.push_back({usedModes[i]->frequency_hz, usedModes[i]->modal_mass,
                                     qRelease, vRelease, alpha, decay.t60_seconds, stereoOutput ? 0.5 * (radL + radR) : centreRadiation});
        }
        if (renderModes.size() < 3) throw std::runtime_error("fewer than three underdamped radiating modes remain");

        std::vector<Audio::ModalPressureTerm> pressureLeftModes, pressureRightModes, pressureMonoModes;
        pressureLeftModes.reserve(renderModes.size()); pressureRightModes.reserve(renderModes.size());
        pressureMonoModes.reserve(renderModes.size());
        for (std::size_t i = 0; i < renderModes.size(); ++i) {
            const auto& m = renderModes[i];
            pressureLeftModes.push_back({m.omega, m.alpha, m.q0, m.v0, radiationLeft.at(i)});
            pressureRightModes.push_back({m.omega, m.alpha, m.q0, m.v0, radiationRight.at(i)});
            pressureMonoModes.push_back({m.omega, m.alpha, m.q0, m.v0, m.radiation});
        }

        const std::size_t count = static_cast<std::size_t>(std::llround(audioDuration * sampleRate));
        std::vector<double> audio, audioLeft, audioRight;
        const double releaseTime = impact.separation_time_seconds;
        double peak = 0.0;
        if (stereoOutput) {
            audioLeft.assign(count, 0.0); audioRight.assign(count, 0.0);
            for (std::size_t n = 0; n < count; ++n) {
                const double t = static_cast<double>(n) / sampleRate;
                const double pL = Audio::evaluateModalPressureAt(
                    pressureLeftModes, t, releaseTime, leftDelay);
                const double pR = Audio::evaluateModalPressureAt(
                    pressureRightModes, t, releaseTime, rightDelay);
                if (!std::isfinite(pL) || !std::isfinite(pR))
                    throw std::runtime_error("non-finite stereo acoustic sample at time " + std::to_string(t));
                audioLeft[n] = pL; audioRight[n] = pR;
                peak = std::max({peak, std::abs(pL), std::abs(pR)});
            }
        } else {
            audio.assign(count, 0.0);
            for (std::size_t n = 0; n < count; ++n) {
                const double p = Audio::evaluateModalPressureAt(
                    pressureMonoModes, static_cast<double>(n) / sampleRate, releaseTime, 0.0);
                if (!std::isfinite(p)) throw std::runtime_error("non-finite acoustic sample at time " + std::to_string(static_cast<double>(n) / sampleRate));
                audio[n] = p; peak = std::max(peak, std::abs(p));
            }
        }
        if (!(peak > 1.0e-30) || !std::isfinite(peak)) throw std::runtime_error("acoustic radiation proxy is numerically zero");
        const double scale = 0.95 / peak;
        if (stereoOutput) {
            for (double& x : audioLeft) x *= scale;
            for (double& x : audioRight) x *= scale;
        } else {
            for (double& x : audio) x *= scale;
        }
        double sumSquares = 0.0;
        if (stereoOutput) {
            for (std::size_t i = 0; i < count; ++i) sumSquares += 0.5 * (audioLeft[i]*audioLeft[i] + audioRight[i]*audioRight[i]);
            if (!Audio::writeStereoPcm16(output, audioLeft, audioRight, sampleRate)) throw std::runtime_error("failed to write stereo WAV file");
        } else {
            for (double x : audio) sumSquares += x * x;
            if (!Audio::writeMonoPcm16(output, audio, sampleRate)) throw std::runtime_error("failed to write mono WAV file");
        }
        const double rms = std::sqrt(sumSquares / static_cast<double>(count));
        const double crestFactor = rms > 0.0 ? 0.95 / rms : 0.0;

        const double topZ = innerProfileSi.front().z;
        const double mouthZ = innerProfileSi.back().z;
        const double axialSpan = mouthZ - topZ;
        if (std::abs(axialSpan) < 1.0e-9) throw std::runtime_error("inner profile has zero axial span");
        const double strikeFraction = (impact.first_contact_surface_z_m - topZ) / axialSpan;
        const bool strikeInTypicalBand = strikeFraction >= 0.50 && strikeFraction <= 0.70;

        std::filesystem::path modesPath(output);
        modesPath.replace_extension(); modesPath += "_modes.csv";
        std::ofstream modesFile(modesPath);
        if (modesFile) {
            modesFile << "frequency_hz,modal_mass_kg,q_at_separation_m,v_at_separation_m_per_s,decay_alpha_per_s,mode_t60_s,radiation_real,radiation_imag,radiation_magnitude\n";
            modesFile << std::setprecision(17);
            for (const auto& row : modeAuditRows)
                modesFile << row.frequency_hz << ',' << row.modal_mass << ',' << row.q_release << ','
                          << row.v_release << ',' << row.alpha << ',' << row.t60_mode << ','
                          << row.radiation.real() << ',' << row.radiation.imag()
                          << ',' << std::abs(row.radiation) << '\n';
        }

        std::filesystem::path metadataPath(output);
        metadataPath.replace_extension(); metadataPath += "_metadata.txt";
        std::ofstream meta(metadataPath);
        if (meta) {
            meta << "BellOptimizer internal clapper physical-pendulum synthesis\n"
                 << "Geometry: " << geometryFile << "\nMesh nodes: " << meshNodes << "\nModes used: " << renderModes.size() << "\n"
                 << "Clapper mass kg: " << clapperMass << "\nPivot r,z m: " << pivotR << ", " << pivotZ << "\n"
                 << "COM length m: " << comLength << "\nTip centre length m: " << tipLength << "\nTip radius m: " << tipRadius << "\n"
                 << "Pivot inertia kg m2: " << inertia << "\nInitial angle deg: " << initialAngleDeg
                 << "\nInitial angular velocity rad/s: " << initialAngularVelocity << "\n"
                 << "Tip steel E Pa: " << tipYoung << "\nTip Poisson ratio: " << tipPoisson << "\n"
                 << "Effective Hertz modulus Pa: " << effectiveModulus << "\nHertz stiffness N/m^1.5: " << hertzStiffness
                 << "\nContact law: spherical Hertz, locally flat wall approximation; no contact damping/plasticity\n"
                 << "Strike azimuth deg: " << strikeAzimuthDeg << "\nFirst contact s: " << impact.first_contact_time_seconds
                 << "\nSeparation s: " << impact.separation_time_seconds << "\nContact duration s: " << impact.contact_duration_seconds
                 << "\nPeak contact force N: " << impact.peak_contact_force_n
                 << "\nFirst contact surface point r,z m: " << impact.first_contact_surface_r_m << ", " << impact.first_contact_surface_z_m
                 << "\nStrike height fraction (0=top, 1=mouth): " << strikeFraction
                 << "\nStrike height in requested typical band 0.50..0.70: " << (strikeInTypicalBand ? "yes" : "no")
                 << "\nConservative energy relative drift: " << impact.maximum_relative_energy_error
                 << "\nMaterial preset: " << materialPreset->id
                 << "\nMaterial preset description: " << materialPreset->description
                 << "\nMaterial damping source: " << materialPreset->damping_source
                 << "\nDecay model: " << decayModelName
                 << "\nFixed T60 s: " << t60
                 << "\nStructural loss factor eta: " << (std::isfinite(structuralLossFactor) ? std::to_string(structuralLossFactor) : "not used")
                 << "\nLoss factor source: " << (lossFactorOverride ? "explicit sensitivity override" : materialPreset->damping_source)
                 << "\nPer-mode decay alpha and T60 are recorded in modal audit CSV"
                 << "\nListener polar/azimuth deg: " << listenerPolarDeg << ", " << listenerAzimuthDeg
                 << "\nOutput channels: " << (stereoOutput ? "2 (stereo)" : "1 (mono)")
                 << "\nListener distance m: " << listenerDistance << "\nEar spacing m: " << earSpacing
                 << "\nHead yaw deg: " << headYawDeg
                 << "\nEar azimuth deg: " << listener.left.azimuth_rad * 180.0 / kPi << ", " << listener.right.azimuth_rad * 180.0 / kPi
                 << "\nEar polar deg (geometric): " << listener.left.polar_rad * 180.0 / kPi << ", " << listener.right.polar_rad * 180.0 / kPi
                 << "\nEar angular separation deg: " << listener.ear_angular_separation_rad * 180.0 / kPi
                 << "\nEar propagation offsets s (relative to centre): " << leftDelay << ", " << rightDelay
                 << "\nSample rate Hz: " << sampleRate << "\nOutput peak normalization: 0.95 (not calibrated SPL)"
                 << "\nOutput RMS: " << rms << "\nOutput crest factor: " << crestFactor
                 << "\nModal audit CSV: " << modesPath.string() << "\n"
                 << "Radiation model: " << (stereoOutput ? stereoRadiationModel : "far-field") << "\n"
                 << "Radiation: " << (stereoOutput && stereoRadiationModel == "finite-distance" ?
                    "outer-shell finite-distance free-space Green-kernel proxy; not a complete exterior acoustic solution" :
                    "outer-surface directional far-field modal approximation") << "; no room acoustics\n"
                 << "Stereo model: shared FEM/contact/modal response; left/right receiver transfers evaluated at geometric ear points; path-distance propagation delay; no per-channel structural solve\n"
                 << "Limitations: no head shadow, pinna/HRTF, binaural filtering or room reflections\n";
        }
        std::cout << std::fixed << std::setprecision(5)
                  << "Generated " << output << " | channels=" << (stereoOutput ? 2 : 1) << " | modes=" << renderModes.size()
                  << " | first contact=" << impact.first_contact_time_seconds << " s"
                  << " | contact=" << impact.contact_duration_seconds << " s"
                  << " | peak force=" << impact.peak_contact_force_n << " N"
                  << " | strike height=" << std::fixed << std::setprecision(1) << (100.0 * strikeFraction) << "%"
                  << " | RMS=" << std::scientific << rms
                  << " | energy drift=" << impact.maximum_relative_energy_error << "\n";
        if (!strikeInTypicalBand)
            std::cout << "WARNING: first contact is outside the typical 50-70% bell-height band; adjust pivot/trajectory if targeting a conventional strike.\n";
        std::cout << "Mode diagnostics: " << modesPath.string() << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "synth-clapper error: " << e.what() << "\n";
        return 1;
    }
}
