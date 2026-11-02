#include "synthesize_contact.hpp"
#include "../core/bell_types.h"
#include "../core/bell_generator.h"
#include "../core/modal_types.hpp"
#include "../core/simple_json.hpp"
#include "../fem/cylindrical_bell_model.hpp"
#include "../fem/cylindrical_modal_evaluator.hpp"
#include "../excitation/hammer_contact.hpp"
#include "../excitation/clapper_dynamics.hpp"
#include "../excitation/meridional_influence.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
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

void writeU16(std::ofstream& f, std::uint16_t v) {
    const char b[2] = {static_cast<char>(v & 255), static_cast<char>((v >> 8) & 255)};
    f.write(b, 2);
}
void writeU32(std::ofstream& f, std::uint32_t v) {
    const char b[4] = {static_cast<char>(v & 255), static_cast<char>((v >> 8) & 255),
                       static_cast<char>((v >> 16) & 255), static_cast<char>((v >> 24) & 255)};
    f.write(b, 4);
}
bool writeMonoPcm16(const std::string& filename, const std::vector<double>& samples, int sampleRate) {
    if (samples.empty() || sampleRate <= 0 || samples.size() > (std::numeric_limits<std::uint32_t>::max() - 36u) / 2u)
        return false;
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;
    const auto bytes = static_cast<std::uint32_t>(samples.size() * 2u);
    f.write("RIFF", 4); writeU32(f, 36u + bytes); f.write("WAVE", 4);
    f.write("fmt ", 4); writeU32(f, 16); writeU16(f, 1); writeU16(f, 1);
    writeU32(f, static_cast<std::uint32_t>(sampleRate)); writeU32(f, static_cast<std::uint32_t>(sampleRate * 2));
    writeU16(f, 2); writeU16(f, 16); f.write("data", 4); writeU32(f, bytes);
    for (double x : samples) {
        if (!std::isfinite(x)) return false;
        const double y = std::clamp(x, -1.0, 1.0);
        const auto v = static_cast<std::int16_t>(std::lrint(y * (y < 0.0 ? 32768.0 : 32767.0)));
        writeU16(f, static_cast<std::uint16_t>(v));
    }
    return static_cast<bool>(f);
}

struct RenderMode { double omega, alpha, q0, v0; std::complex<double> radiation; };

std::size_t basisForLabel(const std::vector<int>& labels, int label) {
    const auto it = std::find(labels.begin(), labels.end(), label);
    return it == labels.end() ? labels.size() : static_cast<std::size_t>(it - labels.begin());
}

double fourierValue(const std::vector<double>& values, std::size_t point,
                    std::size_t pointCount, const std::vector<int>& labels, double angle) {
    if (values.size() != pointCount * labels.size()) return 0.0;
    double value = 0.0;
    for (std::size_t b = 0; b < labels.size(); ++b) {
        const int label = labels[b];
        const double c = values[point * labels.size() + b];
        if (label == 0) value += c;
        else if (label > 0) value += c * std::cos(label * angle);
        else value += c * std::sin((-label) * angle);
    }
    return value;
}

std::complex<double> directionalRadiation(const ModalMode& mode,
        const std::vector<Point2D>& profile, double listenerPolar,
        double listenerAzimuth, unsigned azimuthSamples = 96) {
    if (profile.size() < 2 || mode.outer_surface_harmonics.empty() ||
        mode.outer_surface_fourier_radial.empty() || mode.outer_surface_fourier_axial.empty()) return {};
    constexpr double soundSpeed = 343.0;
    const double omega = 2.0 * kPi * mode.frequency_hz;
    const double waveNumber = omega / soundSpeed;
    const double sinEl = std::sin(listenerPolar), cosEl = std::cos(listenerPolar);
    std::complex<double> integral{};
    for (std::size_t j = 0; j + 1 < profile.size(); ++j) {
        const double dr = profile[j + 1].r - profile[j].r;
        const double dz = profile[j + 1].z - profile[j].z;
        const double ds = std::hypot(dr, dz);
        if (!(ds > 0.0)) continue;
        for (unsigned a = 0; a < azimuthSamples; ++a) {
            const double phi = 2.0 * kPi * (static_cast<double>(a) + 0.5) / azimuthSamples;
            const double nr = dz / ds, nz = -dr / ds;
            const double ur0 = 0.5 * (fourierValue(mode.outer_surface_fourier_radial, j,
                profile.size(), mode.outer_surface_harmonics, phi) +
                fourierValue(mode.outer_surface_fourier_radial, j + 1,
                profile.size(), mode.outer_surface_harmonics, phi));
            const double uz0 = 0.5 * (fourierValue(mode.outer_surface_fourier_axial, j,
                profile.size(), mode.outer_surface_harmonics, phi) +
                fourierValue(mode.outer_surface_fourier_axial, j + 1,
                profile.size(), mode.outer_surface_harmonics, phi));
            const double normal = ur0 * nr + uz0 * nz;
            const double r = 0.5 * (profile[j].r + profile[j + 1].r) * 1.0e-3;
            const double z = 0.5 * (profile[j].z + profile[j + 1].z) * 1.0e-3;
            const double dsMeters = ds * 1.0e-3;
            const double phase = -waveNumber * (r * sinEl * std::cos(phi - listenerAzimuth) + z * cosEl);
            integral += normal * r * dsMeters * (2.0 * kPi / azimuthSamples) *
                        std::complex<double>(std::cos(phase), std::sin(phase));
        }
    }
    return integral;
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
            radiationCoefficients.push_back(directionalRadiation(m, profile, listenerPolar, listenerAzimuth));
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
        if (!writeMonoPcm16(output, audio, sampleRate)) throw std::runtime_error("failed to write WAV file");

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
            std::cout << "Usage: bell synth-clapper [geometry.json] [--output file.wav] [--clapper-mass kg] "
                         "[--clapper-inertia kg-m2] [--com-length m] [--tip-length m] [--tip-radius m] "
                         "[--pivot-z m] [--initial-angle-deg deg] [--initial-angular-velocity rad/s] "
                         "[--contact-step s] [--impact-horizon s] [--t60 s] [--duration s] "
                         "[--strike-azimuth-deg deg] [--listener-polar-deg deg] [--listener-azimuth-deg deg]\n"
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
        else if (a == "--t60") t60 = std::stod(value("--t60"));
        else if (a == "--duration") audioDuration = std::stod(value("--duration"));
        else if (a == "--sample-rate") sampleRate = std::stoi(value("--sample-rate"));
        else if (a == "--strike-azimuth-deg") strikeAzimuthDeg = std::stod(value("--strike-azimuth-deg"));
        else if (a == "--listener-polar-deg") listenerPolarDeg = std::stod(value("--listener-polar-deg"));
        else if (a == "--listener-azimuth-deg") listenerAzimuthDeg = std::stod(value("--listener-azimuth-deg"));
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
            !std::isfinite(t60) || t60 <= 0.0 || !std::isfinite(audioDuration) || audioDuration <= 0.1 ||
            sampleRate < 8000 || sampleRate > 192000 || meshNodes < 11 || meshNodes > 101 ||
            audioDuration * sampleRate > 10000000.0 || !std::isfinite(strikeAzimuthDeg) ||
            !std::isfinite(listenerPolarDeg) || listenerPolarDeg < 0.0 || listenerPolarDeg > 180.0 ||
            !std::isfinite(listenerAzimuthDeg))
            throw std::invalid_argument("invalid clapper, contact, audio, or mesh parameter");

        const BellGeometry geometry = loadGeometry(geometryFile);
        if (geometry.ellipticity < 0.999999)
            throw std::runtime_error("synth-clapper currently requires circular geometry; elliptical impact needs a 3D clapper path");
        BellMaterial material;
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
        std::vector<std::complex<double>> radiation;
        for (const auto& m : modalData.modes) {
            const std::size_t basis = m.inner_surface_harmonics.size();
            if (basis == 0 || m.inner_surface_fourier_radial.size() != innerProfileMm.size() * basis ||
                m.inner_surface_fourier_axial.size() != innerProfileMm.size() * basis ||
                m.outer_surface_harmonics.empty()) continue;
            usedModes.push_back(&m);
            radiation.push_back(directionalRadiation(m, outerProfile,
                listenerPolarDeg * kPi / 180.0, listenerAzimuthDeg * kPi / 180.0));
        }
        if (usedModes.size() < 3) throw std::runtime_error("insufficient modes with both inner and outer boundary fields");

        std::vector<Point2D> innerProfileSi = directModel.mesh_si.inner_profile;
        if (!std::isfinite(pivotZ)) pivotZ = innerProfileSi.front().z - 0.03;
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

        std::vector<RenderMode> renderModes;
        for (std::size_t i = 0; i < usedModes.size(); ++i) {
            const double alpha = std::log(1000.0) / t60;
            const double omega = 2.0 * kPi * usedModes[i]->frequency_hz;
            if (omega * omega <= alpha * alpha || std::abs(radiation[i]) < 1.0e-18) continue;
            renderModes.push_back({omega, alpha, impact.modal_displacement_at_separation[i],
                                   impact.modal_velocity_at_separation[i], radiation[i]});
        }
        if (renderModes.size() < 3) throw std::runtime_error("fewer than three underdamped radiating modes remain");

        const std::size_t count = static_cast<std::size_t>(std::llround(audioDuration * sampleRate));
        std::vector<double> audio(count, 0.0);
        const double releaseTime = impact.separation_time_seconds;
        double peak = 0.0;
        for (std::size_t n = 0; n < count; ++n) {
            const double tau = static_cast<double>(n) / sampleRate - releaseTime;
            if (tau < 0.0) continue;
            double p = 0.0;
            for (const auto& m : renderModes) {
                const double wd = std::sqrt(m.omega * m.omega - m.alpha * m.alpha);
                const double b = (m.v0 + m.alpha * m.q0) / wd;
                const std::complex<double> qdd = std::complex<double>(m.q0, -b) *
                    std::exp(std::complex<double>(-m.alpha, wd) * tau) *
                    std::pow(std::complex<double>(-m.alpha, wd), 2);
                p += std::real(m.radiation * qdd);
            }
            if (!std::isfinite(p)) throw std::runtime_error("non-finite acoustic sample at time " + std::to_string(static_cast<double>(n) / sampleRate));
            audio[n] = p;
            peak = std::max(peak, std::abs(p));
        }
        if (!(peak > 1.0e-30) || !std::isfinite(peak)) throw std::runtime_error("acoustic radiation proxy is numerically zero");
        const double scale = 0.95 / peak;
        for (double& x : audio) x *= scale;
        if (!writeMonoPcm16(output, audio, sampleRate)) throw std::runtime_error("failed to write WAV file");

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
                 << "\nConservative energy relative drift: " << impact.maximum_relative_energy_error
                 << "\nT60 s: " << t60 << "\nListener polar/azimuth deg: " << listenerPolarDeg << ", " << listenerAzimuthDeg
                 << "\nSample rate Hz: " << sampleRate << "\nOutput peak normalization: 0.95 (not calibrated SPL)\n"
                 << "Radiation: outer-surface directional far-field modal approximation; no room acoustics\n";
        }
        std::cout << std::fixed << std::setprecision(5)
                  << "Generated " << output << " | modes=" << renderModes.size()
                  << " | first contact=" << impact.first_contact_time_seconds << " s"
                  << " | contact=" << impact.contact_duration_seconds << " s"
                  << " | peak force=" << impact.peak_contact_force_n << " N"
                  << " | energy drift=" << std::scientific << impact.maximum_relative_energy_error << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "synth-clapper error: " << e.what() << "\n";
        return 1;
    }
}
