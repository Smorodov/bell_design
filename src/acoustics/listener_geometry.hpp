#pragma once

namespace BellOptimizer::Acoustics {

struct ReceiverPoint {
    double x_m = 0.0;
    double y_m = 0.0;
    double z_m = 0.0;
    double distance_m = 0.0;
    double polar_rad = 0.0;
    double azimuth_rad = 0.0;
    double propagation_delay_s = 0.0;
};

struct StereoListenerGeometry {
    ReceiverPoint centre;
    ReceiverPoint left;
    ReceiverPoint right;
    double ear_angular_separation_rad = 0.0;
};

// Ideal point-ear geometry. The head baseline is tangent to the horizontal
// line of sight at zero yaw; positive yaw turns it toward the radial direction.
// The bell is at the origin and +Z is its axis.
StereoListenerGeometry makeStereoListenerGeometry(double centre_distance_m,
                                                   double centre_polar_rad,
                                                   double centre_azimuth_rad,
                                                   double ear_spacing_m,
                                                   double head_yaw_rad,
                                                   double sound_speed_m_s = 343.0);

} // namespace BellOptimizer::Acoustics
