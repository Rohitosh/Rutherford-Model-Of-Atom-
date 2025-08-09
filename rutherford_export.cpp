// rutherford_export.cpp
// Compile: g++ -O2 -std=c++17 rutherford_export.cpp -o rutherford_export
// Run: ./rutherford_export
// Produces: trajectories.csv and angles.csv

#include <bits/stdc++.h>
using namespace std;

// Physical constants
const double e_charge = 1.602176634e-19;     // C
const double epsilon0  = 8.8541878128e-12;   // F/m
const double pi = acos(-1.0);
const double k_coulomb = 1.0 / (4.0 * pi * epsilon0);

// alpha mass (kg)
const double amu = 1.66053906660e-27;
const double mass_alpha = 4.0 * amu;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // --- Parameters (tweak if you like) ---
    int N = 2500;              // number of particles to simulate
    double E_MeV = 5.0;        // alpha energy (MeV)
    double Z1 = 2.0;           // alpha charge number
    double Z2 = 79.0;          // gold nucleus charge number
    double bmax = 1.0e-10;     // maximum impact parameter (m) ~ atomic scale
    int framesBefore = 220;    // frames before hitting foil (pre-foil)
    int framesAfter  = 400;    // frames after scattering
    double startX = -6e-14;    // starting x position (m) left of foil
    double foilX  = 0.0;       // foil plane at x=0
    double exitX  = +6e-14;    // stop recording when x>exitX
    double frame_dt = 1e-17;   // artificial time per frame (seconds) for animation smoothing

    // Derived
    double E_joule = E_MeV * 1.0e6 * e_charge;
    // approximate speed non-relativistic
    double v0 = sqrt(2.0 * E_joule / mass_alpha);

    // Output files
    ofstream traj("trajectories.csv");
    ofstream angs("angles.csv");
    if (!traj || !angs) {
        cerr << "Cannot open output files.\n";
        return 1;
    }

    // CSV headers
    // trajectories.csv: particle,frame,x_m,y_m
    traj << "particle,frame,x_m,y_m\n";
    // angles.csv: particle,theta_deg
    angs << "particle,theta_deg\n";

    // RNG
    mt19937_64 rng((unsigned)chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> unif01(0.0, 1.0);
    uniform_real_distribution<double> unifSign(-1.0, 1.0);

    for (int p = 0; p < N; ++p) {
        // area-weighted impact parameter: b = bmax * sqrt(u)
        double u = unif01(rng);
        double b = bmax * sqrt(u);
        // random sign for y offset
        double sign = (unifSign(rng) >= 0.0 ? 1.0 : -1.0);
        double y0 = sign * b;

        // scattering angle from Rutherford formula
        // tan(theta/2) = k*q1*q2 / (2 E b)
        double q1 = Z1 * e_charge;
        double q2 = Z2 * e_charge;
        double t2 = (k_coulomb * q1 * q2) / (2.0 * E_joule * b);
        double theta; // radians
        if (t2 > 1e300) theta = pi; else theta = 2.0 * atan(t2);
        double theta_deg = theta * 180.0 / pi;
        angs << p << "," << fixed << setprecision(8) << theta_deg << "\n";

        // Build frames before hitting foil: particles move with v0 along +x
        // pre-foil: starting at startX - we will produce framesBefore frames with equal spacing until x reaches foilX
        double pre_dx = (foilX - startX) / double(framesBefore);
        for (int f = 0; f < framesBefore; ++f) {
            double x = startX + pre_dx * f;
            double y = y0;
            traj << p << "," << f << "," << scientific << x << "," << scientific << y << "\n";
        }

        // at impact point x=foilX, the particle instantly has scattering angle theta (classical Rutherford).
        // post-foil velocity components: still speed ~ v0, direction deflected by ±theta around +x.
        // For sign of theta we use original y sign: if y0>0 deflect upward (positive)
        double theta_signed = (y0 >= 0.0 ? theta : -theta);
        double vx_after = v0 * cos(theta_signed);
        double vy_after = v0 * sin(theta_signed);

        // Now produce post-frames: we step in time using frame_dt
        double x = foilX;
        double y = y0;
        for (int f2 = 0; f2 < framesAfter; ++f2) {
            double frameIndex = framesBefore + f2;
            // update position using simple straight-line motion (we already accounted scattering angle)
            x += vx_after * frame_dt;
            y += vy_after * frame_dt;
            traj << p << "," << int(frameIndex) << "," << scientific << x << "," << scientific << y << "\n";
            if (x > exitX) break;
            if (fabs(y) > 1e-11) break; // escape region
        }
    }

    traj.close();
    angs.close();
    cout << "Wrote trajectories.csv and angles.csv  (particles: " << N << ")\n";
    cout << "Parameters: E(MeV)=" << E_MeV << "  bmax(m)=" << bmax << "  v0(m/s)=" << v0 << "\n";
    return 0;
}
