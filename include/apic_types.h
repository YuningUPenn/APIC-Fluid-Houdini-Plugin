#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <cstdint>

namespace apic {

// -------------------------------------------------------
// Basic types
// -------------------------------------------------------
using Scalar  = float;
using Vec3    = Eigen::Matrix<Scalar, 3, 1>;
using Mat3    = Eigen::Matrix<Scalar, 3, 3>;
using Vec3i   = Eigen::Matrix<int, 3, 1>;

// Convenience constructors (real Eigen supports Vec3(x,y,z) directly;
// this helper also works with minimal stubs and unit tests)
inline Vec3 vec3(Scalar x, Scalar y, Scalar z) {
    Vec3 v; v.x() = x; v.y() = y; v.z() = z; return v;
}
inline Vec3i vec3i(int x, int y, int z) {
    Vec3i v; v(0) = x; v(1) = y; v(2) = z; return v;
}

// -------------------------------------------------------
// Transfer method selector
// -------------------------------------------------------
enum class TransferMethod : int {
    APIC   = 0,
    FLIP   = 1,
    PIC    = 2,
    Hybrid = 3   // alpha-blend of APIC and FLIP
};

// -------------------------------------------------------
// Global simulation parameters passed through the HDA
// -------------------------------------------------------
struct SimParams {
    // Grid
    int   gridResX     = 64;
    int   gridResY     = 64;
    int   gridResZ     = 64;
    Scalar gridSpacing = 0.05f;   // world-space cell width (dx)

    // Timing
    Scalar dt          = 1.0f / 24.0f;
    int    substeps    = 2;
    Scalar timeScale   = 1.0f;

    // Physics
    Vec3   gravity     = Vec3::Zero();  // set to (0,-9.8,0) in solver init
    Scalar viscosity   = 0.01f;
    Scalar surfaceTension = 0.0f;
    Scalar density     = 1000.f;  // kg/m^3 (water)

    // Transfer
    TransferMethod method  = TransferMethod::APIC;
    Scalar flipAlpha       = 0.95f;  // used when method == Hybrid

    // Boundary
    Scalar boundaryFriction = 0.0f;

    // Particle init
    Scalar particleSeparation = 0.05f;
};

} // namespace apic
