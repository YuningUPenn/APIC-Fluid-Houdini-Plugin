#include "ApicGrid.h"
#include <cmath>
#include <algorithm>
#include <cassert>

namespace apic {

// -------------------------------------------------------
// Grid initialisation
// -------------------------------------------------------
void ApicGrid::init(int nx, int ny, int nz, Scalar dx, const Vec3& origin) {
    nx_ = nx; ny_ = ny; nz_ = nz;
    dx_ = dx;
    origin_ = origin;

    size_t n = static_cast<size_t>(nx * ny * nz);
    mass_.assign(n, 0.f);
    velocity_.assign(n, Vec3::Zero());
    force_.assign(n, Vec3::Zero());
    pressure_.assign(n, 0.f);
    divergence_.assign(n, 0.f);
    solid_.assign(n, 0);
    fluid_.assign(n, 0);
    velocityOld_.assign(n, Vec3::Zero());

    markBoundary();
}

void ApicGrid::clear() {
    std::fill(mass_.begin(),       mass_.end(),       0.f);
    std::fill(velocity_.begin(),   velocity_.end(),   Vec3::Zero());
    std::fill(force_.begin(),      force_.end(),      Vec3::Zero());
    std::fill(pressure_.begin(),   pressure_.end(),   0.f);
    std::fill(divergence_.begin(), divergence_.end(), 0.f);
    std::fill(fluid_.begin(),      fluid_.end(),      0);
    // solid_ is NOT cleared -- it's set once and stays
}

void ApicGrid::markBoundary() {
    for (int k = 0; k < nz_; ++k)
    for (int j = 0; j < ny_; ++j)
    for (int i = 0; i < nx_; ++i) {
        bool isBound = (i == 0 || i == nx_-1 ||
                        j == 0 || j == ny_-1 ||
                        k == 0 || k == nz_-1);
        solid(i, j, k) = isBound;
    }
}

Vec3i ApicGrid::idx3(int flat) const {
    int k = flat / (nx_ * ny_);
    int rem = flat % (nx_ * ny_);
    int j = rem / nx_;
    int i = rem % nx_;
    return {i, j, k};
}

Vec3 ApicGrid::nodePos(int flat) const {
    Vec3i v = idx3(flat);
    Vec3 p; p.x() = v[0]; p.y() = v[1]; p.z() = v[2];
    return origin_ + dx_ * p;
}

void ApicGrid::snapshotVelocity() {
    velocityOld_ = velocity_;
}

// -------------------------------------------------------
// Normalise momentum->velocity after P2G accumulation
// -------------------------------------------------------
void ApicGrid::normaliseMomentum() {
    size_t n = numNodes();
    for (size_t i = 0; i < n; ++i) {
        if (mass_[i] > 1e-10f) {
            velocity_[i] *= (1.0f / mass_[i]);
            fluid_[i] = 1;
        } else {
            velocity_[i].setZero();
        }
    }
}

// -------------------------------------------------------
// Quadratic B-spline interpolation kernel
// N(x)  for |x| < 1.5
//
// Same kernel used in Jiang et al. / standard MPM literature:
//
//   N(x) = 0.75 - x^2            for |x| <= 0.5
//   N(x) = 0.5*(1.5 - |x|)^2     for 0.5 < |x| <= 1.5
//   N(x) = 0                      otherwise
// -------------------------------------------------------
Scalar ApicGrid::kernel(Scalar x) {
    Scalar ax = std::abs(x);
    if (ax < 0.5f)  return 0.75f - ax * ax;
    if (ax < 1.5f)  { Scalar t = 1.5f - ax; return 0.5f * t * t; }
    return 0.f;
}

Scalar ApicGrid::kernelGrad(Scalar x) {
    // dN/dx  (not divided by dx -- caller divides by dx to get world-space grad)
    Scalar ax = std::abs(x);
    if (ax < 0.5f)  return -2.f * x;
    if (ax < 1.5f)  return (x > 0.f ? -1.f : 1.f) * (1.5f - ax);
    return 0.f;
}

// -------------------------------------------------------
// gatherWeights
//
// For a particle at world position 'pos', computes the
// 3x3x3 stencil of grid nodes and their weights.
//
// Particle coordinates in grid space:
//   s = (pos - origin) / dx
//
// Base node:
//   base = floor(s - 0.5) + 1   (for quadratic kernel, stencil is
//                                 centred on the nearest node)
//
// For each offset (di, dj, dk) in {-1, 0, 1}:
//   x_frac = s.x - (base.x + di)
//   w_x    = N(x_frac)
//   gw_x   = N'(x_frac) / dx   (world-space gradient component)
// -------------------------------------------------------
void ApicGrid::gatherWeights(const Vec3& pos,
                              std::vector<WeightEntry>& entries) const {
    entries.clear();
    entries.reserve(27);

    // Particle in grid coordinates
    Vec3 s = (pos - origin_) * (1.0f / dx_);

    // Base node (closest node to the left for 3-node stencil)
    int bx = static_cast<int>(std::floor(s.x() - 0.5f)) + 1;
    int by = static_cast<int>(std::floor(s.y() - 0.5f)) + 1;
    int bz = static_cast<int>(std::floor(s.z() - 0.5f)) + 1;

    for (int dk = -1; dk <= 1; ++dk)
    for (int dj = -1; dj <= 1; ++dj)
    for (int di = -1; di <= 1; ++di) {
        int ni = bx + di;
        int nj = by + dj;
        int nk = bz + dk;

        if (!inside(ni, nj, nk)) continue;

        // Fractional distances (in grid units)
        Scalar fx = s.x() - ni;
        Scalar fy = s.y() - nj;
        Scalar fz = s.z() - nk;

        Scalar wx = kernel(fx);
        Scalar wy = kernel(fy);
        Scalar wz = kernel(fz);
        Scalar w  = wx * wy * wz;

        if (w < 1e-12f) continue;

        // World-space gradient: grad_w = (dw/dx, dw/dy, dw/dz)
        //   dw/dx = (dN_x/df_x * 1/dx) * N_y * N_z
        Vec3 grad;
        grad.x() = (kernelGrad(fx) / dx_) * wy * wz;
        grad.y() = wx * (kernelGrad(fy) / dx_) * wz;
        grad.z() = wx * wy * (kernelGrad(fz) / dx_);

        // Offset x_i - x_p in world space
        Vec3 xi = nodePos(ni, nj, nk);
        Vec3 offset = xi - pos;

        entries.push_back({idx(ni, nj, nk), w, grad, offset});
    }
}

} // namespace apic
