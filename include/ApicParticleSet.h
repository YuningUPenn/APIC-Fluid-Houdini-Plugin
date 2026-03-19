#pragma once
#include "apic_types.h"

namespace apic {

// -------------------------------------------------------
// ApicParticleSet
//
// Stores all per-particle data required by the APIC method.
// Each particle p carries:
//   x_p  -- position       (Vec3)
//   v_p  -- velocity       (Vec3)
//   m_p  -- mass           (Scalar, constant)
//   B_p  -- affine matrix  (Mat3)  <-- the key APIC state
//
// The affine matrix B encodes the locally-affine velocity
// field around the particle.  After a G2P step:
//   B_p^{n+1} = sum_i  w_ip * v_i^{n+1} * (x_i - x_p)^T
//
// P2G uses B to transfer sub-grid velocity information:
//   m_i * v_i += w_ip * m_p * (v_p + B_p * D_p^{-1} * (x_i - x_p))
//
// For the quadratic interpolation kernel, D_p = (1/4)*dx^2*I,
// so D_p^{-1} = (4/dx^2)*I  (just a scalar factor).
// -------------------------------------------------------
class ApicParticleSet {
public:
    // ---- Lifecycle ----
    ApicParticleSet() = default;
    void reserve(size_t n);
    void clear();
    void addParticle(const Vec3& pos, const Vec3& vel, Scalar mass);

    // ---- Accessors ----
    size_t size() const { return positions_.size(); }

    // Per-particle read/write
    Vec3&       position(size_t i)       { return positions_[i]; }
    const Vec3& position(size_t i) const { return positions_[i]; }

    Vec3&       velocity(size_t i)       { return velocities_[i]; }
    const Vec3& velocity(size_t i) const { return velocities_[i]; }

    Scalar      mass(size_t i)     const { return masses_[i]; }
    Scalar&     mass(size_t i)           { return masses_[i]; }

    Mat3&       affineB(size_t i)        { return affineBs_[i]; }
    const Mat3& affineB(size_t i)  const { return affineBs_[i]; }

    // Bulk SoA access (useful for vectorized loops)
    std::vector<Vec3>&   positions()  { return positions_; }
    std::vector<Vec3>&   velocities() { return velocities_; }
    std::vector<Scalar>& masses()     { return masses_; }
    std::vector<Mat3>&   affineBs()   { return affineBs_; }

    // ---- Angular momentum diagnostic ----
    // Computes total angular momentum: sum_p m_p * (x_p x v_p) + tr(B_p^T * epsilon)
    // The second term captures the angular momentum stored in B.
    Vec3 totalAngularMomentum() const;

    // ---- Reset affine matrices (used for PIC/FLIP fallback) ----
    void zeroAffineBs();

private:
    std::vector<Vec3>   positions_;
    std::vector<Vec3>   velocities_;
    std::vector<Scalar> masses_;
    std::vector<Mat3>   affineBs_;   // APIC affine state
};

} // namespace apic
