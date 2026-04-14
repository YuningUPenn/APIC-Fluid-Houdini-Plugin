#include "ApicParticleSet.h"
#include <cassert>

namespace apic {

    void ApicParticleSet::reserve(size_t n) {
        positions_.reserve(n);
        velocities_.reserve(n);
        masses_.reserve(n);
        affineBs_.reserve(n);
    }

    void ApicParticleSet::clear() {
        positions_.clear();
        velocities_.clear();
        masses_.clear();
        affineBs_.clear();
    }

    void ApicParticleSet::addParticle(const Vec3& pos, const Vec3& vel, Scalar mass) {
        positions_.push_back(pos);
        velocities_.push_back(vel);
        masses_.push_back(mass);
        affineBs_.push_back(Mat3::Zero());
    }

    void ApicParticleSet::zeroAffineBs() {
        for (auto& B : affineBs_) {
            B.setZero();
        }
    }

    // -------------------------------------------------------
    // Angular momentum computation
    //
    // From Jiang et al. Eq (11):
    //   L_tot = sum_p  m_p * (x_p x v_p)  +  (B_p)^T : epsilon
    //
    // The second term is the angular momentum encoded in B_p.
    // For a skew-symmetric part of B_p:
    //   (B^T : epsilon)_k = sum_{alpha,beta} B_{beta,alpha} * epsilon_{alpha,beta,k}
    //
    // Simpler form: the angular momentum from B is
    //   L_from_B = [B_{32}-B_{23}, B_{13}-B_{31}, B_{21}-B_{12}] * m_p
    // (the axial vector of the skew part of B_p, times m_p)
    // -------------------------------------------------------
    Vec3 ApicParticleSet::totalAngularMomentum() const {
        Vec3 L = Vec3::Zero();

        for (size_t p = 0; p < positions_.size(); ++p) {
            // Orbital angular momentum
            L += masses_[p] * positions_[p].cross(velocities_[p]);

            // Angular momentum stored in affine matrix (skew-symmetric part)
            const Mat3& B = affineBs_[p];
            Vec3 Lspin;

            Lspin.x() = B(2, 1) - B(1, 2);
            Lspin.y() = B(0, 2) - B(2, 0);
            Lspin.z() = B(1, 0) - B(0, 1);

            L += masses_[p] * Lspin;
        }

        return L;
    }

} // namespace apic