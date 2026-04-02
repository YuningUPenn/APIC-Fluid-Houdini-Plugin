#include "CollisionHandler.h"
#include <algorithm>
#include <cmath>

namespace apic {

void CollisionHandler::addObstacle(SDFFunc sdf) {
    obstacles_.push_back(std::move(sdf));
}

Scalar CollisionHandler::evalSDF(const Vec3& pos) const {
    Scalar d = std::numeric_limits<Scalar>::max();
    for (const auto& sdf : obstacles_) {
        d = std::min(d, sdf(pos));
    }
    return d;
}

void CollisionHandler::markSolidNodes(ApicGrid& grid) const {
    if (obstacles_.empty()) return;
    const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
    for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
        Vec3 pos = grid.nodePos(i, j, k);
        if (evalSDF(pos) <= 0.f)
            grid.solid(i, j, k) = true;
    }
}

// -------------------------------------------------------
// Enforce grid boundary conditions on MAC faces.
// For the basic domain-wall case, enforceBoundaryFaces()
// already zeros the appropriate faces.
// SDF obstacle support is deferred (future work).
// -------------------------------------------------------
void CollisionHandler::enforceGridBoundary(ApicGrid& grid,
                                            const SimParams& /*params*/) const {
    // Domain walls: handled by grid.enforceBoundaryFaces()
    // which is called from normaliseFaceMomentum() and gridUpdate().
    // SDF obstacles would zero face velocities between solid and fluid nodes here.
    grid.enforceBoundaryFaces();
}

// -------------------------------------------------------
// Resolve particle-solid collisions
// -------------------------------------------------------
void CollisionHandler::resolveParticleCollisions(ApicParticleSet& particles,
                                                   const SimParams& params) const {
    const Scalar dx = params.gridSpacing;
    const size_t np = particles.size();

    for (size_t p = 0; p < np; ++p) {
        Vec3& xp = particles.position(p);
        Vec3& vp = particles.velocity(p);

        for (const auto& sdf : obstacles_) {
            Scalar d = sdf(xp);
            if (d < 0.f) {
                Vec3 n;
                {Vec3 e=Vec3::Zero(); e.x()=dx*0.01f;
                 n.x() = sdf(xp+e) - sdf(xp-e);}
                {Vec3 e=Vec3::Zero(); e.y()=dx*0.01f;
                 n.y() = sdf(xp+e) - sdf(xp-e);}
                {Vec3 e=Vec3::Zero(); e.z()=dx*0.01f;
                 n.z() = sdf(xp+e) - sdf(xp-e);}
                if (n.squaredNorm() > 1e-8f) n.normalize();
                xp -= d * n;
                Scalar vn = vp.dot(n);
                if (vn < 0.f) {
                    vp -= vn * n;
                    vp *= (1.f - params.boundaryFriction);
                }
            }
        }
    }
}

CollisionHandler::SDFFunc CollisionHandler::makeBoxSDF(const Vec3& minPt,
                                                         const Vec3& maxPt) {
    return [minPt, maxPt](const Vec3& p) -> Scalar {
        Vec3 d;
        for (int i = 0; i < 3; ++i) {
            Scalar lo = p(i,0) - minPt(i,0);
            Scalar hi = maxPt(i,0) - p(i,0);
            d(i,0) = std::min(lo, hi);
        }
        return std::min({d.x(), d.y(), d.z()});
    };
}

} // namespace apic
