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

// -------------------------------------------------------
// Mark solid nodes
// -------------------------------------------------------
void CollisionHandler::markSolidNodes(ApicGrid& grid) const {
    if (obstacles_.empty()) return;

    const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
    for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
        Vec3 pos = grid.nodePos(i, j, k);
        if (evalSDF(pos) <= 0.f) {
            grid.solid(i, j, k) = true;
        }
    }
}

// -------------------------------------------------------
// Enforce grid boundary conditions
// Zero velocity component into solid nodes.
// -------------------------------------------------------
void CollisionHandler::enforceGridBoundary(ApicGrid& grid,
                                            const SimParams& params) const {
    const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
    const Scalar dx = grid.dx();

    for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
        if (!grid.solid(i, j, k)) continue;
        Vec3& v = grid.velocity(i, j, k);

        if (obstacles_.empty()) {
            // Domain walls: zero the velocity components pointing into the wall
            if (i == 0 && v.x() < 0.f)      v.x() = 0.f;
            if (i == nx-1 && v.x() > 0.f)   v.x() = 0.f;
            if (j == 0 && v.y() < 0.f)      v.y() = 0.f;
            if (j == ny-1 && v.y() > 0.f)   v.y() = 0.f;
            if (k == 0 && v.z() < 0.f)      v.z() = 0.f;
            if (k == nz-1 && v.z() > 0.f)   v.z() = 0.f;
        } else {
            // SDF-based: compute normal from SDF gradient
            Vec3 pos = grid.nodePos(i, j, k);
            Vec3 n;
            {Vec3 d=Vec3::Zero(); d.x()=dx;
             n.x() = evalSDF(pos+d) - evalSDF(pos-d);}
            {Vec3 d=Vec3::Zero(); d.y()=dx;
             n.y() = evalSDF(pos+d) - evalSDF(pos-d);}
            {Vec3 d=Vec3::Zero(); d.z()=dx;
             n.z() = evalSDF(pos+d) - evalSDF(pos-d);}
            if (n.squaredNorm() > 1e-8f) n.normalize();

            // Remove normal component (no-penetration)
            Scalar vn = v.dot(n);
            if (vn < 0.f) {
                v -= vn * n;
                // Apply friction to tangential component
                v *= (1.f - params.boundaryFriction);
            }
        }
    }
}

// -------------------------------------------------------
// Resolve particle-solid collisions
// Push particles out and zero normal velocity.
// -------------------------------------------------------
void CollisionHandler::resolveParticleCollisions(ApicParticleSet& particles,
                                                   const SimParams& params) const {
    const Scalar dx = params.gridSpacing;
    const size_t np = particles.size();

    for (size_t p = 0; p < np; ++p) {
        Vec3& xp = particles.position(p);
        Vec3& vp = particles.velocity(p);

        // Domain boundary — use grid origin and size from params
        Vec3 lo, hi;
        lo.x() = lo.y() = lo.z() = -static_cast<float>(params.gridResX) * dx * 0.5f + dx;
        hi.x() = hi.y() = hi.z() = static_cast<float>(params.gridResX) * dx * 0.5f - dx;

        for (int axis = 0; axis < 3; ++axis) {
            Scalar pv = xp(axis);
            Scalar lv = lo(axis);
            Scalar hv = hi(axis);
            if (pv < lv) { xp(axis) = lv; if (vp(axis) < 0.f) vp(axis) = 0.f; }
            if (pv > hv) { xp(axis) = hv; if (vp(axis) > 0.f) vp(axis) = 0.f; }
        }

        // SDF obstacles
        for (const auto& sdf : obstacles_) {
            Scalar d = sdf(xp);
            if (d < 0.f) {
                // Estimate normal from finite difference
                Vec3 n;
                {Vec3 d=Vec3::Zero(); d.x()=dx*0.01f;
                 n.x() = sdf(xp+d) - sdf(xp-d);}
                {Vec3 d=Vec3::Zero(); d.y()=dx*0.01f;
                 n.y() = sdf(xp+d) - sdf(xp-d);}
                {Vec3 d=Vec3::Zero(); d.z()=dx*0.01f;
                 n.z() = sdf(xp+d) - sdf(xp-d);}
                if (n.squaredNorm() > 1e-8f) n.normalize();

                // Push out
                xp -= d * n;

                // Remove penetrating velocity
                Scalar vn = vp.dot(n);
                if (vn < 0.f) {
                    vp -= vn * n;
                    vp *= (1.f - params.boundaryFriction);
                }
            }
        }
    }
}

// -------------------------------------------------------
// Factory: box SDF
// -------------------------------------------------------
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
