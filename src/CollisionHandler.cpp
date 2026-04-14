#include "CollisionHandler.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

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

        const int nx = grid.nx();
        const int ny = grid.ny();
        const int nz = grid.nz();

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    Vec3 pos = grid.nodePos(i, j, k);
                    if (evalSDF(pos) <= 0.f) {
                        grid.solid(i, j, k) = 1;
                    }
                }
            }
        }
    }

    // -------------------------------------------------------
    // Enforce obstacle-aware grid boundary conditions on MAC faces.
    // Zero face velocity if either adjacent node is solid.
    // -------------------------------------------------------
    void CollisionHandler::enforceGridBoundary(ApicGrid& grid,
        const SimParams& /*params*/) const {

        // u-faces lie between node (i,j,k) and (i+1,j,k)
        for (int k = 0; k < grid.nUz(); ++k) {
            for (int j = 0; j < grid.nUy(); ++j) {
                for (int i = 0; i < grid.nUx(); ++i) {
                    const int n0 = grid.idx(i, j, k);
                    const int n1 = grid.idx(i + 1, j, k);
                    if (grid.solid(n0) || grid.solid(n1)) {
                        grid.uFaceVel(i, j, k) = 0.f;
                    }
                }
            }
        }

        // v-faces lie between node (i,j,k) and (i,j+1,k)
        for (int k = 0; k < grid.nVz(); ++k) {
            for (int j = 0; j < grid.nVy(); ++j) {
                for (int i = 0; i < grid.nVx(); ++i) {
                    const int n0 = grid.idx(i, j, k);
                    const int n1 = grid.idx(i, j + 1, k);
                    if (grid.solid(n0) || grid.solid(n1)) {
                        grid.vFaceVel(i, j, k) = 0.f;
                    }
                }
            }
        }

        // w-faces lie between node (i,j,k) and (i,j,k+1)
        for (int k = 0; k < grid.nWz(); ++k) {
            for (int j = 0; j < grid.nWy(); ++j) {
                for (int i = 0; i < grid.nWx(); ++i) {
                    const int n0 = grid.idx(i, j, k);
                    const int n1 = grid.idx(i, j, k + 1);
                    if (grid.solid(n0) || grid.solid(n1)) {
                        grid.wFaceVel(i, j, k) = 0.f;
                    }
                }
            }
        }

        // Also keep the outer domain walls enforced
        grid.enforceBoundaryFaces();
    }

    // -------------------------------------------------------
    // Resolve particle-solid collisions
    // -------------------------------------------------------
    void CollisionHandler::resolveParticleCollisions(ApicParticleSet& particles,
        const SimParams& params) const {
        if (obstacles_.empty()) return;

        const Scalar dx = params.gridSpacing;
        const Scalar eps = std::max(dx * 0.01f, 1e-4f);
        const size_t np = particles.size();

        for (size_t p = 0; p < np; ++p) {
            Vec3& xp = particles.position(p);
            Vec3& vp = particles.velocity(p);

            for (const auto& sdf : obstacles_) {
                const Scalar d = sdf(xp);
                if (d < 0.f) {
                    Vec3 n = Vec3::Zero();

                    {
                        Vec3 e = Vec3::Zero();
                        e.x() = eps;
                        n.x() = sdf(xp + e) - sdf(xp - e);
                    }
                    {
                        Vec3 e = Vec3::Zero();
                        e.y() = eps;
                        n.y() = sdf(xp + e) - sdf(xp - e);
                    }
                    {
                        Vec3 e = Vec3::Zero();
                        e.z() = eps;
                        n.z() = sdf(xp + e) - sdf(xp - e);
                    }

                    const Scalar n2 = n.squaredNorm();
                    if (n2 > 1e-12f) {
                        n /= std::sqrt(n2);
                    }
                    else {
                        n = Vec3(0.f, 1.f, 0.f);
                    }

                    // Push particle out of obstacle
                    xp -= d * n;

                    // Remove inward normal velocity, keep tangential motion with friction
                    const Scalar vn = vp.dot(n);
                    if (vn < 0.f) {
                        const Vec3 vNormal = vn * n;
                        const Vec3 vTangent = vp - vNormal;
                        vp = vTangent * (1.f - params.boundaryFriction);
                    }
                }
            }
        }
    }

    // -------------------------------------------------------
    // Proper signed distance to an axis-aligned box.
    // Negative inside, positive outside, 0 on the surface.
    // -------------------------------------------------------
    CollisionHandler::SDFFunc CollisionHandler::makeBoxSDF(const Vec3& minPt,
        const Vec3& maxPt) {
        const Vec3 center = Scalar(0.5f) * (minPt + maxPt);
        const Vec3 halfExtent = Scalar(0.5f) * (maxPt - minPt);

        return [center, halfExtent](const Vec3& p) -> Scalar {
            const Vec3 q = (p - center).cwiseAbs() - halfExtent;
            const Vec3 qMax = q.cwiseMax(Scalar(0));

            const Scalar outsideDist = qMax.norm();
            const Scalar insideDist = std::min(
                std::max(q.x(), std::max(q.y(), q.z())),
                Scalar(0)
            );

            return outsideDist + insideDist;
            };
    }

} // namespace apic