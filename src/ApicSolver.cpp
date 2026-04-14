#include "ApicSolver.h"
#include <cmath>
#include <cassert>

namespace apic {

    // -------------------------------------------------------
    // Construction / setup
    // -------------------------------------------------------
    ApicSolver::ApicSolver(const SimParams& params) : params_(params) {
        grid_.init(params_.gridResX, params_.gridResY, params_.gridResZ,
            params_.gridSpacing, params_.gridOrigin);
    }

    void ApicSolver::setParams(const SimParams& params) {
        params_ = params;
        grid_.init(params_.gridResX, params_.gridResY, params_.gridResZ,
            params_.gridSpacing, params_.gridOrigin);
    }

    void ApicSolver::reset() {
        grid_.clear();
        particles_.clear();
        pressureResidual_ = 0.f;
    }

    // -------------------------------------------------------
    // Main step
    // -------------------------------------------------------
    void ApicSolver::step() {
        Scalar substepDt = (params_.dt * params_.timeScale) / params_.substeps;

        grid_.clear();

        // 1. P2G
        p2gTransfer();

        // Snapshot BEFORE gravity for FLIP delta
        grid_.snapshotVelocity();

        // 2. Grid update (gravity, BC)
        gridUpdate(substepDt);

        // 3. Pressure projection
        pressureSolve(substepDt);

        // 4. G2P
        g2pTransfer();

        // 5. Advect
        advectParticles(substepDt);

        // Hard clamp: keep particles inside domain
        {
            const float lo = grid_.origin().x() + grid_.dx();
            const float hi = grid_.origin().x() + (grid_.nx() - 2) * grid_.dx();
            for (size_t p = 0; p < particles_.size(); ++p) {
                Vec3& xp = particles_.position(p);
                Vec3& vp = particles_.velocity(p);
                for (int axis = 0; axis < 3; ++axis) {
                    if (xp(axis) < lo) { xp(axis) = lo; if (vp(axis) < 0) vp(axis) = 0; }
                    if (xp(axis) > hi) { xp(axis) = hi; if (vp(axis) > 0) vp(axis) = 0; }
                }
            }
        }

        if (collision_)
            collision_->resolveParticleCollisions(particles_, params_);
    }

    void ApicSolver::stepFrame() {
        for (int s = 0; s < params_.substeps; ++s) step();
    }

    // -------------------------------------------------------
    // Step 1: P2G  -- MAC staggered transfer
    //
    // MAC APIC (paper Section 6):
    //   For each axis a in {x,y,z}:
    //     m_a[face] += w_aip * mp
    //     m_a[face] * v_a[face] += w_aip * mp * (vp[a] + C_pa . (xface - xp))
    //   where C_pa = affineB.col(a)  (stored directly, no Dp_inv)
    // -------------------------------------------------------
    void ApicSolver::p2gTransfer() {
        switch (params_.method) {
        case TransferMethod::APIC:
        case TransferMethod::Hybrid:
            p2gAPIC(); break;
        case TransferMethod::FLIP:
        case TransferMethod::PIC:
        default:
            p2gPIC(); break;
        }
        grid_.normaliseFaceMomentum();
    }

    void ApicSolver::p2gAPIC() {
        std::vector<ApicGrid::WeightEntry> weights;
        const size_t np = particles_.size();

        for (size_t p = 0; p < np; ++p) {
            const Vec3& xp = particles_.position(p);
            const Vec3& vp = particles_.velocity(p);
            const Scalar mp = particles_.mass(p);
            const Mat3& Cp = particles_.affineB(p);

            for (int axis = 0; axis < 3; ++axis) {
                grid_.gatherWeightsFace(xp, axis, weights);
                auto& faceVel = (axis == 0) ? grid_.uFaceVelArr() : (axis == 1) ? grid_.vFaceVelArr() : grid_.wFaceVelArr();
                auto& faceMass = (axis == 0) ? grid_.uFaceMassArr() : (axis == 1) ? grid_.vFaceMassArr() : grid_.wFaceMassArr();
                for (const auto& we : weights) {
                    Scalar affineContrib = Cp.col(axis).dot(we.offset);
                    Scalar vel_a = vp(axis) + affineContrib;
                    faceMass[we.flatIdx] += we.weight * mp;
                    faceVel[we.flatIdx] += we.weight * mp * vel_a;
                }
            }
        }
    }

    void ApicSolver::p2gPIC() {
        std::vector<ApicGrid::WeightEntry> weights;
        const size_t np = particles_.size();

        for (size_t p = 0; p < np; ++p) {
            const Vec3& xp = particles_.position(p);
            const Vec3& vp = particles_.velocity(p);
            const Scalar mp = particles_.mass(p);

            for (int axis = 0; axis < 3; ++axis) {
                grid_.gatherWeightsFace(xp, axis, weights);
                auto& faceVel = (axis == 0) ? grid_.uFaceVelArr() : (axis == 1) ? grid_.vFaceVelArr() : grid_.wFaceVelArr();
                auto& faceMass = (axis == 0) ? grid_.uFaceMassArr() : (axis == 1) ? grid_.vFaceMassArr() : grid_.wFaceMassArr();
                for (const auto& we : weights) {
                    faceMass[we.flatIdx] += we.weight * mp;
                    faceVel[we.flatIdx] += we.weight * mp * vp(axis);
                }
            }
        }
    }

    // -------------------------------------------------------
    // Step 2: Grid update
    // Apply gravity and viscosity damping to face velocities
    // -------------------------------------------------------
    void ApicSolver::gridUpdate(Scalar substepDt) {
        // Apply gravity to each face axis
        const Vec3& g = params_.gravity;
        const Scalar visc = 1.f - params_.viscosity * substepDt;

        // u-faces (x-component)
        for (int k = 0; k < grid_.nUz(); ++k)
            for (int j = 0; j < grid_.nUy(); ++j)
                for (int i = 0; i < grid_.nUx(); ++i) {
                    Scalar& u = grid_.uFaceVel(i, j, k);
                    u += substepDt * g.x();
                    u *= visc;
                }
        // v-faces (y-component)
        for (int k = 0; k < grid_.nVz(); ++k)
            for (int j = 0; j < grid_.nVy(); ++j)
                for (int i = 0; i < grid_.nVx(); ++i) {
                    Scalar& v = grid_.vFaceVel(i, j, k);
                    v += substepDt * g.y();
                    v *= visc;
                }
        // w-faces (z-component)
        for (int k = 0; k < grid_.nWz(); ++k)
            for (int j = 0; j < grid_.nWy(); ++j)
                for (int i = 0; i < grid_.nWx(); ++i) {
                    Scalar& w = grid_.wFaceVel(i, j, k);
                    w += substepDt * g.z();
                    w *= visc;
                }

        // Enforce solid wall BC (zero normal velocity at boundaries)
        grid_.enforceBoundaryFaces();

        if (collision_)
            collision_->enforceGridBoundary(grid_, params_);
    }

    // -------------------------------------------------------
    // Step 3: Pressure solve
    // -------------------------------------------------------
    void ApicSolver::pressureSolve(Scalar substepDt) {
        SimParams pparams = params_;
        pparams.dt = substepDt;
        pressureSolver_.setMethod(PressureSolver::Method::PCG);
        pressureResidual_ = pressureSolver_.solve(grid_, pparams);

        // Diagnostic print
        const int nx = grid_.nx(), ny = grid_.ny(), nz = grid_.nz();
        const Scalar inv_dx = 1.f / grid_.dx();
        Scalar maxDivBefore = 0.f, maxDivAfter = 0.f;
        size_t nFluid = 0;
        bool printed = false;

        for (int k = 1; k < nz - 1; ++k)
            for (int j = 1; j < ny - 1; ++j)
                for (int i = 1; i < nx - 1; ++i) {
                    int flat = grid_.idx(i, j, k);
                    if (!grid_.fluid(flat) || grid_.solid(flat)) continue;
                    ++nFluid;
                    maxDivBefore = std::max(maxDivBefore, std::abs(grid_.divergence(flat)));

                    // Re-compute divergence after projection
                    Scalar u_r = grid_.uInside(i, j, k) ? grid_.uFaceVel(i, j, k) : 0.f;
                    Scalar u_l = grid_.uInside(i - 1, j, k) ? grid_.uFaceVel(i - 1, j, k) : 0.f;
                    Scalar v_t = grid_.vInside(i, j, k) ? grid_.vFaceVel(i, j, k) : 0.f;
                    Scalar v_b = grid_.vInside(i, j - 1, k) ? grid_.vFaceVel(i, j - 1, k) : 0.f;
                    Scalar w_f = grid_.wInside(i, j, k) ? grid_.wFaceVel(i, j, k) : 0.f;
                    Scalar w_k = grid_.wInside(i, j, k - 1) ? grid_.wFaceVel(i, j, k - 1) : 0.f;
                    Scalar da = (u_r - u_l + v_t - v_b + w_f - w_k) * inv_dx;
                    maxDivAfter = std::max(maxDivAfter, std::abs(da));

                    if (!printed) {
                        printf("  sample: divBefore=%.4f  p=%.4f\n",
                            grid_.divergence(flat), grid_.pressure(flat));
                        printed = true;
                    }
                }
        printf("before=%.4f  after=%.4f  residual=%.2e  nFluid=%zu\n",
            maxDivBefore, maxDivAfter, pressureResidual_, nFluid);
    }

    // -------------------------------------------------------
    // Step 4: G2P  -- MAC staggered transfer
    //
    // MAC APIC (paper Eq.14):
    //   vp[a] = sum_i w_aip * v_ai
    //   C_pa  = sum_i grad(w_aip) * v_ai   (no Dp_inv)
    // -------------------------------------------------------
    void ApicSolver::g2pTransfer() {
        switch (params_.method) {
        case TransferMethod::APIC:    g2pAPIC(); break;
        case TransferMethod::PIC:     g2pPIC();  break;
        case TransferMethod::FLIP:    g2pFLIP(); break;
        case TransferMethod::Hybrid:  g2pAPIC(); break;
        }
    }

    void ApicSolver::g2pAPIC() {
        std::vector<ApicGrid::WeightEntry> weights;
        const size_t np = particles_.size();

        for (size_t p = 0; p < np; ++p) {
            const Vec3& xp = particles_.position(p);
            Vec3 vNew = Vec3::Zero();
            Mat3 Cnew = Mat3::Zero();  // columns: c_px, c_py, c_pz

            for (int axis = 0; axis < 3; ++axis) {
                grid_.gatherWeightsFace(xp, axis, weights);
                for (const auto& we : weights) {
                    Scalar faceVel = 0.f;
                    if (axis == 0) faceVel = grid_.uFaceVelArr()[we.flatIdx];
                    else if (axis == 1) faceVel = grid_.vFaceVelArr()[we.flatIdx];
                    else              faceVel = grid_.wFaceVelArr()[we.flatIdx];

                    vNew(axis) += we.weight * faceVel;
                    Cnew.col(axis) += we.gradient * faceVel;
                }
            }

            particles_.velocity(p) = vNew;
            particles_.affineB(p) = Cnew;
        }
    }

    void ApicSolver::g2pPIC() {
        std::vector<ApicGrid::WeightEntry> weights;
        const size_t np = particles_.size();

        for (size_t p = 0; p < np; ++p) {
            const Vec3& xp = particles_.position(p);
            Vec3 vNew = Vec3::Zero();

            for (int axis = 0; axis < 3; ++axis) {
                grid_.gatherWeightsFace(xp, axis, weights);
                for (const auto& we : weights) {
                    Scalar faceVel = 0.f;
                    if (axis == 0) faceVel = grid_.uFaceVelArr()[we.flatIdx];
                    else if (axis == 1) faceVel = grid_.vFaceVelArr()[we.flatIdx];
                    else              faceVel = grid_.wFaceVelArr()[we.flatIdx];
                    vNew(axis) += we.weight * faceVel;
                }
            }

            particles_.velocity(p) = vNew;
            particles_.affineB(p).setZero();
        }
    }

    void ApicSolver::g2pFLIP() {
        std::vector<ApicGrid::WeightEntry> weights;
        const size_t np = particles_.size();

        for (size_t p = 0; p < np; ++p) {
            const Vec3& xp = particles_.position(p);
            Vec3 dv = Vec3::Zero();
            Mat3 Cnew = Mat3::Zero();

            for (int axis = 0; axis < 3; ++axis) {
                grid_.gatherWeightsFace(xp, axis, weights);
                const auto& faceNew = (axis == 0) ? grid_.uFaceVelArr() : (axis == 1) ? grid_.vFaceVelArr() : grid_.wFaceVelArr();
                const auto& faceOld = (axis == 0) ? grid_.uFaceOldArr() : (axis == 1) ? grid_.vFaceOldArr() : grid_.wFaceOldArr();
                for (const auto& we : weights) {
                    Scalar fNew = faceNew[we.flatIdx];
                    Scalar fOld = faceOld[we.flatIdx];
                    dv(axis) += we.weight * (fNew - fOld);
                    Cnew.col(axis) += we.gradient * fNew;
                }
            }

            particles_.velocity(p) += dv;
            particles_.affineB(p) = Cnew;
        }
    }

    // -------------------------------------------------------
    // Step 5: Advect
    // -------------------------------------------------------
    void ApicSolver::advectParticles(Scalar substepDt) {
        const Scalar dx = grid_.dx();
        const Scalar maxSpeed = 0.5f * dx / substepDt;
        const size_t np = particles_.size();
        for (size_t p = 0; p < np; ++p) {
            Vec3& vp = particles_.velocity(p);
            Scalar speed = vp.norm();
            if (speed > maxSpeed) vp *= (maxSpeed / speed);
            particles_.position(p) += substepDt * vp;
        }
    }

    // -------------------------------------------------------
    // Diagnostics
    // -------------------------------------------------------
    Scalar ApicSolver::totalKineticEnergy() const {
        Scalar E = 0.f;
        for (size_t p = 0; p < particles_.size(); ++p)
            E += 0.5f * particles_.mass(p) * particles_.velocity(p).squaredNorm();
        return E;
    }

    Vec3 ApicSolver::totalAngularMomentum() const {
        return particles_.totalAngularMomentum();
    }

} // namespace apic
