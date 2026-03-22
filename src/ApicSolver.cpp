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
// Public step APIs
// -------------------------------------------------------
void ApicSolver::step() {
    Scalar substepDt = (params_.dt * params_.timeScale) / params_.substeps;

    grid_.clear();

    // 1. P2G
    p2gTransfer();

    // 2. Grid update (forces + boundary)
    gridUpdate(substepDt);

    // Snapshot velocity for FLIP increment (before pressure solve)
    grid_.snapshotVelocity();

    // 3. Pressure projection
    pressureSolve(substepDt);

    // 4. G2P
    g2pTransfer();

    // 5. Advect
    advectParticles(substepDt);

    // Hard clamp: keep particles inside grid domain
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

    // Resolve particle collisions
    if (collision_) {
        collision_->resolveParticleCollisions(particles_, params_);
    }
}

void ApicSolver::stepFrame() {
    for (int s = 0; s < params_.substeps; ++s) {
        step();
    }
}

// -------------------------------------------------------
// Step 1: P2G transfer (dispatch by method)
// -------------------------------------------------------
void ApicSolver::p2gTransfer() {
    switch (params_.method) {
    case TransferMethod::APIC:
    case TransferMethod::Hybrid:
        p2gAPIC();  break;
    case TransferMethod::FLIP:
    case TransferMethod::PIC:
    default:
        p2gPIC();   break;
    }
    grid_.normaliseMomentum();
}

// APIC P2G  (Eq. 8 in Jiang et al.)
// m_i += w_ip * m_p
// m_i * v_i += w_ip * m_p * (v_p + B_p * D_p^{-1} * (x_i - x_p))
//
// For quadratic kernel: D_p = (1/4) * dx^2 * I
// => D_p^{-1} = (4/dx^2) * I   (scalar scaling)
void ApicSolver::p2gAPIC() {
    const Scalar dx = grid_.dx();
    const Scalar Dp_inv = 4.0f / (dx * dx);  // D_p^{-1} for quadratic kernel

    std::vector<ApicGrid::WeightEntry> weights;

    const size_t np = particles_.size();
    for (size_t p = 0; p < np; ++p) {
        const Vec3&   xp = particles_.position(p);
        const Vec3&   vp = particles_.velocity(p);
        const Scalar  mp = particles_.mass(p);
        const Mat3&   Bp = particles_.affineB(p);

        grid_.gatherWeights(xp, weights);

        for (const auto& we : weights) {
            // Affine contribution: C_p * (x_i - x_p) = B_p * D_p^{-1} * offset
            Vec3 affineContrib = Bp * (Dp_inv * we.offset);

            grid_.mass(we.flatIdx)     += we.weight * mp;
            grid_.velocity(we.flatIdx) += we.weight * mp * (vp + affineContrib);
        }
    }
}

// PIC P2G -- pure weighted mass/momentum transfer
void ApicSolver::p2gPIC() {
    std::vector<ApicGrid::WeightEntry> weights;

    const size_t np = particles_.size();
    for (size_t p = 0; p < np; ++p) {
        const Vec3&   xp = particles_.position(p);
        const Vec3&   vp = particles_.velocity(p);
        const Scalar  mp = particles_.mass(p);

        grid_.gatherWeights(xp, weights);
        for (const auto& we : weights) {
            grid_.mass(we.flatIdx)     += we.weight * mp;
            grid_.velocity(we.flatIdx) += we.weight * mp * vp;
        }
    }
}

// -------------------------------------------------------
// Step 2: Grid update -- apply forces, enforce BC
// -------------------------------------------------------
void ApicSolver::gridUpdate(Scalar substepDt) {
    const size_t n = grid_.numNodes();
    for (size_t i = 0; i < n; ++i) {
        if (!grid_.fluid(i) || grid_.solid(i)) continue;
        // Apply gravity
        grid_.velocity(i) += substepDt * params_.gravity;
        // Viscosity is applied implicitly in the pressure solve for now;
        // explicit viscosity could be added here as a Laplacian term.
    }

    // Enforce velocity BC on boundaries
    if (collision_) {
        collision_->enforceGridBoundary(grid_, params_);
    }
}

// -------------------------------------------------------
// Step 3: Pressure solve
// -------------------------------------------------------
void ApicSolver::pressureSolve(Scalar substepDt) {
    // Pass effective dt to the solver (used to scale divergence)
    SimParams pparams = params_;
    pparams.dt = substepDt;
    pressureResidual_ = pressureSolver_.solve(grid_, pparams);
    printf("pressure residual=%.4e  substepDt=%.4f  density=%.1f\n",
        pressureResidual_, substepDt, params_.density);
}

// -------------------------------------------------------
// Step 4: G2P transfer
// -------------------------------------------------------
void ApicSolver::g2pTransfer() {
    switch (params_.method) {
    case TransferMethod::APIC:
        g2pAPIC();  break;
    case TransferMethod::PIC:
        g2pPIC();   break;
    case TransferMethod::FLIP:
        g2pFLIP();  break;
    case TransferMethod::Hybrid:
        // Hybrid: alpha * APIC + (1-alpha) * PIC -- we achieve this
        // by blending the resulting velocities after both G2P calls.
        // Simpler: just run APIC (the B matrix encodes APIC blend).
        g2pAPIC();  break;
    }
}

// APIC G2P  (Eq. 10 in Jiang et al.)
// v_p^{n+1} = sum_i w_ip * v_i^{n+1}
// B_p^{n+1} = sum_i w_ip * v_i^{n+1} * (x_i - x_p)^T
void ApicSolver::g2pAPIC() {
    std::vector<ApicGrid::WeightEntry> weights;
    const size_t np = particles_.size();

    for (size_t p = 0; p < np; ++p) {
        const Vec3& xp = particles_.position(p);
        grid_.gatherWeights(xp, weights);

        Vec3 vNew  = Vec3::Zero();
        Mat3 Bnew  = Mat3::Zero();

        for (const auto& we : weights) {
            const Vec3& vi = grid_.velocity(we.flatIdx);
            vNew  += we.weight * vi;
            Bnew  += we.weight * vi * we.offset.transpose();
        }

        particles_.velocity(p) = vNew;
        particles_.affineB(p)  = Bnew;
    }
}

// PIC G2P -- interpolate velocity, zero out B
void ApicSolver::g2pPIC() {
    std::vector<ApicGrid::WeightEntry> weights;
    const size_t np = particles_.size();

    for (size_t p = 0; p < np; ++p) {
        const Vec3& xp = particles_.position(p);
        grid_.gatherWeights(xp, weights);

        Vec3 vNew = Vec3::Zero();
        for (const auto& we : weights) {
            vNew += we.weight * grid_.velocity(we.flatIdx);
        }
        particles_.velocity(p) = vNew;
        particles_.affineB(p).setZero();
    }
}

// FLIP G2P -- increment v_p by delta v from grid, preserve B
void ApicSolver::g2pFLIP() {
    std::vector<ApicGrid::WeightEntry> weights;
    const size_t np = particles_.size();

    for (size_t p = 0; p < np; ++p) {
        const Vec3& xp = particles_.position(p);
        grid_.gatherWeights(xp, weights);

        Vec3 dv = Vec3::Zero();
        for (const auto& we : weights) {
            // Delta velocity = v^{n+1} - v^{n} (both on grid)
            Vec3 delta = grid_.velocity(we.flatIdx) - grid_.velocityOld(we.flatIdx);
            dv += we.weight * delta;
        }
        particles_.velocity(p) += dv;
        // B is not updated in FLIP mode
    }
}

// -------------------------------------------------------
// Step 5: Advect particles (forward Euler)
// -------------------------------------------------------
void ApicSolver::advectParticles(Scalar substepDt) {
    const Scalar dx = grid_.dx();
    const Scalar maxSpeed = 0.5f * dx / substepDt;  // CFL < 0.5
    const size_t np = particles_.size();
    for (size_t p = 0; p < np; ++p) {
        // Clamp velocity to CFL limit
        Vec3& vp = particles_.velocity(p);
        Scalar speed = vp.norm();
        if (speed > maxSpeed)
            vp *= (maxSpeed / speed);
        particles_.position(p) += substepDt * vp;
    }
}

// -------------------------------------------------------
// Diagnostics
// -------------------------------------------------------
Scalar ApicSolver::totalKineticEnergy() const {
    Scalar E = 0.f;
    const size_t np = particles_.size();
    for (size_t p = 0; p < np; ++p) {
        E += 0.5f * particles_.mass(p) * particles_.velocity(p).squaredNorm();
    }
    return E;
}

Vec3 ApicSolver::totalAngularMomentum() const {
    return particles_.totalAngularMomentum();
}

} // namespace apic
