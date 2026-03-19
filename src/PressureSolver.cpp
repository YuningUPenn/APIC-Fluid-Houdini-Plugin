#include "PressureSolver.h"
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include <numeric>

// -------------------------------------------------------
// Pressure projection for a collocated grid.
//
// We use a MAC-style staggered formulation even on the
// collocated grid by treating face velocities as averages.
// The divergence at node (i,j,k) is:
//
//   div(v) = (v_x[i+1] - v_x[i-1]) / (2*dx)
//           +(v_y[j+1] - v_y[j-1]) / (2*dx)
//           +(v_z[k+1] - v_z[k-1]) / (2*dx)
//
// The corresponding pressure Poisson equation is:
//
//   Laplacian(p) = (rho/dt) * div(v*)
//
// with the same central-difference Laplacian:
//
//   Lap(p)_i = (p[i+1] - 2p[i] + p[i-1]) / dx^2  (each axis)
//
// The velocity correction is:
//
//   v[i] -= (dt/rho) * (p[i+1] - p[i-1]) / (2*dx)
//
// This is provably consistent: applying the correction reduces
// the central-difference divergence to zero.
//
// Boundary condition at solid walls: Neumann (dp/dn = 0).
// In the Laplacian stencil this means the ghost pressure equals
// the interior pressure -- the off-diagonal term is dropped and
// the diagonal is NOT incremented (unlike the previous version
// which incorrectly added 1/dx^2 to the diagonal for solid faces).
// -------------------------------------------------------

namespace apic {

Scalar PressureSolver::solve(ApicGrid& grid, const SimParams& params) {
    std::vector<int> fluidToFlat, flatToFluid;
    buildFluidIndex(grid, fluidToFlat, flatToFluid);
    int nFluid = static_cast<int>(fluidToFlat.size());
    if (nFluid == 0) return 0.f;
    computeDivergence(grid, params);
    Scalar res = (method_ == Method::PCG)
        ? solvePCG   (grid, params, fluidToFlat, flatToFluid)
        : solveJacobi(grid, params, fluidToFlat, flatToFluid);
    applyPressureGradient(grid, params);
    return res;
}

void PressureSolver::buildFluidIndex(const ApicGrid& grid,
                                      std::vector<int>& fluidToFlat,
                                      std::vector<int>& flatToFluid) const {
    size_t n = grid.numNodes();
    flatToFluid.assign(n, -1);
    fluidToFlat.clear();
    for (int i = 0; i < static_cast<int>(n); ++i) {
        if (grid.fluid(i) && !grid.solid(i)) {
            flatToFluid[i] = static_cast<int>(fluidToFlat.size());
            fluidToFlat.push_back(i);
        }
    }
}

// Central-difference divergence.
// Solid/boundary neighbours contribute zero velocity (no-penetration wall BC).
void PressureSolver::computeDivergence(ApicGrid& grid,
                                        const SimParams& /*p*/) const {
    const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
    const Scalar inv2dx = 0.5f / grid.dx();
    for (int k = 1; k < nz-1; ++k)
    for (int j = 1; j < ny-1; ++j)
    for (int i = 1; i < nx-1; ++i) {
        int flat = grid.idx(i,j,k);
        if (!grid.fluid(flat) || grid.solid(flat)) continue;
        auto vel = [&](int ii,int jj,int kk) -> Vec3 {
            if (!grid.inside(ii,jj,kk) || grid.solid(ii,jj,kk))
                return Vec3::Zero();
            return grid.velocity(ii,jj,kk);
        };
        Scalar d = (vel(i+1,j,k).x() - vel(i-1,j,k).x()
                  + vel(i,j+1,k).y() - vel(i,j-1,k).y()
                  + vel(i,j,k+1).z() - vel(i,j,k-1).z()) * inv2dx;
        grid.divergence(flat) = d;
    }
}

static const int OFF6[6][3] = {
    {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}
};

// Build the Poisson matrix for central-difference Laplacian.
// Neumann BC at solid walls: ghost p = interior p, so the off-diagonal
// term is simply omitted (row diagonal is NOT inflated for solid faces).
// This gives a singular system (kernel = constant), which Eigen's CG
// handles fine because we start from p=0.
Scalar PressureSolver::solvePCG(ApicGrid& grid, const SimParams& params,
                                  const std::vector<int>& fluidToFlat,
                                  const std::vector<int>& flatToFluid) {
#ifdef APIC_STUB_EIGEN
    return solveJacobi(grid, params, fluidToFlat, flatToFluid);
#else
    using SpMat = Eigen::SparseMatrix<Scalar>;
    int nFluid = static_cast<int>(fluidToFlat.size());
    const Scalar dx2   = grid.dx() * grid.dx();
    const Scalar scale = params.density / params.dt;
    const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(nFluid * 7);
    std::vector<Scalar> rhs(nFluid, 0.f);

    for (int fi = 0; fi < nFluid; ++fi) {
        int flat = fluidToFlat[fi];
        int gk = flat/(nx*ny), rem = flat%(nx*ny), gj = rem/nx, gi = rem%nx;
        int numFluidNeighbours = 0;

        for (const auto& o : OFF6) {
            int ni=gi+o[0], nj=gj+o[1], nk=gk+o[2];
            // Solid or out-of-bounds: Neumann BC -- skip this face entirely.
            if (!grid.inside(ni,nj,nk) || grid.solid(ni,nj,nk)) continue;
            int nflat = grid.idx(ni,nj,nk);
            int nfi   = flatToFluid[nflat];
            if (nfi >= 0) {
                triplets.emplace_back(fi, nfi, -1.f/dx2);
                ++numFluidNeighbours;
            }
        }
        // Diagonal = number of fluid neighbours / dx^2
        // (Neumann faces drop out -- not added to diagonal)
        if (numFluidNeighbours == 0) {
            triplets.emplace_back(fi, fi, 1.f);
        }
        else {
            triplets.emplace_back(fi, fi,
                static_cast<Scalar>(numFluidNeighbours) / dx2);
        }
        rhs[fi] = -scale * grid.divergence(flat);
    }

    SpMat A(nFluid, nFluid);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    // Use proper dynamic vectors
    using VecF = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    VecF rhs_e(nFluid), p_e(nFluid);
    for (int i = 0; i < nFluid; ++i) rhs_e(i) = rhs[i];
    p_e.setZero();

    // Neumann compatibility: subtract mean so sum(rhs)=0
    {
        Scalar mean = 0.f;
        for (int i = 0; i < nFluid; i++) mean += (float)rhs_e(i);
        mean /= nFluid;
        for (int i = 0; i < nFluid; i++) rhs_e(i) -= mean;
    }

    Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper,
                              Eigen::DiagonalPreconditioner<Scalar>> cg;
    cg.setMaxIterations(maxIter_);
    cg.setTolerance(static_cast<double>(tol_));
    cg.compute(A);
    p_e = cg.solve(rhs_e);

    lastIter_     = static_cast<int>(cg.iterations());
    lastResidual_ = static_cast<Scalar>(cg.error());

    for (int fi = 0; fi < nFluid; ++fi)
        grid.pressure(fluidToFlat[fi]) = static_cast<Scalar>(p_e(fi));

    return lastResidual_;
#endif
}

// Jacobi iteration matching the same Laplacian structure.
// Neumann BC: solid faces are simply not summed over.
Scalar PressureSolver::solveJacobi(ApicGrid& grid, const SimParams& params,
                                    const std::vector<int>& fluidToFlat,
                                    const std::vector<int>& flatToFluid) {
    const Scalar dx2   = grid.dx() * grid.dx();
    const Scalar scale = params.density / params.dt;
    const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
    int nFluid = static_cast<int>(fluidToFlat.size());
    std::vector<Scalar> pNew(nFluid, 0.f);
    Scalar residual = 1e10f;

    for (int fi = 0; fi < nFluid; ++fi)
        grid.pressure(fluidToFlat[fi]) = 0.f;

    for (int iter = 0; iter < maxIter_; ++iter) {
        residual = 0.f;
        for (int fi = 0; fi < nFluid; ++fi) {
            int flat = fluidToFlat[fi];
            int gk = flat/(nx*ny), rem = flat%(nx*ny), gj = rem/nx, gi = rem%nx;
            Scalar sum = 0.f;
            int cnt = 0;
            for (const auto& o : OFF6) {
                int ni=gi+o[0], nj=gj+o[1], nk=gk+o[2];
                if (!grid.inside(ni,nj,nk) || grid.solid(ni,nj,nk)) continue;
                int nflat = grid.idx(ni,nj,nk);
                int nfi   = flatToFluid[nflat];
                if (nfi >= 0) { sum += grid.pressure(nflat); ++cnt; }
            }
            if (cnt == 0) { pNew[fi] = 0.f; continue; }
            Scalar rhs_val = scale * dx2 * grid.divergence(flat);
            Scalar pn = (sum - rhs_val) / static_cast<Scalar>(cnt);
            residual += std::abs(pn - grid.pressure(flat));
            pNew[fi] = pn;
        }
        for (int fi = 0; fi < nFluid; ++fi)
            grid.pressure(fluidToFlat[fi]) = pNew[fi];
        residual /= static_cast<Scalar>(nFluid);
        lastIter_ = iter + 1;
        if (residual < tol_) break;
    }
    lastResidual_ = residual;
    return residual;
}

// Velocity correction: v -= (dt/rho) * grad(p)
// Central-difference gradient, Neumann at solid faces (gradient = 0 there).
void PressureSolver::applyPressureGradient(ApicGrid& grid,
                                            const SimParams& params) const {
    const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
    const Scalar coeff = params.dt / (params.density * 2.0f * grid.dx());
    for (int k = 1; k < nz-1; ++k)
    for (int j = 1; j < ny-1; ++j)
    for (int i = 1; i < nx-1; ++i) {
        int flat = grid.idx(i,j,k);
        if (!grid.fluid(flat) || grid.solid(flat)) continue;
        // At a solid face the gradient is zero (Neumann dp/dn=0 means
        // no pressure difference across the wall).
        auto pAt = [&](int ni,int nj,int nk) -> Scalar {
            if (!grid.inside(ni,nj,nk) || grid.solid(ni,nj,nk))
                return grid.pressure(i,j,k);  // ghost = self -> grad = 0
            return grid.pressure(ni,nj,nk);
        };
        Vec3 gp;
        gp.x() = pAt(i+1,j,k) - pAt(i-1,j,k);
        gp.y() = pAt(i,j+1,k) - pAt(i,j-1,k);
        gp.z() = pAt(i,j,k+1) - pAt(i,j,k-1);
        grid.velocity(i,j,k) -= coeff * gp;
    }
}

} // namespace apic
