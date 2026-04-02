#include "PressureSolver.h"
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include <numeric>

// -------------------------------------------------------
// MAC-grid pressure projection.
//
// Divergence at node (i,j,k) -- forward difference of face velocities:
//   div = (u[i,j,k] - u[i-1,j,k] + v[i,j,k] - v[i,j-1,k]
//          + w[i,j,k] - w[i,j,k-1]) / dx
//
// Poisson equation  A*p = (rho/dt)*div
// where A is the positive-definite discrete Laplacian.
//
// Velocity correction:
//   u[i,j,k] -= (dt/rho)*(p[i+1,j,k]-p[i,j,k])/dx
//   v[i,j,k] -= (dt/rho)*(p[i,j+1,k]-p[i,j,k])/dx
//   w[i,j,k] -= (dt/rho)*(p[i,j,k+1]-p[i,j,k])/dx
//
// (forward gradient is the adjoint of backward divergence => exact cancellation)
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

// -------------------------------------------------------
// MAC divergence: backward differences from face velocities.
// u[i,j,k] is the x-velocity at the face between nodes (i,j,k)
// and (i+1,j,k).  (u[-1,j,k] at solid wall = 0 by BC.)
// -------------------------------------------------------
void PressureSolver::computeDivergence(ApicGrid& grid,
                                        const SimParams& /*p*/) const {
    const int nx=grid.nx(), ny=grid.ny(), nz=grid.nz();
    const Scalar inv_dx = 1.f / grid.dx();

    for (int k=1;k<nz-1;++k)
    for (int j=1;j<ny-1;++j)
    for (int i=1;i<nx-1;++i) {
        int flat=grid.idx(i,j,k);
        if (!grid.fluid(flat)||grid.solid(flat)) continue;

        // u-face at (i,j,k) is between nodes i and i+1; face at (i-1,j,k) is to the left.
        Scalar u_right = grid.uInside(i,   j, k) ? grid.uFaceVel(i,  j, k) : 0.f;
        Scalar u_left  = grid.uInside(i-1, j, k) ? grid.uFaceVel(i-1,j, k) : 0.f;
        Scalar v_top   = grid.vInside(i, j,   k) ? grid.vFaceVel(i, j,  k) : 0.f;
        Scalar v_bot   = grid.vInside(i, j-1, k) ? grid.vFaceVel(i,j-1, k) : 0.f;
        Scalar w_front = grid.wInside(i, j, k  ) ? grid.wFaceVel(i, j,  k) : 0.f;
        Scalar w_back  = grid.wInside(i, j, k-1) ? grid.wFaceVel(i, j,k-1) : 0.f;

        grid.divergence(flat) = (u_right - u_left
                               + v_top   - v_bot
                               + w_front - w_back) * inv_dx;
    }
}

static const int OFF6[6][3] = {
    {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}
};

// -------------------------------------------------------
// Build the Laplacian matrix A (same 7-point stencil as before).
// The system is A*p = (rho/dt)*div.
// Neumann BC at solid/boundary faces: those faces are excluded
// from the stencil (pressure gradient at wall = 0).
// -------------------------------------------------------
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
        int gk=flat/(nx*ny), rem=flat%(nx*ny), gj=rem/nx, gi=rem%nx;
        int numNonSolidNeighbours = 0;  // counts both fluid AND air neighbors

        for (const auto& o : OFF6) {
            int ni=gi+o[0], nj=gj+o[1], nk=gk+o[2];
            if (!grid.inside(ni,nj,nk)||grid.solid(ni,nj,nk)) continue;
            ++numNonSolidNeighbours;  // count ALL non-solid (fluid + air)
            int nflat=grid.idx(ni,nj,nk);
            int nfi=flatToFluid[nflat];
            if (nfi >= 0) {
                // Fluid neighbor: off-diagonal entry
                triplets.emplace_back(fi, nfi, -1.f/dx2);
            }
            // Air neighbor (nfi < 0): Dirichlet p_air=0
            // contributes 1/dx2 to diagonal (already counted via numNonSolidNeighbours)
            // but NOT to rhs (since p_air=0, rhs contribution = 0)
        }
        if (numNonSolidNeighbours == 0)
            triplets.emplace_back(fi, fi, 1.f);
        else
            triplets.emplace_back(fi, fi,
                static_cast<Scalar>(numNonSolidNeighbours)/dx2);

        // A*p = -(rho/dt)*div  =>  rhs = -scale * div
        // (derived from: div_new = div - (dt/rho)*A*p = 0  =>  A*p = -(rho/dt)*div)
        rhs[fi] = -scale * grid.divergence(flat);
    }

    SpMat A(nFluid, nFluid);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    using VecF = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    VecF rhs_e(nFluid), p_e(nFluid);
    for (int i=0;i<nFluid;++i) rhs_e(i)=rhs[i];
    p_e.setZero();

    // Neumann compatibility: subtract mean so CG converges.
    // Constant pressure shift does not affect velocity gradients.
    {
        Scalar mean = 0.f;
        for (int i=0;i<nFluid;++i) mean += rhs_e(i);
        mean /= nFluid;
        for (int i=0;i<nFluid;++i) rhs_e(i) -= mean;
    }

    Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper,
                              Eigen::DiagonalPreconditioner<Scalar>> cg;
    cg.setMaxIterations(maxIter_);
    cg.setTolerance(static_cast<double>(tol_));
    cg.compute(A);
    p_e = cg.solve(rhs_e);

    lastIter_     = static_cast<int>(cg.iterations());
    lastResidual_ = static_cast<Scalar>(cg.error());

    for (int fi=0;fi<nFluid;++fi)
        grid.pressure(fluidToFlat[fi]) = static_cast<Scalar>(p_e(fi));

    return lastResidual_;
#endif
}

// -------------------------------------------------------
// Jacobi fallback
// -------------------------------------------------------
Scalar PressureSolver::solveJacobi(ApicGrid& grid, const SimParams& params,
                                    const std::vector<int>& fluidToFlat,
                                    const std::vector<int>& flatToFluid) {
    const Scalar dx2   = grid.dx() * grid.dx();
    const Scalar scale = params.density / params.dt;
    const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
    int nFluid = static_cast<int>(fluidToFlat.size());
    std::vector<Scalar> pNew(nFluid, 0.f);
    Scalar residual = 1e10f;

    for (int fi=0;fi<nFluid;++fi)
        grid.pressure(fluidToFlat[fi]) = 0.f;

    for (int iter=0;iter<maxIter_;++iter) {
        residual = 0.f;
        for (int fi=0;fi<nFluid;++fi) {
            int flat=fluidToFlat[fi];
            int gk=flat/(nx*ny), rem=flat%(nx*ny), gj=rem/nx, gi=rem%nx;
            Scalar sum=0.f; int cnt=0;
            for (const auto& o : OFF6) {
                int ni=gi+o[0], nj=gj+o[1], nk=gk+o[2];
                if (!grid.inside(ni,nj,nk)||grid.solid(ni,nj,nk)) continue;
                ++cnt;  // count ALL non-solid neighbors (fluid + air)
                int nflat=grid.idx(ni,nj,nk);
                int nfi=flatToFluid[nflat];
                if (nfi>=0) sum+=grid.pressure(nflat);
                // air neighbor: p=0, contributes 0 to sum
            }
            if (cnt==0) { pNew[fi]=0.f; continue; }
            // A*p = -scale*div  =>  cnt*p[i]/dx2 - sum/dx2 = -scale*div
            // p[i] = (sum - scale*dx2*div) / cnt
            Scalar rhs_val = scale * dx2 * grid.divergence(flat);
            Scalar pn = (sum - rhs_val) / static_cast<Scalar>(cnt);
            residual += std::abs(pn - grid.pressure(flat));
            pNew[fi] = pn;
        }
        for (int fi=0;fi<nFluid;++fi)
            grid.pressure(fluidToFlat[fi]) = pNew[fi];
        residual /= static_cast<Scalar>(nFluid);
        lastIter_ = iter+1;
        if (residual < tol_) break;
    }
    lastResidual_ = residual;
    return residual;
}

// -------------------------------------------------------
// Apply pressure gradient to face velocities.
// Forward differences: grad at u-face (i,j,k) = (p[i+1,j,k]-p[i,j,k])/dx
// Only update faces between two fluid (non-solid) nodes.
// -------------------------------------------------------
void PressureSolver::applyPressureGradient(ApicGrid& grid,
                                            const SimParams& params) const {
    const Scalar coeff = params.dt / (params.density * grid.dx());

    // u-faces: between node (i,j,k) and (i+1,j,k)
    for (int k=0;k<grid.nUz();++k)
    for (int j=0;j<grid.nUy();++j)
    for (int i=0;i<grid.nUx();++i) {
        int n0=grid.idx(i,  j,k);
        int n1=grid.idx(i+1,j,k);
        // Skip if either adjacent node is solid or neither is fluid
        if (grid.solid(n0)||grid.solid(n1)) continue;
        if (!grid.fluid(n0)&&!grid.fluid(n1)) continue;
        grid.uFaceVel(i,j,k) -= coeff*(grid.pressure(n1)-grid.pressure(n0));
    }

    // v-faces: between node (i,j,k) and (i,j+1,k)
    for (int k=0;k<grid.nVz();++k)
    for (int j=0;j<grid.nVy();++j)
    for (int i=0;i<grid.nVx();++i) {
        int n0=grid.idx(i,j,  k);
        int n1=grid.idx(i,j+1,k);
        if (grid.solid(n0)||grid.solid(n1)) continue;
        if (!grid.fluid(n0)&&!grid.fluid(n1)) continue;
        grid.vFaceVel(i,j,k) -= coeff*(grid.pressure(n1)-grid.pressure(n0));
    }

    // w-faces: between node (i,j,k) and (i,j,k+1)
    for (int k=0;k<grid.nWz();++k)
    for (int j=0;j<grid.nWy();++j)
    for (int i=0;i<grid.nWx();++i) {
        int n0=grid.idx(i,j,k  );
        int n1=grid.idx(i,j,k+1);
        if (grid.solid(n0)||grid.solid(n1)) continue;
        if (!grid.fluid(n0)&&!grid.fluid(n1)) continue;
        grid.wFaceVel(i,j,k) -= coeff*(grid.pressure(n1)-grid.pressure(n0));
    }

    // Re-enforce boundary faces after gradient
    grid.enforceBoundaryFaces();
}

} // namespace apic
