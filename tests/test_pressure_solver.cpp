// test_pressure_solver.cpp  (MAC grid version)
//
// Sets up a sinusoidal face velocity field on a MAC grid,
// runs the pressure solver, and verifies divergence is reduced.

#include "ApicGrid.h"
#include "PressureSolver.h"
#include "apic_types.h"
#include <cstdio>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace apic;

// MAC divergence: backward diff of face velocities
static Scalar computeMaxDiv(const ApicGrid& grid) {
    const int nx=grid.nx(), ny=grid.ny(), nz=grid.nz();
    const Scalar inv_dx = 1.f / grid.dx();
    Scalar maxDiv = 0.f;
    for (int k=1;k<nz-1;++k)
    for (int j=1;j<ny-1;++j)
    for (int i=1;i<nx-1;++i) {
        int flat = grid.idx(i,j,k);
        if (!grid.fluid(flat)||grid.solid(flat)) continue;

        Scalar u_r = grid.uInside(i,  j,k) ? grid.uFaceVel(i,  j,k) : 0.f;
        Scalar u_l = grid.uInside(i-1,j,k) ? grid.uFaceVel(i-1,j,k) : 0.f;
        Scalar v_t = grid.vInside(i,j,  k) ? grid.vFaceVel(i,j,  k) : 0.f;
        Scalar v_b = grid.vInside(i,j-1,k) ? grid.vFaceVel(i,j-1,k) : 0.f;
        Scalar w_f = grid.wInside(i,j,k  ) ? grid.wFaceVel(i,j,k  ) : 0.f;
        Scalar w_k = grid.wInside(i,j,k-1) ? grid.wFaceVel(i,j,k-1) : 0.f;

        Scalar d = (u_r-u_l + v_t-v_b + w_f-w_k) * inv_dx;
        maxDiv = std::max(maxDiv, std::abs(d));
    }
    return maxDiv;
}

// Fill face velocities with a sinusoidal divergent field.
// u(i,j,k) = A * sin(pi * (i+0.5) / nx)  etc.
// This has non-zero divergence which the solver should remove.
static void fillFaceVelocity(ApicGrid& grid) {
    const int nx=grid.nx(), ny=grid.ny(), nz=grid.nz();
    const Scalar A = 0.5f;

    // u-faces (x-velocity at (i+0.5, j, k))
    for (int k=0;k<grid.nUz();++k)
    for (int j=0;j<grid.nUy();++j)
    for (int i=0;i<grid.nUx();++i)
        grid.uFaceVel(i,j,k) = A*(Scalar)std::sin(M_PI*(i+0.5f)/nx);

    // v-faces (y-velocity at (i, j+0.5, k))
    for (int k=0;k<grid.nVz();++k)
    for (int j=0;j<grid.nVy();++j)
    for (int i=0;i<grid.nVx();++i)
        grid.vFaceVel(i,j,k) = A*(Scalar)std::sin(M_PI*(j+0.5f)/ny);

    // w-faces (z-velocity at (i, j, k+0.5))
    for (int k=0;k<grid.nWz();++k)
    for (int j=0;j<grid.nWy();++j)
    for (int i=0;i<grid.nWx();++i)
        grid.wFaceVel(i,j,k) = A*(Scalar)std::sin(M_PI*(k+0.5f)/nz);

    // Enforce boundary: zero normal faces at walls
    grid.enforceBoundaryFaces();

    // Mark all interior cells as fluid
    for (int k=1;k<nz-1;++k)
    for (int j=1;j<ny-1;++j)
    for (int i=1;i<nx-1;++i)
        if (!grid.solid(i,j,k)) grid.fluid(i,j,k) = 1;
}

int main() {
    int failures = 0;

    SimParams params;
    params.gridResX = params.gridResY = params.gridResZ = 16;
    params.gridSpacing = 1.0f;
    params.density = 1.0f;
    params.dt = 1.0f;

    printf("params: density=%.1f  dt=%.4f  dx=%.1f  scale=%.4f\n",
           params.density, params.dt, params.gridSpacing,
           params.density/params.dt);

    // --- PCG test ---
    ApicGrid grid;
    grid.init(params.gridResX, params.gridResY, params.gridResZ, params.gridSpacing);
    fillFaceVelocity(grid);

    Scalar divBefore = computeMaxDiv(grid);
    printf("Max divergence before pressure solve: %.6f\n", divBefore);

    PressureSolver solver;
    solver.setMethod(PressureSolver::Method::PCG);
    solver.setMaxIterations(400);
    solver.setTolerance(1e-6f);
    Scalar residual = solver.solve(grid, params);

    Scalar divAfter = computeMaxDiv(grid);
    printf("Max divergence after  pressure solve: %.6f\n", divAfter);
    printf("Solver residual: %.6e  iterations: %d\n", residual, solver.lastIterations());

    printf("\n=== Pressure Solver Convergence Test ===\n");
    if (divAfter > 0.5f * divBefore && divBefore > 1e-6f) {
        printf("FAIL: divergence not sufficiently reduced (before=%.4f after=%.4f)\n",
               divBefore, divAfter);
        ++failures;
    } else {
        printf("PASS: divergence reduced from %.4f to %.4f\n", divBefore, divAfter);
    }
    if (residual > 1e-3f) {
        printf("FAIL: solver residual too large (%.6e)\n", residual);
        ++failures;
    } else {
        printf("PASS: solver converged (residual=%.6e)\n", residual);
    }

    // --- Jacobi test ---
    grid.init(params.gridResX, params.gridResY, params.gridResZ, params.gridSpacing);
    fillFaceVelocity(grid);
    PressureSolver jacobi;
    jacobi.setMethod(PressureSolver::Method::Jacobi);
    jacobi.setMaxIterations(1000);
    jacobi.setTolerance(1e-5f);
    Scalar jRes = jacobi.solve(grid, params);
    printf("\n--- Jacobi fallback ---\n");
    printf("Max divergence after Jacobi: %.6f  (residual=%.6e)\n",
           computeMaxDiv(grid), jRes);

    return failures > 0 ? 1 : 0;
}
