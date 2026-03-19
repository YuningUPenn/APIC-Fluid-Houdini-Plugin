#include "ApicGrid.h"
#include "PressureSolver.h"
#include "apic_types.h"
#include <cstdio>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace apic;

// Central-difference divergence matching PressureSolver internals.
static Scalar computeMaxDiv(const ApicGrid& grid) {
    const int nx=grid.nx(), ny=grid.ny(), nz=grid.nz();
    const Scalar inv2dx = 0.5f / grid.dx();
    Scalar maxDiv = 0.f;
    for (int k=1;k<nz-1;++k)
    for (int j=1;j<ny-1;++j)
    for (int i=1;i<nx-1;++i) {
        int flat=grid.idx(i,j,k);
        if (!grid.fluid(flat)||grid.solid(flat)) continue;
        auto vAt=[&](int ii,int jj,int kk,int ax)->Scalar{
            if(!grid.inside(ii,jj,kk)||grid.solid(grid.idx(ii,jj,kk)))return 0.f;
            return grid.velocity(ii,jj,kk)(ax);};
        Scalar d=(vAt(i+1,j,k,0)-vAt(i-1,j,k,0)
                 +vAt(i,j+1,k,1)-vAt(i,j-1,k,1)
                 +vAt(i,j,k+1,2)-vAt(i,j,k-1,2))*inv2dx;
        maxDiv=std::max(maxDiv,std::abs(d));
    }
    return maxDiv;
}

static void fillVelocity(ApicGrid& grid) {
    const int nx=grid.nx(),ny=grid.ny(),nz=grid.nz();
    const Scalar dx=grid.dx();
    const Scalar A=0.5f;
    // Sinusoidal field: zero normal velocity at all walls (natural Neumann)
    for (int k=1;k<nz-1;++k)
    for (int j=1;j<ny-1;++j)
    for (int i=1;i<nx-1;++i) {
        Vec3 v;
        v.x()=A*(Scalar)std::sin(M_PI*i/(nx-1));
        v.y()=A*(Scalar)std::sin(M_PI*j/(ny-1));
        v.z()=A*(Scalar)std::sin(M_PI*k/(nz-1));
        grid.velocity(i,j,k)=v;
        grid.mass(i,j,k)=1.f;
        grid.fluid(grid.idx(i,j,k))=1;
    }
}

int main() {
    int failures=0;

    SimParams params;
    params.gridResX=params.gridResY=params.gridResZ=16;
    params.gridSpacing=1.0f;
    params.density=1.0f;
    params.dt=1.0f;

    printf("params: density=%.1f  dt=%.4f  dx=%.1f  scale=%.4f\n",
           params.density,params.dt,params.gridSpacing,
           params.density/params.dt);

    // --- PCG test ---
    ApicGrid grid;
    grid.init(params.gridResX,params.gridResY,params.gridResZ,params.gridSpacing);
    fillVelocity(grid);

    Scalar divBefore=computeMaxDiv(grid);
    printf("Max divergence before pressure solve: %.6f\n",divBefore);

    PressureSolver solver;
    solver.setMethod(PressureSolver::Method::PCG);
    solver.setMaxIterations(400);
    solver.setTolerance(1e-6f);
    Scalar residual=solver.solve(grid,params);

    Scalar divAfter=computeMaxDiv(grid);
    printf("Max divergence after  pressure solve: %.6f\n",divAfter);
    printf("Solver residual: %.6e  iterations: %d\n",residual,solver.lastIterations());

    printf("\n=== Pressure Solver Convergence Test ===\n");
    if (divAfter > 0.5f * divBefore && divBefore > 1e-6f) {
        printf("FAIL: divergence not sufficiently reduced (before=%.4f after=%.4f)\n",
               divBefore,divAfter);
        ++failures;
    } else {
        printf("PASS: divergence reduced from %.4f to %.4f\n",divBefore,divAfter);
    }
    if (residual>1e-3f) {
        printf("FAIL: solver residual too large (%.6e)\n",residual);
        ++failures;
    } else {
        printf("PASS: solver converged (residual=%.6e)\n",residual);
    }

    // --- Jacobi test ---
    grid.init(params.gridResX,params.gridResY,params.gridResZ,params.gridSpacing);
    fillVelocity(grid);
    PressureSolver jacobi;
    jacobi.setMethod(PressureSolver::Method::Jacobi);
    jacobi.setMaxIterations(1000);
    jacobi.setTolerance(1e-5f);
    Scalar jRes=jacobi.solve(grid,params);
    printf("\n--- Jacobi fallback ---\n");
    printf("Max divergence after Jacobi: %.6f  (residual=%.6e)\n",
           computeMaxDiv(grid),jRes);

    return failures>0?1:0;
}
