#pragma once
#include "apic_types.h"
#include "ApicGrid.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

namespace apic {

// -------------------------------------------------------
// PressureSolver
//
// Enforces incompressibility by solving the Poisson equation:
//
//   Laplacian(p) = divergence(v*) / dt
//
// Then corrects velocities:
//   v = v* - dt * grad(p)
//
// Uses Eigen's ConjugateGradient (PCG) solver with an
// Incomplete Cholesky preconditioner for the sparse system.
//
// A simpler Jacobi fallback is also provided for debugging.
// -------------------------------------------------------
class PressureSolver {
public:
    enum class Method { PCG, Jacobi };

    PressureSolver() = default;

    // Configure
    void setMethod(Method m) { method_ = m; }
    void setMaxIterations(int n) { maxIter_ = n; }
    void setTolerance(Scalar tol) { tol_ = tol; }

    // Main entry point.
    // Reads velocity from grid, builds the Poisson system,
    // solves for pressure, writes corrected velocity back.
    // Returns residual after solve (for convergence monitoring).
    Scalar solve(ApicGrid& grid, const SimParams& params);

    // Diagnostics: residual from last solve
    Scalar lastResidual() const { return lastResidual_; }
    int    lastIterations() const { return lastIter_; }

private:
    Method  method_  = Method::PCG;
    int     maxIter_ = 200;
    Scalar  tol_     = 1e-5f;
    Scalar  lastResidual_ = 0.f;
    int     lastIter_     = 0;

    // Build the flat index mapping for fluid (non-solid) nodes
    void buildFluidIndex(const ApicGrid& grid,
                         std::vector<int>& fluidToFlat,
                         std::vector<int>& flatToFluid) const;

    // Compute velocity divergence at each grid node
    void computeDivergence(ApicGrid& grid, const SimParams& params) const;

    // Apply pressure gradient correction to velocities
    void applyPressureGradient(ApicGrid& grid, const SimParams& params) const;

    // PCG solve
    Scalar solvePCG(ApicGrid& grid, const SimParams& params,
                    const std::vector<int>& fluidToFlat,
                    const std::vector<int>& flatToFluid);

    // Jacobi solve (simpler, slower -- useful for unit tests)
    Scalar solveJacobi(ApicGrid& grid, const SimParams& params,
                       const std::vector<int>& fluidToFlat,
                       const std::vector<int>& flatToFluid);
};

} // namespace apic
