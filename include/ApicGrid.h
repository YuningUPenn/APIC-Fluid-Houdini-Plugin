#pragma once
#include "apic_types.h"
#include <functional>

namespace apic {

// -------------------------------------------------------
// ApicGrid
//
// Uniform collocated Cartesian grid.
// Stores per-node: mass, velocity, force.
// Also holds temporary pressure and divergence buffers
// needed during the projection step.
//
// Indexing:  flat index = i + nx*(j + ny*k)
//            node world position = origin + dx * Vec3(i, j, k)
//
// NOTE: The design doc specifies a collocated grid (all
// quantities at cell centres/nodes) rather than a staggered
// MAC grid, for simplicity.  This simplifies P2G/G2P but
// means the pressure projection uses a collocated Poisson
// solve (standard central differences).
// -------------------------------------------------------
class ApicGrid {
public:
    // ---- Construction ----
    ApicGrid() = default;
    void init(int nx, int ny, int nz, Scalar dx, const Vec3& origin = Vec3::Zero());

    // ---- Accessors ----
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    Scalar dx() const { return dx_; }
    const Vec3& origin() const { return origin_; }
    size_t numNodes() const { return static_cast<size_t>(nx_ * ny_ * nz_); }

    // Flat index from 3D index
    int idx(int i, int j, int k) const { return i + nx_ * (j + ny_ * k); }
    // 3D index from flat
    Vec3i idx3(int flat) const;
    // World position of node
    Vec3 nodePos(int i, int j, int k) const {
        Vec3 p; p.x() = i; p.y() = j; p.z() = k;
        return origin_ + dx_ * p;
    }
    Vec3 nodePos(int flat) const;

    // Is a node inside the domain?
    bool inside(int i, int j, int k) const {
        return i>=0 && i<nx_ && j>=0 && j<ny_ && k>=0 && k<nz_;
    }

    // ---- Grid data buffers ----
    Scalar& mass(int flat)     { return mass_[flat]; }
    Vec3&   velocity(int flat) { return velocity_[flat]; }
    Vec3&   force(int flat)    { return force_[flat]; }
    Scalar& pressure(int flat) { return pressure_[flat]; }
    Scalar& divergence(int flat){ return divergence_[flat]; }

    // Convenience 3D overloads
    Scalar& mass(int i,int j,int k)      { return mass_[idx(i,j,k)]; }
    Vec3&   velocity(int i,int j,int k)  { return velocity_[idx(i,j,k)]; }
    Vec3&   force(int i,int j,int k)     { return force_[idx(i,j,k)]; }
    Scalar& pressure(int i,int j,int k)  { return pressure_[idx(i,j,k)]; }

    // Const overloads (flat index)
    Scalar mass(int flat) const       { return mass_[flat]; }
    Vec3   velocity(int flat) const   { return velocity_[flat]; }
    Scalar pressure(int flat) const   { return pressure_[flat]; }
    Scalar divergence(int flat) const { return divergence_[flat]; }
    // Const overloads (3D index)
    Scalar mass(int i,int j,int k)     const { return mass_[idx(i,j,k)]; }
    Vec3   velocity(int i,int j,int k) const { return velocity_[idx(i,j,k)]; }
    Scalar pressure(int i,int j,int k) const { return pressure_[idx(i,j,k)]; }

    // ---- Lifecycle ----
    // Zero all grid buffers (mass, velocity, force, pressure, divergence)
    void clear();

    // Normalise accumulated momentum by mass to get velocity
    // (call after P2G accumulation)
    void normaliseMomentum();

    // ---- Interpolation kernel ----
    // Quadratic B-spline N(x) used for weights w_ip = N((x_p - x_i)/dx)
    // Returns weight and optionally its gradient.
    static Scalar kernel(Scalar x);
    static Scalar kernelGrad(Scalar x); // dN/dx

    // Compute interpolation weights and indices for a world-space particle position.
    // Fills 'nodes' with up to 27 (3x3x3) node indices and corresponding weights.
    struct WeightEntry {
        int   flatIdx;
        Scalar weight;
        Vec3  gradient;   // grad of weight in world space (= dN/dx / dx)
        Vec3  offset;     // x_i - x_p (world space)
    };
    void gatherWeights(const Vec3& pos,
                       std::vector<WeightEntry>& entries) const;

    // ---- SDF boundary ----
    // Marks nodes outside the fluid domain as solid.
    // solid_[i] = true means the node is a solid boundary cell.
    uint8_t& solid(int flat)       { return solid_[flat]; }
    bool     solid(int flat) const { return solid_[flat] != 0; }
    uint8_t& solid(int i,int j,int k)       { return solid_[idx(i,j,k)]; }
    bool     solid(int i,int j,int k) const { return solid_[idx(i,j,k)] != 0; }

    // Mark domain boundary nodes as solid
    void markBoundary();

    // ---- Fluid marker ----
    // fluid_[i] = true once mass has been rasterised to that node
    uint8_t& fluid(int flat)       { return fluid_[flat]; }
    bool     fluid(int flat) const { return fluid_[flat] != 0; }
    void     clearFluidMarkers()   { std::fill(fluid_.begin(), fluid_.end(), 0); }

private:
    int nx_ = 0, ny_ = 0, nz_ = 0;
    Scalar dx_ = 1.0f;
    Vec3   origin_ = Vec3::Zero();

    std::vector<Scalar> mass_;
    std::vector<Vec3>   velocity_;
    std::vector<Vec3>   force_;
    std::vector<Scalar> pressure_;
    std::vector<Scalar> divergence_;
    std::vector<uint8_t> solid_;
    std::vector<uint8_t> fluid_;

    // Velocity snapshot before grid update (used by FLIP/Hybrid)
    std::vector<Vec3>   velocityOld_;

public:
    // Save/restore velocity snapshot for FLIP increment
    void snapshotVelocity();
    Vec3 velocityOld(int flat) const { return velocityOld_[flat]; }
};

} // namespace apic
