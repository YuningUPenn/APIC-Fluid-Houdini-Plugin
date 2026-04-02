#pragma once
#include "apic_types.h"
#include <functional>

namespace apic {

// -------------------------------------------------------
// ApicGrid -- MAC staggered grid
//
// Pressure and solid/fluid markers live at nodes (i,j,k).
// Velocity lives on staggered faces:
//   u-faces: (nx-1)*ny*nz, between nodes (i,j,k) and (i+1,j,k)
//             world pos = origin + ((i+0.5)*dx, j*dx, k*dx)
//   v-faces: nx*(ny-1)*nz, between nodes (i,j,k) and (i,j+1,k)
//             world pos = origin + (i*dx, (j+0.5)*dx, k*dx)
//   w-faces: nx*ny*(nz-1), between nodes (i,j,k) and (i,j,k+1)
//             world pos = origin + (i*dx, j*dx, (k+0.5)*dx)
//
// Pressure Poisson equation (standard MAC):
//   div_face(v*) = (u[i,j,k]-u[i-1,j,k] + v[i,j,k]-v[i,j-1,k]
//                 + w[i,j,k]-w[i,j,k-1]) / dx
//   -Lap(p) = (rho/dt)*div    =>  A*p = (rho/dt)*div  (A pos-def)
//   u[i,j,k] -= (dt/rho)*(p[i+1,j,k]-p[i,j,k])/dx
// -------------------------------------------------------
class ApicGrid {
public:
    ApicGrid() = default;
    void init(int nx, int ny, int nz, Scalar dx,
              const Vec3& origin = Vec3::Zero());

    // ---- Node accessors (pressure, solid, fluid) ----
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    Scalar dx() const { return dx_; }
    const Vec3& origin() const { return origin_; }
    size_t numNodes() const { return static_cast<size_t>(nx_*ny_*nz_); }

    int idx(int i, int j, int k) const { return i + nx_*(j + ny_*k); }
    Vec3i idx3(int flat) const;
    Vec3 nodePos(int i, int j, int k) const {
        return origin_ + dx_*vec3((Scalar)i,(Scalar)j,(Scalar)k);
    }
    Vec3 nodePos(int flat) const;
    bool inside(int i, int j, int k) const {
        return i>=0&&i<nx_&&j>=0&&j<ny_&&k>=0&&k<nz_;
    }

    Scalar& mass(int flat)       { return mass_[flat]; }
    Scalar  mass(int flat) const { return mass_[flat]; }
    Scalar& mass(int i,int j,int k)      { return mass_[idx(i,j,k)]; }
    Scalar  mass(int i,int j,int k) const{ return mass_[idx(i,j,k)]; }

    Scalar& pressure(int flat)       { return pressure_[flat]; }
    Scalar  pressure(int flat) const { return pressure_[flat]; }
    Scalar& pressure(int i,int j,int k)      { return pressure_[idx(i,j,k)]; }
    Scalar  pressure(int i,int j,int k) const{ return pressure_[idx(i,j,k)]; }

    Scalar& divergence(int flat)       { return divergence_[flat]; }
    Scalar  divergence(int flat) const { return divergence_[flat]; }

    uint8_t& solid(int flat)       { return solid_[flat]; }
    bool     solid(int flat) const { return solid_[flat]!=0; }
    uint8_t& solid(int i,int j,int k)       { return solid_[idx(i,j,k)]; }
    bool     solid(int i,int j,int k) const { return solid_[idx(i,j,k)]!=0; }

    uint8_t& fluid(int flat)       { return fluid_[flat]; }
    bool     fluid(int flat) const { return fluid_[flat]!=0; }
    uint8_t& fluid(int i,int j,int k)       { return fluid_[idx(i,j,k)]; }
    bool     fluid(int i,int j,int k) const { return fluid_[idx(i,j,k)]!=0; }
    void     clearFluidMarkers()   { std::fill(fluid_.begin(),fluid_.end(),0); }

    // ---- MAC face dimensions ----
    int nUx() const { return nx_-1; }
    int nUy() const { return ny_; }
    int nUz() const { return nz_; }
    int nVx() const { return nx_; }
    int nVy() const { return ny_-1; }
    int nVz() const { return nz_; }
    int nWx() const { return nx_; }
    int nWy() const { return ny_; }
    int nWz() const { return nz_-1; }

    int numUFaces() const { return nUx()*nUy()*nUz(); }
    int numVFaces() const { return nVx()*nVy()*nVz(); }
    int numWFaces() const { return nWx()*nWy()*nWz(); }

    // ---- Face index ----
    int uIdx(int i,int j,int k) const { return i + nUx()*(j + nUy()*k); }
    int vIdx(int i,int j,int k) const { return i + nVx()*(j + nVy()*k); }
    int wIdx(int i,int j,int k) const { return i + nWx()*(j + nWy()*k); }

    bool uInside(int i,int j,int k) const {
        return i>=0&&i<nUx()&&j>=0&&j<nUy()&&k>=0&&k<nUz();
    }
    bool vInside(int i,int j,int k) const {
        return i>=0&&i<nVx()&&j>=0&&j<nVy()&&k>=0&&k<nVz();
    }
    bool wInside(int i,int j,int k) const {
        return i>=0&&i<nWx()&&j>=0&&j<nWy()&&k>=0&&k<nWz();
    }

    // ---- Face velocity accessors ----
    Scalar& uFaceVel(int i,int j,int k)       { return uFaceVel_[uIdx(i,j,k)]; }
    Scalar  uFaceVel(int i,int j,int k) const { return uFaceVel_[uIdx(i,j,k)]; }
    Scalar& vFaceVel(int i,int j,int k)       { return vFaceVel_[vIdx(i,j,k)]; }
    Scalar  vFaceVel(int i,int j,int k) const { return vFaceVel_[vIdx(i,j,k)]; }
    Scalar& wFaceVel(int i,int j,int k)       { return wFaceVel_[wIdx(i,j,k)]; }
    Scalar  wFaceVel(int i,int j,int k) const { return wFaceVel_[wIdx(i,j,k)]; }

    Scalar& uFaceMass(int i,int j,int k)       { return uFaceMass_[uIdx(i,j,k)]; }
    Scalar  uFaceMass(int i,int j,int k) const { return uFaceMass_[uIdx(i,j,k)]; }
    Scalar& vFaceMass(int i,int j,int k)       { return vFaceMass_[vIdx(i,j,k)]; }
    Scalar  vFaceMass(int i,int j,int k) const { return vFaceMass_[vIdx(i,j,k)]; }
    Scalar& wFaceMass(int i,int j,int k)       { return wFaceMass_[wIdx(i,j,k)]; }
    Scalar  wFaceMass(int i,int j,int k) const { return wFaceMass_[wIdx(i,j,k)]; }

    // FLIP snapshot accessors
    Scalar uFaceOld(int i,int j,int k) const { return uFaceOld_[uIdx(i,j,k)]; }
    Scalar vFaceOld(int i,int j,int k) const { return vFaceOld_[vIdx(i,j,k)]; }
    Scalar wFaceOld(int i,int j,int k) const { return wFaceOld_[wIdx(i,j,k)]; }

    // Raw face array access (for pressure solver and P2G)
    const std::vector<Scalar>& uFaceVelArr()  const { return uFaceVel_; }
    const std::vector<Scalar>& vFaceVelArr()  const { return vFaceVel_; }
    const std::vector<Scalar>& wFaceVelArr()  const { return wFaceVel_; }
    std::vector<Scalar>& uFaceVelArr()  { return uFaceVel_; }
    std::vector<Scalar>& vFaceVelArr()  { return vFaceVel_; }
    std::vector<Scalar>& wFaceVelArr()  { return wFaceVel_; }
    std::vector<Scalar>& uFaceMassArr() { return uFaceMass_; }
    std::vector<Scalar>& vFaceMassArr() { return vFaceMass_; }
    std::vector<Scalar>& wFaceMassArr() { return wFaceMass_; }
    const std::vector<Scalar>& uFaceOldArr() const { return uFaceOld_; }
    const std::vector<Scalar>& vFaceOldArr() const { return vFaceOld_; }
    const std::vector<Scalar>& wFaceOldArr() const { return wFaceOld_; }

    // ---- Lifecycle ----
    void clear();
    void markBoundary();

    // Normalize face momentum -> velocity, mark fluid cells, enforce BC
    void normaliseFaceMomentum();
    // Enforce zero normal velocity at boundary faces
    void enforceBoundaryFaces();
    // Snapshot face velocities for FLIP
    void snapshotVelocity();

    // ---- Kernel ----
    static Scalar kernel(Scalar x);
    static Scalar kernelGrad(Scalar x);

    // ---- Weight entries (shared for node and face) ----
    struct WeightEntry {
        int    flatIdx;
        Scalar weight;
        Vec3   gradient;   // world-space gradient of weight
        Vec3   offset;     // face_pos - particle_pos (world space)
    };

    // Node-based weights (kept for compatibility)
    void gatherWeights(const Vec3& pos,
                       std::vector<WeightEntry>& entries) const;

    // Face-based weights for MAC transfers
    // axis=0: u-faces, axis=1: v-faces, axis=2: w-faces
    void gatherWeightsFace(const Vec3& pos, int axis,
                           std::vector<WeightEntry>& entries) const;

private:
    int    nx_=0, ny_=0, nz_=0;
    Scalar dx_=1.f;
    Vec3   origin_=Vec3::Zero();

    // Node arrays (pressure, markers)
    std::vector<Scalar>  mass_;        // node mass (for fluid marking)
    std::vector<Scalar>  pressure_;
    std::vector<Scalar>  divergence_;
    std::vector<uint8_t> solid_;
    std::vector<uint8_t> fluid_;

    // MAC face velocity/mass arrays
    std::vector<Scalar> uFaceVel_;
    std::vector<Scalar> vFaceVel_;
    std::vector<Scalar> wFaceVel_;
    std::vector<Scalar> uFaceMass_;
    std::vector<Scalar> vFaceMass_;
    std::vector<Scalar> wFaceMass_;

    // FLIP face velocity snapshots
    std::vector<Scalar> uFaceOld_;
    std::vector<Scalar> vFaceOld_;
    std::vector<Scalar> wFaceOld_;
};

} // namespace apic
