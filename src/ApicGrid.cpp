#include "ApicGrid.h"
#include <cmath>
#include <algorithm>
#include <cassert>

namespace apic {

void ApicGrid::init(int nx, int ny, int nz, Scalar dx, const Vec3& origin) {
    nx_=nx; ny_=ny; nz_=nz;
    dx_=dx;
    origin_=origin;

    size_t n = static_cast<size_t>(nx*ny*nz);
    mass_.assign(n, 0.f);
    pressure_.assign(n, 0.f);
    divergence_.assign(n, 0.f);
    solid_.assign(n, 0);
    fluid_.assign(n, 0);

    // Allocate face arrays
    uFaceVel_.assign(numUFaces(), 0.f);
    vFaceVel_.assign(numVFaces(), 0.f);
    wFaceVel_.assign(numWFaces(), 0.f);
    uFaceMass_.assign(numUFaces(), 0.f);
    vFaceMass_.assign(numVFaces(), 0.f);
    wFaceMass_.assign(numWFaces(), 0.f);
    uFaceOld_.assign(numUFaces(), 0.f);
    vFaceOld_.assign(numVFaces(), 0.f);
    wFaceOld_.assign(numWFaces(), 0.f);

    markBoundary();
}

void ApicGrid::clear() {
    std::fill(mass_.begin(),       mass_.end(),       0.f);
    std::fill(pressure_.begin(),   pressure_.end(),   0.f);
    std::fill(divergence_.begin(), divergence_.end(), 0.f);
    std::fill(fluid_.begin(),      fluid_.end(),      0);

    std::fill(uFaceVel_.begin(),  uFaceVel_.end(),  0.f);
    std::fill(vFaceVel_.begin(),  vFaceVel_.end(),  0.f);
    std::fill(wFaceVel_.begin(),  wFaceVel_.end(),  0.f);
    std::fill(uFaceMass_.begin(), uFaceMass_.end(), 0.f);
    std::fill(vFaceMass_.begin(), vFaceMass_.end(), 0.f);
    std::fill(wFaceMass_.begin(), wFaceMass_.end(), 0.f);
    // solid_ is NOT cleared
}

void ApicGrid::markBoundary() {
    for (int k=0;k<nz_;++k)
    for (int j=0;j<ny_;++j)
    for (int i=0;i<nx_;++i) {
        solid(i,j,k) = (i==0||i==nx_-1||j==0||j==ny_-1||k==0||k==nz_-1) ? 1 : 0;
    }
}

Vec3i ApicGrid::idx3(int flat) const {
    int k=flat/(nx_*ny_), rem=flat%(nx_*ny_);
    return {rem%nx_, rem/nx_, k};
}

Vec3 ApicGrid::nodePos(int flat) const {
    Vec3i v=idx3(flat);
    return origin_+dx_*vec3((Scalar)v[0],(Scalar)v[1],(Scalar)v[2]);
}

// -------------------------------------------------------
// Normalise face momentum -> velocity
// Mark cells as fluid based on face mass
// Enforce boundary faces
// -------------------------------------------------------
void ApicGrid::normaliseFaceMomentum() {
    // Normalize u-faces
    for (int k=0;k<nUz();++k)
    for (int j=0;j<nUy();++j)
    for (int i=0;i<nUx();++i) {
        Scalar m = uFaceMass(i,j,k);
        if (m > 1e-10f) uFaceVel(i,j,k) /= m;
        else             uFaceVel(i,j,k) = 0.f;
    }
    // Normalize v-faces
    for (int k=0;k<nVz();++k)
    for (int j=0;j<nVy();++j)
    for (int i=0;i<nVx();++i) {
        Scalar m = vFaceMass(i,j,k);
        if (m > 1e-10f) vFaceVel(i,j,k) /= m;
        else             vFaceVel(i,j,k) = 0.f;
    }
    // Normalize w-faces
    for (int k=0;k<nWz();++k)
    for (int j=0;j<nWy();++j)
    for (int i=0;i<nWx();++i) {
        Scalar m = wFaceMass(i,j,k);
        if (m > 1e-10f) wFaceVel(i,j,k) /= m;
        else             wFaceVel(i,j,k) = 0.f;
    }

    // Mark cells as fluid: a cell is fluid if any adjacent face has mass
    clearFluidMarkers();
    for (int k=1;k<nz_-1;++k)
    for (int j=1;j<ny_-1;++j)
    for (int i=1;i<nx_-1;++i) {
        if (solid(i,j,k)) continue;
        bool fl = false;
        // Check adjacent u-faces (i-1 and i)
        if (i>0     && uFaceMass(i-1,j,k)>1e-10f) fl=true;
        if (i<nx_-1 && uFaceMass(i,  j,k)>1e-10f) fl=true;
        // Check adjacent v-faces (j-1 and j)
        if (j>0     && vFaceMass(i,j-1,k)>1e-10f) fl=true;
        if (j<ny_-1 && vFaceMass(i,j,  k)>1e-10f) fl=true;
        // Check adjacent w-faces (k-1 and k)
        if (k>0     && wFaceMass(i,j,k-1)>1e-10f) fl=true;
        if (k<nz_-1 && wFaceMass(i,j,k  )>1e-10f) fl=true;
        fluid(i,j,k) = fl ? 1 : 0;
    }

    enforceBoundaryFaces();
}

void ApicGrid::enforceBoundaryFaces() {
    // Zero out u-faces adjacent to solid walls in x
    for (int k=0;k<nUz();++k)
    for (int j=0;j<nUy();++j) {
        uFaceVel(0,       j,k) = 0.f;   // left wall (between node 0 and 1)
        uFaceVel(nUx()-1, j,k) = 0.f;   // right wall (between node nx-2 and nx-1)
    }
    // Zero out v-faces adjacent to solid walls in y
    for (int k=0;k<nVz();++k)
    for (int i=0;i<nVx();++i) {
        vFaceVel(i,0,       k) = 0.f;
        vFaceVel(i,nVy()-1, k) = 0.f;
    }
    // Zero out w-faces adjacent to solid walls in z
    for (int j=0;j<nWy();++j)
    for (int i=0;i<nWx();++i) {
        wFaceVel(i,j,0      ) = 0.f;
        wFaceVel(i,j,nWz()-1) = 0.f;
    }
}

void ApicGrid::snapshotVelocity() {
    uFaceOld_ = uFaceVel_;
    vFaceOld_ = vFaceVel_;
    wFaceOld_ = wFaceVel_;
}

// -------------------------------------------------------
// Quadratic B-spline kernel (same as before)
// -------------------------------------------------------
Scalar ApicGrid::kernel(Scalar x) {
    Scalar ax=std::abs(x);
    if (ax<0.5f) return 0.75f-ax*ax;
    if (ax<1.5f) { Scalar t=1.5f-ax; return 0.5f*t*t; }
    return 0.f;
}

Scalar ApicGrid::kernelGrad(Scalar x) {
    Scalar ax=std::abs(x);
    if (ax<0.5f) return -2.f*x;
    if (ax<1.5f) return (x>0.f?-1.f:1.f)*(1.5f-ax);
    return 0.f;
}

// -------------------------------------------------------
// gatherWeights -- node-based (kept for compatibility)
// -------------------------------------------------------
void ApicGrid::gatherWeights(const Vec3& pos,
                              std::vector<WeightEntry>& entries) const {
    entries.clear();
    entries.reserve(27);

    Vec3 s=(pos-origin_)*(1.f/dx_);
    int bx=(int)std::floor(s.x()-0.5f)+1;
    int by=(int)std::floor(s.y()-0.5f)+1;
    int bz=(int)std::floor(s.z()-0.5f)+1;

    for (int dk=-1;dk<=1;++dk)
    for (int dj=-1;dj<=1;++dj)
    for (int di=-1;di<=1;++di) {
        int ni=bx+di, nj=by+dj, nk=bz+dk;
        if (!inside(ni,nj,nk)) continue;
        Scalar fx=s.x()-ni, fy=s.y()-nj, fz=s.z()-nk;
        Scalar wx=kernel(fx),wy=kernel(fy),wz=kernel(fz);
        Scalar w=wx*wy*wz;
        if (w<1e-12f) continue;
        Vec3 grad;
        grad.x()=(kernelGrad(fx)/dx_)*wy*wz;
        grad.y()=wx*(kernelGrad(fy)/dx_)*wz;
        grad.z()=wx*wy*(kernelGrad(fz)/dx_);
        Vec3 xi=nodePos(ni,nj,nk);
        entries.push_back({idx(ni,nj,nk), w, grad, xi-pos});
    }
}

// -------------------------------------------------------
// gatherWeightsFace -- MAC face-based weights
//
// For axis=0 (u-faces at (i+0.5, j, k) in grid units):
//   Particle coord in u-space: (sx-0.5, sy, sz)
//
// For axis=1 (v-faces at (i, j+0.5, k)):
//   Particle coord in v-space: (sx, sy-0.5, sz)
//
// For axis=2 (w-faces at (i, j, k+0.5)):
//   Particle coord in w-space: (sx, sy, sz-0.5)
// -------------------------------------------------------
void ApicGrid::gatherWeightsFace(const Vec3& pos, int axis,
                                  std::vector<WeightEntry>& entries) const {
    entries.clear();
    entries.reserve(27);

    Vec3 s=(pos-origin_)*(1.f/dx_);

    // Shift particle coordinate by -0.5 in the face axis direction
    Scalar sx=s.x(), sy=s.y(), sz=s.z();
    if (axis==0) sx-=0.5f;
    if (axis==1) sy-=0.5f;
    if (axis==2) sz-=0.5f;

    // Base index for 3-point stencil
    int bx=(int)std::floor(sx-0.5f)+1;
    int by=(int)std::floor(sy-0.5f)+1;
    int bz=(int)std::floor(sz-0.5f)+1;

    for (int dk=-1;dk<=1;++dk)
    for (int dj=-1;dj<=1;++dj)
    for (int di=-1;di<=1;++di) {
        int fi=bx+di, fj=by+dj, fk=bz+dk;

        // Bounds check per axis
        bool ok=false;
        if      (axis==0) ok=uInside(fi,fj,fk);
        else if (axis==1) ok=vInside(fi,fj,fk);
        else              ok=wInside(fi,fj,fk);
        if (!ok) continue;

        // Fractional distances from face to particle (in shifted coords)
        Scalar fx=sx-fi, fy=sy-fj, fz=sz-fk;

        Scalar wx=kernel(fx), wy=kernel(fy), wz=kernel(fz);
        Scalar w=wx*wy*wz;
        if (w<1e-12f) continue;

        Vec3 grad;
        grad.x()=(kernelGrad(fx)/dx_)*wy*wz;
        grad.y()=wx*(kernelGrad(fy)/dx_)*wz;
        grad.z()=wx*wy*(kernelGrad(fz)/dx_);

        // World-space face position
        Vec3 facePos=origin_;
        if (axis==0) {
            facePos.x()+=((Scalar)fi+0.5f)*dx_;
            facePos.y()+=((Scalar)fj)*dx_;
            facePos.z()+=((Scalar)fk)*dx_;
        } else if (axis==1) {
            facePos.x()+=((Scalar)fi)*dx_;
            facePos.y()+=((Scalar)fj+0.5f)*dx_;
            facePos.z()+=((Scalar)fk)*dx_;
        } else {
            facePos.x()+=((Scalar)fi)*dx_;
            facePos.y()+=((Scalar)fj)*dx_;
            facePos.z()+=((Scalar)fk+0.5f)*dx_;
        }

        int flatF=0;
        if      (axis==0) flatF=uIdx(fi,fj,fk);
        else if (axis==1) flatF=vIdx(fi,fj,fk);
        else              flatF=wIdx(fi,fj,fk);

        entries.push_back({flatF, w, grad, facePos-pos});
    }
}

} // namespace apic
