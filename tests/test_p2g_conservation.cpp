// test_p2g_conservation.cpp  (MAC grid version)
//
// Tests: after a MAC P2G transfer, total mass and linear momentum
// on the faces must equal the total mass and momentum of particles.

#include "ApicSolver.h"
#include <cstdio>
#include <cmath>

using namespace apic;

static bool nearEqual(Scalar a, Scalar b, Scalar tol = 1e-3f) {
    return std::abs(a - b) < tol * (std::abs(a) + std::abs(b) + 1.f);
}

int main() {
    int failures = 0;

    SimParams params;
    params.gridResX = params.gridResY = params.gridResZ = 16;
    params.gridSpacing  = 1.0f / 16.0f;
    params.dt           = 1.0f / 24.0f;
    params.substeps     = 1;
    params.method       = TransferMethod::APIC;
    params.gravity      = Vec3::Zero();

    ApicSolver solver(params);
    ApicParticleSet& particles = solver.particles();
    ApicGrid& grid = solver.grid();

    const Scalar dx = params.gridSpacing;
    const Scalar mass_per_particle = 1.0f;

    // Place particles with known velocities
    for (int k = 4; k <= 11; ++k)
    for (int j = 4; j <= 11; ++j)
    for (int i = 4; i <= 11; ++i) {
        Vec3 pos = vec3(i + 0.5f, j + 0.5f, k + 0.5f) * dx;
        Vec3 vel = vec3(
            0.1f * std::sin((float)(i + j)),
            0.1f * std::cos((float)(j + k)),
            0.1f * std::sin((float)(i * k))
        );
        particles.addParticle(pos, vel, mass_per_particle);
    }

    // Expected totals from particles
    Scalar pMassTotal = 0.f;
    Vec3   pMomTotal  = Vec3::Zero();
    for (size_t p = 0; p < particles.size(); ++p) {
        pMassTotal += particles.mass(p);
        pMomTotal  += particles.mass(p) * particles.velocity(p);
    }

    // Run MAC P2G (replicate what ApicSolver::p2gAPIC does)
    grid.clear();
    {
        std::vector<ApicGrid::WeightEntry> weights;
        for (size_t p = 0; p < particles.size(); ++p) {
            const Vec3&  xp = particles.position(p);
            const Vec3&  vp = particles.velocity(p);
            const Scalar mp = particles.mass(p);
            const Mat3&  Cp = particles.affineB(p);  // zero initially

            for (int axis = 0; axis < 3; ++axis) {
                grid.gatherWeightsFace(xp, axis, weights);
                auto& faceVel  = (axis==0)?grid.uFaceVelArr():(axis==1)?grid.vFaceVelArr():grid.wFaceVelArr();
                auto& faceMass = (axis==0)?grid.uFaceMassArr():(axis==1)?grid.vFaceMassArr():grid.wFaceMassArr();
                for (const auto& we : weights) {
                    Scalar affineContrib = Cp.col(axis).dot(we.offset);
                    Scalar vel_a = vp(axis) + affineContrib;
                    faceMass[we.flatIdx] += we.weight * mp;
                    faceVel [we.flatIdx] += we.weight * mp * vel_a;
                }
            }
        }
    }

    // Sum face mass and momentum (before normalisation, faceVel stores momentum)
    Scalar gMassTotal = 0.f;
    Vec3   gMomTotal  = Vec3::Zero();

    // u-faces contribute to x-momentum
    for (int k=0;k<grid.nUz();++k)
    for (int j=0;j<grid.nUy();++j)
    for (int i=0;i<grid.nUx();++i) {
        gMassTotal  += grid.uFaceMass(i,j,k);
        gMomTotal.x() += grid.uFaceVel(i,j,k);   // momentum
    }
    // v-faces contribute to y-momentum
    for (int k=0;k<grid.nVz();++k)
    for (int j=0;j<grid.nVy();++j)
    for (int i=0;i<grid.nVx();++i) {
        // mass is summed per axis; divide by 3 at the end to avoid triple counting
        gMomTotal.y() += grid.vFaceVel(i,j,k);
    }
    // w-faces contribute to z-momentum
    for (int k=0;k<grid.nWz();++k)
    for (int j=0;j<grid.nWy();++j)
    for (int i=0;i<grid.nWx();++i) {
        gMomTotal.z() += grid.wFaceVel(i,j,k);
    }
    // Mass was summed only over u-faces above; sum v and w too
    for (int k=0;k<grid.nVz();++k)
    for (int j=0;j<grid.nVy();++j)
    for (int i=0;i<grid.nVx();++i)
        gMassTotal += grid.vFaceMass(i,j,k);
    for (int k=0;k<grid.nWz();++k)
    for (int j=0;j<grid.nWy();++j)
    for (int i=0;i<grid.nWx();++i)
        gMassTotal += grid.wFaceMass(i,j,k);
    // Each particle mass was scattered to 3 axes -> divide by 3
    gMassTotal /= 3.f;

    printf("=== P2G Conservation Test ===\n");
    printf("Particle mass:  %.6f\n", pMassTotal);
    printf("Grid mass:      %.6f\n", gMassTotal);
    printf("Particle mom:   [%.4f, %.4f, %.4f]\n",
           pMomTotal.x(), pMomTotal.y(), pMomTotal.z());
    printf("Grid mom:       [%.4f, %.4f, %.4f]\n",
           gMomTotal.x(), gMomTotal.y(), gMomTotal.z());

    if (!nearEqual(pMassTotal, gMassTotal)) {
        printf("FAIL: mass not conserved (diff = %.6f)\n",
               std::abs(pMassTotal - gMassTotal));
        ++failures;
    } else {
        printf("PASS: mass conserved\n");
    }

    bool momOk = true;
    for (int axis = 0; axis < 3; ++axis) {
        if (!nearEqual(pMomTotal(axis), gMomTotal(axis))) {
            printf("FAIL: momentum[%d] not conserved (diff = %.6f)\n",
                   axis, std::abs(pMomTotal(axis) - gMomTotal(axis)));
            ++failures;
            momOk = false;
        }
    }
    if (momOk) printf("PASS: linear momentum conserved\n");

    return failures > 0 ? 1 : 0;
}
