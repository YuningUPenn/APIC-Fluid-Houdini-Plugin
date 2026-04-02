// test_angular_momentum.cpp  (MAC grid version)
//
// Tests angular momentum preservation for APIC vs PIC
// using MAC face-based transfers.

#include "ApicSolver.h"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace apic;
static constexpr float PI = 3.14159265f;

int main() {
    int failures = 0;

    SimParams params;
    params.gridResX = params.gridResY = params.gridResZ = 32;
    params.gridSpacing = 1.0f / 32.0f;
    params.dt = 1.0f / 24.0f;
    params.substeps = 1;
    params.gravity = Vec3::Zero();
    params.viscosity = 0.f;  // no damping for conservation test

    const Scalar dx = params.gridSpacing;
    const Scalar cx = 0.5f, cy = 0.5f;
    const Scalar r = 0.25f;
    const int    N = 64;
    const Scalar omega = 2.0f * PI;

    auto runTest = [&](TransferMethod method, const char* name) -> Scalar {
        params.method = method;
        ApicSolver solver(params);
        ApicParticleSet& P = solver.particles();
        ApicGrid& G = solver.grid();

        // Init: rigid-body rotation in XY plane
        for (int n = 0; n < N; ++n) {
            Scalar theta = 2.0f * PI * n / N;
            Vec3 pos = vec3(cx + r*std::cos(theta), cy + r*std::sin(theta), 0.5f);
            Vec3 vel = vec3(-omega*r*std::sin(theta), omega*r*std::cos(theta), 0.f);
            P.addParticle(pos, vel, 1.f);
        }

        std::vector<ApicGrid::WeightEntry> wts;

        auto doRoundTrip = [&]() {
            // P2G (MAC)
            G.clear();
            for (size_t p = 0; p < P.size(); ++p) {
                const Vec3&  xp = P.position(p);
                const Vec3&  vp = P.velocity(p);
                const Scalar mp = P.mass(p);
                const Mat3&  Cp = P.affineB(p);

                for (int axis = 0; axis < 3; ++axis) {
                    G.gatherWeightsFace(xp, axis, wts);
                    auto& faceVel  = (axis==0)?G.uFaceVelArr():(axis==1)?G.vFaceVelArr():G.wFaceVelArr();
                    auto& faceMass = (axis==0)?G.uFaceMassArr():(axis==1)?G.vFaceMassArr():G.wFaceMassArr();
                    for (const auto& we : wts) {
                        Scalar affine = (method == TransferMethod::APIC)
                                        ? Cp.col(axis).dot(we.offset) : 0.f;
                        faceMass[we.flatIdx] += we.weight * mp;
                        faceVel [we.flatIdx] += we.weight * mp * (vp(axis) + affine);
                    }
                }
            }
            G.normaliseFaceMomentum();

            // G2P (MAC)
            for (size_t p = 0; p < P.size(); ++p) {
                const Vec3& xp = P.position(p);
                Vec3 vNew = Vec3::Zero();
                Mat3 Cnew = Mat3::Zero();

                for (int axis = 0; axis < 3; ++axis) {
                    G.gatherWeightsFace(xp, axis, wts);
                    const auto& faceVel = (axis==0)?G.uFaceVelArr():(axis==1)?G.vFaceVelArr():G.wFaceVelArr();
                    for (const auto& we : wts) {
                        Scalar fv = faceVel[we.flatIdx];
                        vNew(axis) += we.weight * fv;
                        if (method == TransferMethod::APIC)
                            Cnew.col(axis) += we.gradient * fv;
                    }
                }
                P.velocity(p) = vNew;
                if (method == TransferMethod::APIC)
                    P.affineB(p) = Cnew;
                else
                    P.affineB(p).setZero();
            }
        };

        // Angular momentum (z-component)
        auto angMom = [&]() {
            Scalar L = 0.f;
            for (size_t p = 0; p < P.size(); ++p) {
                Vec3 rv = P.position(p) - vec3(cx, cy, 0.5f);
                L += P.mass(p)*(rv.x()*P.velocity(p).y() - rv.y()*P.velocity(p).x());
            }
            return L;
        };

        Scalar L0 = angMom();

        const int steps = 20;
        for (int i = 0; i < steps; ++i)
            doRoundTrip();

        Scalar Lf = angMom();
        Scalar relError = std::abs(Lf - L0) / (std::abs(L0) + 1e-8f);
        printf("[%s] L0=%.5f  Lf=%.5f  rel_error=%.5e\n", name, L0, Lf, relError);
        return relError;
    };

    Scalar apicErr = runTest(TransferMethod::APIC, "APIC");
    Scalar picErr  = runTest(TransferMethod::PIC,  "PIC ");

    printf("\n=== Angular Momentum Preservation Test ===\n");
    if (apicErr > 1e-2f) {
        printf("FAIL: APIC angular momentum error too large (%.5e)\n", apicErr);
        ++failures;
    } else {
        printf("PASS: APIC preserves angular momentum (error=%.5e)\n", apicErr);
    }

    if (apicErr < picErr)
        printf("PASS: APIC (%.5e) better than PIC (%.5e)\n", apicErr, picErr);
    else
        printf("WARNING: APIC not better than PIC\n");

    return failures > 0 ? 1 : 0;
}
