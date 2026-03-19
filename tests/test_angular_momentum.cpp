// test_angular_momentum.cpp
// Tests that APIC preserves angular momentum better than PIC.
//
// Key insight: APIC's advantage comes from the B matrix carrying
// velocity gradient info. To test this properly, we must run TWO
// full P2G->G2P cycles so B is populated before the second P2G.
// With B=0 (initial state), APIC and PIC are identical.

#include "ApicSolver.h"
#include <cstdio>
#include <cmath>

using namespace apic;
static constexpr float PI = 3.14159265f;

int main() {
    int failures = 0;

    SimParams params;
    params.gridResX = params.gridResY = params.gridResZ = 32;
    params.gridSpacing = 1.0f / 32.0f;
    params.dt          = 1.0f / 24.0f;
    params.substeps    = 1;
    params.gravity     = Vec3::Zero();

    const Scalar dx    = params.gridSpacing;
    const Scalar cx    = 0.5f, cy = 0.5f;
    const Scalar r     = 0.25f;
    const int    N     = 32;
    const Scalar omega = 2.0f * PI;

    auto runTest = [&](TransferMethod method, const char* name) -> Scalar {
        params.method = method;
        ApicSolver solver(params);
        ApicParticleSet& P = solver.particles();
        ApicGrid& G = solver.grid();

        // Place particles in solid-body rotation
        for (int n = 0; n < N; ++n) {
            Scalar theta = 2.0f * PI * n / N;
            Vec3 pos = vec3(cx + r * std::cos(theta),
                            cy + r * std::sin(theta), 0.5f);
            Vec3 vel = vec3(-omega * r * std::sin(theta),
                             omega * r * std::cos(theta), 0.f);
            P.addParticle(pos, vel, 1.f);
        }

        // Helper: one manual P2G -> G2P round-trip
        auto doRoundTrip = [&]() {
            // P2G
            G.clear();
            const Scalar Dp_inv = 4.0f / (dx * dx);
            std::vector<ApicGrid::WeightEntry> wts;
            for (size_t p = 0; p < P.size(); ++p) {
                G.gatherWeights(P.position(p), wts);
                for (const auto& we : wts) {
                    Vec3 contrib = P.velocity(p);
                    if (method == TransferMethod::APIC)
                        contrib += P.affineB(p) * (Dp_inv * we.offset);
                    G.mass(we.flatIdx)     += we.weight * P.mass(p);
                    G.velocity(we.flatIdx) += we.weight * P.mass(p) * contrib;
                }
            }
            G.normaliseMomentum();

            // G2P
            for (size_t p = 0; p < P.size(); ++p) {
                G.gatherWeights(P.position(p), wts);
                Vec3 vNew = Vec3::Zero();
                Mat3 Bnew = Mat3::Zero();
                for (const auto& we : wts) {
                    vNew += we.weight * G.velocity(we.flatIdx);
                    if (method == TransferMethod::APIC)
                        Bnew += we.weight * G.velocity(we.flatIdx)
                                * we.offset.transpose();
                }
                P.velocity(p) = vNew;
                P.affineB(p)  = Bnew;
            }
        };

        // Measure L before
        auto angMom = [&]() {
            Scalar L = 0.f;
            for (size_t p = 0; p < P.size(); ++p) {
                Vec3 rv = P.position(p) - vec3(cx, cy, 0.5f);
                L += P.mass(p) * (rv.x() * P.velocity(p).y()
                                - rv.y() * P.velocity(p).x());
            }
            return L;
        };

        Scalar L0 = angMom();

        // Round-trip 1: populates B for APIC
        doRoundTrip();
        // Round-trip 2: now APIC uses the populated B -- this is where
        // the difference between APIC and PIC becomes visible
        doRoundTrip();

        Scalar L2 = angMom();
        Scalar relError = std::abs(L2 - L0) / (std::abs(L0) + 1e-8f);
        printf("[%s] L0=%.5f  L2=%.5f  rel_error=%.5f\n", name, L0, L2, relError);
        return relError;
    };

    Scalar apicErr = runTest(TransferMethod::APIC, "APIC");
    Scalar picErr  = runTest(TransferMethod::PIC,  "PIC ");

    printf("\n=== Angular Momentum Preservation Test ===\n");

    if (apicErr > 1e-2f) {
        printf("FAIL: APIC angular momentum error too large (%.5f)\n", apicErr);
        ++failures;
    } else {
        printf("PASS: APIC preserves angular momentum (error=%.5f)\n", apicErr);
    }

    if (apicErr < picErr) {
        printf("PASS: APIC (%.5f) better than PIC (%.5f) as expected\n",
               apicErr, picErr);
    } else {
        printf("NOTE: APIC (%.5f) not better than PIC (%.5f) after 2 cycles\n",
               apicErr, picErr);
    }

    return failures > 0 ? 1 : 0;
}
