#pragma once
#include <SOP/SOP_Node.h>
#include "ApicSolver.h"
#include "CollisionHandler.h"
#include <memory>

// -------------------------------------------------------
// SOP_ApicFluidSolver
//
// The Houdini-facing SOP node that wraps the APIC solver.
// On each cook:
//   1. Reads particle geometry from input 0
//   2. Reads (optional) collision SDF geometry from input 1
//   3. Reads HDA parameters and populates SimParams
//   4. Calls ApicSolver::step() / stepFrame()
//   5. Writes updated particle positions/velocities/Bmatrix
//      back to the output GU_Detail
// -------------------------------------------------------
class SOP_ApicFluidSolver : public SOP_Node {
public:
    static OP_Node*  myConstructor(OP_Network*, const char*, OP_Operator*);
    static PRM_Template myTemplateList[];

    // Returns the operator table name used to register this SOP
    static const char* theSOPTypeName;

protected:
    SOP_ApicFluidSolver(OP_Network* net, const char* name, OP_Operator* op);
    ~SOP_ApicFluidSolver() override = default;

    // Called by Houdini when the node needs to re-cook
    OP_ERROR cookMySop(OP_Context& context) override;

    // Input labels shown in the network editor
    const char* inputLabel(unsigned idx) const;

private:
    // ---- Houdini <-> SimParams bridge ----
    apic::SimParams readParams(fpreal t) const;

    // ---- Geometry I/O ----
    // Load particles from input GU_Detail into the solver
    bool loadParticles(const GU_Detail* src);

    // Write solver particle state back to output GU_Detail
    void writeParticles(GU_Detail* dst);

    // Build collision SDF from optional input geometry (input 1)
    void buildCollisionSDF(const GU_Detail* collGeo);

    // ---- State ----
    std::unique_ptr<apic::ApicSolver>      solver_;
    std::unique_ptr<apic::CollisionHandler> collision_;

    // Track the last cooked time to detect resets
    fpreal lastCookedTime_ = -1e30;

    // Attribute handles cached across cooks
    GA_RWHandleV3  velAttrib_;
    GA_RWHandleFA  bmatAttrib_;   // APIC affine matrix stored as 9 floats

    // ---- Parameter name constants ----
    static const char* PARM_PARTICLE_SEP;
    static const char* PARM_TIME_SCALE;
    static const char* PARM_SUBSTEPS;
    static const char* PARM_TRANSFER_METHOD;
    static const char* PARM_VISCOSITY;
    static const char* PARM_SURFACE_TENSION;
    static const char* PARM_GRAVITY;
    static const char* PARM_GRID_RES;
    static const char* PARM_SHOW_PARTICLES;
    static const char* PARM_SHOW_GRID;
    static const char* PARM_SHOW_VECTORS;
    static const char* PARM_COLOR_BY_VEL;
};
