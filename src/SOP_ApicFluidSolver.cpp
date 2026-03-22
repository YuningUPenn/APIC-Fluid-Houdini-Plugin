#include "SOP_ApicFluidSolver.h"
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <UT/UT_DSOVersion.h>
#include <GA/GA_Handle.h>
#include <GA/GA_Iterator.h>
#include <GU/GU_Detail.h>

// -------------------------------------------------------
// Parameter name constants
// -------------------------------------------------------
const char* SOP_ApicFluidSolver::PARM_PARTICLE_SEP    = "particlesep";
const char* SOP_ApicFluidSolver::PARM_TIME_SCALE      = "timescale";
const char* SOP_ApicFluidSolver::PARM_SUBSTEPS        = "substeps";
const char* SOP_ApicFluidSolver::PARM_TRANSFER_METHOD = "transfermethod";
const char* SOP_ApicFluidSolver::PARM_VISCOSITY       = "viscosity";
const char* SOP_ApicFluidSolver::PARM_SURFACE_TENSION = "surfacetension";
const char* SOP_ApicFluidSolver::PARM_GRAVITY         = "gravity";
const char* SOP_ApicFluidSolver::PARM_GRID_RES        = "gridres";
const char* SOP_ApicFluidSolver::PARM_SHOW_PARTICLES  = "showparticles";
const char* SOP_ApicFluidSolver::PARM_SHOW_GRID       = "showgrid";
const char* SOP_ApicFluidSolver::PARM_SHOW_VECTORS    = "showvectors";
const char* SOP_ApicFluidSolver::PARM_COLOR_BY_VEL    = "colorbyvel";
const char* SOP_ApicFluidSolver::theSOPTypeName       = "apic_fluid_solver";

// -------------------------------------------------------
// PRM_Name definitions
// -------------------------------------------------------
static PRM_Name prmParticleSep (SOP_ApicFluidSolver::PARM_PARTICLE_SEP,  "Particle Separation");
static PRM_Name prmTimeScale   (SOP_ApicFluidSolver::PARM_TIME_SCALE,    "Time Scale");
static PRM_Name prmSubsteps    (SOP_ApicFluidSolver::PARM_SUBSTEPS,      "Substeps");
static PRM_Name prmMethod      (SOP_ApicFluidSolver::PARM_TRANSFER_METHOD,"Transfer Method");
static PRM_Name prmViscosity   (SOP_ApicFluidSolver::PARM_VISCOSITY,     "Viscosity");
static PRM_Name prmSurfTens    (SOP_ApicFluidSolver::PARM_SURFACE_TENSION,"Surface Tension");
static PRM_Name prmGravity     (SOP_ApicFluidSolver::PARM_GRAVITY,       "Gravity");
static PRM_Name prmGridRes     (SOP_ApicFluidSolver::PARM_GRID_RES,      "Grid Resolution");
static PRM_Name prmShowPart    (SOP_ApicFluidSolver::PARM_SHOW_PARTICLES,"Show Particles");
static PRM_Name prmShowGrid    (SOP_ApicFluidSolver::PARM_SHOW_GRID,     "Show Grid");
static PRM_Name prmShowVec     (SOP_ApicFluidSolver::PARM_SHOW_VECTORS,  "Show Vectors");
static PRM_Name prmColorVel    (SOP_ApicFluidSolver::PARM_COLOR_BY_VEL,  "Color by Vel");

// Default values
static PRM_Default defParticleSep(0.05);
static PRM_Default defTimeScale  (1.0);
static PRM_Default defSubsteps   (2);
static PRM_Default defViscosity  (0.01);
static PRM_Default defSurfTens   (0.0);
static PRM_Default defGravity0(0);
static PRM_Default defGravity1(-9.8);
static PRM_Default defGravity2(0);
static PRM_Default defGravity[3] = { defGravity0, defGravity1, defGravity2 };
static PRM_Default defGridRes    (64);
static PRM_Default defShowPart   (1);
static PRM_Default defShowGrid   (0);
static PRM_Default defShowVec    (0);
static PRM_Default defColorVel   (0);
static PRM_Default defMethod     (0);

// -------------------------------------------------------
// Parameter template list
// -------------------------------------------------------
PRM_Template SOP_ApicFluidSolver::myTemplateList[] = {
    PRM_Template(PRM_FLT_J, 1, &prmParticleSep, &defParticleSep),
    PRM_Template(PRM_FLT_J, 1, &prmTimeScale,   &defTimeScale),
    PRM_Template(PRM_INT_J, 1, &prmSubsteps,    &defSubsteps),
    PRM_Template(PRM_INT_J, 1, &prmGridRes,     &defGridRes),
    PRM_Template(PRM_INT_J, 1, &prmMethod,      &defMethod),
    PRM_Template(PRM_FLT_J, 1, &prmViscosity,   &defViscosity),
    PRM_Template(PRM_FLT_J, 1, &prmSurfTens,    &defSurfTens),
    PRM_Template(PRM_XYZ_J, 3, &prmGravity,      defGravity),
    PRM_Template(PRM_TOGGLE, 1, &prmShowPart,   &defShowPart),
    PRM_Template(PRM_TOGGLE, 1, &prmShowGrid,   &defShowGrid),
    PRM_Template(PRM_TOGGLE, 1, &prmShowVec,    &defShowVec),
    PRM_Template(PRM_TOGGLE, 1, &prmColorVel,   &defColorVel),
    PRM_Template()
};

// -------------------------------------------------------
// Constructor / factory
// -------------------------------------------------------
OP_Node* SOP_ApicFluidSolver::myConstructor(OP_Network* net,
                                              const char*  name,
                                              OP_Operator* op) {
    return new SOP_ApicFluidSolver(net, name, op);
}

SOP_ApicFluidSolver::SOP_ApicFluidSolver(OP_Network* net,
                                          const char*  name,
                                          OP_Operator* op)
    : SOP_Node(net, name, op)
{
    // setInputDescription does not exist in SOP_Node; input labels come from
    // the inputLabel() virtual method below.
}

// HDK signature: const char* inputLabel(unsigned) const override
const char* SOP_ApicFluidSolver::inputLabel(unsigned idx) const {
    switch (idx) {
        case 0:  return "Fluid Source";
        case 1:  return "Collision Geometry";
        default: return "Input";
    }
}

// -------------------------------------------------------
// cookMySop
// -------------------------------------------------------
OP_ERROR SOP_ApicFluidSolver::cookMySop(OP_Context& context) {
    flags().setTimeDep(true);
    if (lockInputs(context) >= UT_ERROR_ABORT)
        return error();

    fpreal t = context.getTime();
    apic::SimParams params = readParams(t);

    if (!solver_ || t < 1e-5f) {
        solver_.reset();
        collision_.reset();
        collision_ = std::make_unique<apic::CollisionHandler>();
        solver_ = std::make_unique<apic::ApicSolver>(params);
        solver_->setCollisionHandler(collision_.get());
        buildCollisionSDF(nullptr);
        const GU_Detail* inputGeo0 = inputGeo(0);
        if (inputGeo0 && inputGeo0->getNumPoints() > 0)
            loadParticles(inputGeo0);
        lastCookedTime_ = t;
    }
    else {
        solver_->setParams(params);
    }

    solver_->stepFrame();
    lastCookedTime_ = t;

    gdp->clearAndDestroy();
    writeParticles(gdp);

    unlockInputs();
    return error();
}

// -------------------------------------------------------
// Parameter reading
// -------------------------------------------------------
apic::SimParams SOP_ApicFluidSolver::readParams(fpreal t) const {
    apic::SimParams p;
    p.particleSeparation = static_cast<float>(evalFloat(PARM_PARTICLE_SEP, 0, t));
    p.timeScale          = static_cast<float>(evalFloat(PARM_TIME_SCALE,   0, t));
    p.substeps           = evalInt(PARM_SUBSTEPS, 0, t);
    p.gridSpacing = p.particleSeparation * 2.0f;
    p.gridResX = p.gridResY = p.gridResZ = evalInt(PARM_GRID_RES, 0, t);
    float halfSize = static_cast<float>(p.gridResX) * p.gridSpacing * 0.5f;
    p.gridOrigin = apic::vec3(-halfSize, -halfSize, -halfSize);

    switch (evalInt(PARM_TRANSFER_METHOD, 0, t)) {
        case 0:  p.method = apic::TransferMethod::APIC;   break;
        case 1:  p.method = apic::TransferMethod::FLIP;   break;
        case 2:  p.method = apic::TransferMethod::PIC;    break;
        default: p.method = apic::TransferMethod::Hybrid; break;
    }

    p.viscosity      = static_cast<float>(evalFloat(PARM_VISCOSITY,       0, t));
    p.surfaceTension = static_cast<float>(evalFloat(PARM_SURFACE_TENSION, 0, t));

    p.gravity = apic::vec3(
        static_cast<float>(evalFloat(PARM_GRAVITY, 0, t)),
        static_cast<float>(evalFloat(PARM_GRAVITY, 1, t)),
        static_cast<float>(evalFloat(PARM_GRAVITY, 2, t))
    );

    p.dt = 1.0f / 24.0f;
    return p;
}

// -------------------------------------------------------
// Geometry I/O
// -------------------------------------------------------
bool SOP_ApicFluidSolver::loadParticles(const GU_Detail* src) {
    solver_->particles().clear();

    GA_ROHandleV3 velH (src->findAttribute(GA_ATTRIB_POINT, "v"));
    GA_ROHandleF  massH(src->findAttribute(GA_ATTRIB_POINT, "mass"));
    const float   defaultMass = 1.0f;

    // Iterate over all points
    for (GA_Offset ptoff : src->getPointRange()) {
        UT_Vector3F pos = src->getPos3(ptoff);
        UT_Vector3F vel(0.f, 0.f, 0.f);
        float       m = defaultMass;

        if (velH.isValid())  vel = velH.get(ptoff);
        if (massH.isValid()) m   = massH.get(ptoff);

        solver_->particles().addParticle(
            apic::vec3(pos.x(), pos.y(), pos.z()),
            apic::vec3(vel.x(), vel.y(), vel.z()),
            m
        );
    }
    return true;
}

void SOP_ApicFluidSolver::writeParticles(GU_Detail* dst) {
    const auto& P  = solver_->particles();
    const size_t np = P.size();
    if (np == 0) return;

    GA_Offset ptStart = dst->appendPointBlock(static_cast<GA_Size>(np));

    GA_RWHandleV3 velH (dst->addFloatTuple(GA_ATTRIB_POINT, "v",       3));
    GA_RWHandleF  massH(dst->addFloatTuple(GA_ATTRIB_POINT, "mass",    1));
    GA_RWHandleFA bmH  (dst->addFloatTuple(GA_ATTRIB_POINT, "bmatrix", 9));

    for (size_t p = 0; p < np; ++p) {
        GA_Offset ptoff = ptStart + static_cast<GA_Offset>(p);

        const apic::Vec3& xp = P.position(p);
        dst->setPos3(ptoff, UT_Vector3F(xp.x(), xp.y(), xp.z()));

        if (velH.isValid()) {
            const apic::Vec3& vp = P.velocity(p);
            velH.set(ptoff, UT_Vector3F(vp.x(), vp.y(), vp.z()));
        }
        if (massH.isValid())
            massH.set(ptoff, P.mass(p));

        if (bmH.isValid()) {
            const apic::Mat3& B = P.affineB(p);
            UT_FloatArray bvals(9, 9);
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    bvals(r*3+c) = B(r, c);
            bmH.set(ptoff, bvals);
        }
    }

    dst->bumpDataIdsForAddOrRemove(true, true, true);
}

void SOP_ApicFluidSolver::buildCollisionSDF(const GU_Detail* /*collGeo*/) {
    // Domain boundary is handled by ApicGrid::markBoundary() (solid_ flags).
    // No SDF obstacle needed for the basic box domain.
    // Future: add internal obstacles here from collGeo VDB.
}

// -------------------------------------------------------
// DSO entry point
// -------------------------------------------------------
void newSopOperator(OP_OperatorTable* table) {
    table->addOperator(new OP_Operator(
        SOP_ApicFluidSolver::theSOPTypeName,
        "APIC Fluid Solver",
        SOP_ApicFluidSolver::myConstructor,
        SOP_ApicFluidSolver::myTemplateList,
        1,       // min inputs
        2,       // max inputs
        nullptr  // variables
    ));
}
