"""
create_hda.py
Run this script inside a Houdini Python Shell to create the
APIC Fluid Solver HDA from scratch.

Usage:
  In Houdini → Windows → Python Shell:
    exec(open('/path/to/apic_fluid/python/create_hda.py').read())
"""

import hou
import os

HDA_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "hda", "apic_fluid_solver.hda"
)
os.makedirs(os.path.dirname(HDA_PATH), exist_ok=True)


def create_hda():
    # ---- Create a temporary SOP node to work from ----
    obj = hou.node("/obj")
    geo = obj.createNode("geo", "apic_build_temp")
    sop = geo.createNode("null", "apic_fluid_solver_base")

    # ---- Convert to HDA ----
    hda_node = sop.createDigitalAsset(
        name="apic_fluid_solver",
        hda_file_name=HDA_PATH,
        description="APIC Fluid Solver",
        min_num_inputs=1,
        max_num_inputs=2,
        version="0.1"
    )

    hda_def  = hda_node.type().definition()
    parm_grp = hda_def.parmTemplateGroup()

    # ---- Clear defaults ----
    for t in parm_grp.parmTemplates():
        parm_grp.remove(t.name())

    # ---- Simulation Parameters folder ----
    sim_folder = hou.FolderParmTemplate(
        "sim_params", "Simulation Parameters",
        folder_type=hou.folderType.Simple
    )
    sim_folder.addParmTemplate(hou.FloatParmTemplate(
        "particlesep", "Particle Separation", 1,
        default_value=(0.05,),
        min=0.005, max=0.5,
        help="Controls resolution. Smaller = more detail."
    ))
    sim_folder.addParmTemplate(hou.FloatParmTemplate(
        "timescale", "Time Scale", 1,
        default_value=(1.0,),
        min=0.01, max=5.0
    ))
    sim_folder.addParmTemplate(hou.IntParmTemplate(
        "substeps", "Substeps", 1,
        default_value=(2,),
        min=1, max=16,
        help="More substeps = more stable but slower."
    ))
    sim_folder.addParmTemplate(hou.IntParmTemplate(
        "gridres", "Grid Resolution", 1,
        default_value=(64,),
        min=8, max=256
    ))
    parm_grp.append(sim_folder)

    # ---- Physical Properties folder ----
    phys_folder = hou.FolderParmTemplate(
        "phys_props", "Physical Properties",
        folder_type=hou.folderType.Simple
    )
    phys_folder.addParmTemplate(hou.MenuParmTemplate(
        "transfermethod", "Transfer Method",
        menu_items=("apic", "flip", "pic", "hybrid"),
        menu_labels=("APIC", "FLIP", "PIC", "Hybrid"),
        default_value=0
    ))
    phys_folder.addParmTemplate(hou.FloatParmTemplate(
        "viscosity", "Viscosity", 1,
        default_value=(0.01,),
        min=0.0, max=1.0
    ))
    phys_folder.addParmTemplate(hou.FloatParmTemplate(
        "surfacetension", "Surface Tension", 1,
        default_value=(0.0,),
        min=0.0, max=1.0
    ))
    phys_folder.addParmTemplate(hou.FloatParmTemplate(
        "gravity", "Gravity", 3,
        default_value=(0.0, -9.8, 0.0),
        naming_scheme=hou.parmNamingScheme.XYZW
    ))
    parm_grp.append(phys_folder)

    # ---- Visualization folder ----
    viz_folder = hou.FolderParmTemplate(
        "visualization", "Visualization",
        folder_type=hou.folderType.Simple
    )
    viz_folder.addParmTemplate(hou.ToggleParmTemplate(
        "showparticles", "Show Particles", default_value=True
    ))
    viz_folder.addParmTemplate(hou.ToggleParmTemplate(
        "showgrid", "Show Grid", default_value=False
    ))
    viz_folder.addParmTemplate(hou.ToggleParmTemplate(
        "showvectors", "Show Vectors", default_value=False
    ))
    viz_folder.addParmTemplate(hou.ToggleParmTemplate(
        "colorbyvel", "Color by Velocity", default_value=False
    ))
    parm_grp.append(viz_folder)

    # ---- Action buttons ----
    buttons_folder = hou.FolderParmTemplate(
        "actions", "Actions",
        folder_type=hou.folderType.Simple
    )
    buttons_folder.addParmTemplate(hou.ButtonParmTemplate(
        "btn_reset", "Reset",
        script_callback="hou.phm().on_reset(kwargs)",
        script_callback_language=hou.scriptLanguage.Python
    ))
    buttons_folder.addParmTemplate(hou.ButtonParmTemplate(
        "btn_simulate", "Simulate",
        script_callback="hou.phm().on_simulate(kwargs)",
        script_callback_language=hou.scriptLanguage.Python
    ))
    buttons_folder.addParmTemplate(hou.ButtonParmTemplate(
        "btn_cache", "Cache to Disk",
        script_callback="hou.phm().on_cache_to_disk(kwargs)",
        script_callback_language=hou.scriptLanguage.Python
    ))
    parm_grp.append(buttons_folder)

    # ---- Hidden internal state ----
    parm_grp.append(hou.IntParmTemplate(
        "_reset_flag", "Reset Flag", 1,
        default_value=(0,),
        is_hidden=True
    ))

    hda_def.setParmTemplateGroup(parm_grp)

    # ---- Embed the Python module ----
    python_src_path = os.path.join(os.path.dirname(__file__), "apic_fluid_hda.py")
    with open(python_src_path, "r") as f:
        python_src = f.read()
    hda_def.addSection("PythonModule", python_src)

    # ---- Node icon / metadata ----
    hda_def.setIcon("SOP_fluid")
    hda_def.setComment(
        "APIC Fluid Solver HDA\n"
        "Implements the Affine Particle-In-Cell method\n"
        "(Jiang et al. 2015) for fluid simulation."
    )

    hda_def.save(HDA_PATH)
    print("HDA saved to: {}".format(HDA_PATH))

    # Clean up temp nodes
    geo.destroy()

    return hda_node


if __name__ == "__main__" or True:
    result = create_hda()
    print("HDA creation complete:", result)
