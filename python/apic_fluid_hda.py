"""
apic_fluid_hda.py
Python scripts embedded in the APIC Fluid Solver HDA.

These callbacks are attached to the HDA's Python module tab
(Edit > Type Properties > Scripts > Python Module).

Usage inside Houdini:
  - OnCreated : hou.phm().on_created(kwargs)
  - OnInputChanged : hou.phm().on_input_changed(kwargs)
  - Button callbacks are wired via the parameter's callback script.
"""

import hou
import os


# -------------------------------------------------------
# HDA lifecycle callbacks
# -------------------------------------------------------

def on_created(kwargs):
    """Called once when the HDA node is first placed in the network."""
    node = kwargs["node"]
    _set_default_gravity(node)
    _init_progress_bar(node)
    hou.ui.displayMessage(
        "APIC Fluid Solver ready.\n"
        "Connect fluid source geometry to input 0, "
        "then press Simulate.",
        title="APIC Fluid Solver",
        severity=hou.severityType.Message
    )


def on_input_changed(kwargs):
    """Called when an input connection changes."""
    node = kwargs["node"]
    input_index = kwargs.get("input_index", -1)
    if input_index == 1:
        # Collision geometry connected/disconnected
        has_collision = node.input(1) is not None
        node.parm("collision_info").set(
            "Collision geometry connected." if has_collision
            else "No collision geometry."
        )


# -------------------------------------------------------
# Button callbacks (wired in parameter panel)
# -------------------------------------------------------

def on_simulate(kwargs):
    """Run simulation for the current frame range."""
    node = kwargs["node"]
    with _progress_context(node, "Simulating..."):
        _cook_node(node)


def on_reset(kwargs):
    """Reset the simulation state."""
    node = kwargs["node"]
    # Force re-initialise by setting a special flag attribute
    # that the C++ SOP reads to clear its internal state.
    node.parm("_reset_flag").set(
        node.parm("_reset_flag").eval() + 1
    )
    node.cook(force=True)
    hou.ui.displayMessage("Simulation reset.", title="APIC Fluid")


def on_cache_to_disk(kwargs):
    """Bake the simulation to disk as a BGEO sequence."""
    node = kwargs["node"]

    cache_path = hou.ui.selectFile(
        title="Choose cache directory",
        file_type=hou.fileType.Directory,
        chooser_mode=hou.chooserMode.Write
    )
    if not cache_path:
        return

    cache_path = cache_path.rstrip("/\\")
    os.makedirs(hou.expandString(cache_path), exist_ok=True)

    frame_range = hou.playbar.frameRange()
    start, end = int(frame_range[0]), int(frame_range[1])

    with hou.InterruptableOperation(
        "Caching APIC fluid",
        long_operation_name="Writing BGEO frames",
        open_interrupt_dialog=True
    ) as op:
        for frame in range(start, end + 1):
            if op.isInterrupted():
                break

            hou.setFrame(frame)
            node.cook(force=True)

            out_path = "{}/apic_fluid.{}.bgeo.sc".format(cache_path, str(frame).zfill(4))
            gdp = node.geometry()
            gdp.saveToFile(out_path)

            progress = (frame - start) / max(end - start, 1)
            op.updateLongProgress(
                progress,
                "Frame {} / {}".format(frame, end)
            )

    hou.ui.displayMessage(
        "Cache written to:\n{}".format(cache_path),
        title="APIC Fluid — Cache Complete"
    )


# -------------------------------------------------------
# Visualization helpers
# -------------------------------------------------------

def update_visualization(kwargs):
    """Called when visualization toggles change."""
    node = kwargs["node"]
    show_particles = node.parm("showparticles").eval()
    show_grid      = node.parm("showgrid").eval()
    show_vectors   = node.parm("showvectors").eval()
    color_by_vel   = node.parm("colorbyvel").eval()

    # Update display flags on downstream visualization nodes if they exist
    # (This is a stub — in a full HDA, downstream viz SOPs would be
    #  wired to these parameters directly via channel references.)
    _update_color_viz(node, color_by_vel)


def _update_color_viz(node, color_by_vel):
    """Toggle Cd attribute generation on/off."""
    # In production: set a parameter on an inner AttributeWrangle node
    # that sets Cd = normalize(v) * speed or Cd = {1,1,1}
    pass


# -------------------------------------------------------
# Transfer method comparison helper
# -------------------------------------------------------

def compare_methods(kwargs):
    """
    Utility: duplicate the node three times (APIC / FLIP / PIC)
    and lay them out side by side for comparison.
    """
    node = kwargs["node"]
    parent = node.parent()

    methods = [("APIC", 0), ("FLIP", 1), ("PIC", 2)]
    offset_x = 3.0

    created = []
    for name, method_idx in methods:
        dup = parent.createNode("apic_fluid_solver", "{}_{}".format(
            node.name(), name
        ))
        dup.setInput(0, node.input(0))
        if node.input(1):
            dup.setInput(1, node.input(1))

        # Copy all parameters then override method
        for parm in node.parms():
            try:
                dup.parm(parm.name()).set(parm.eval())
            except Exception:
                pass
        dup.parm("transfermethod").set(method_idx)

        # Layout
        pos = node.position()
        dup.setPosition(hou.Vector2(pos[0] + offset_x * (method_idx - 1), pos[1] - 2))
        created.append(dup)

    hou.ui.displayMessage(
        "Created {} comparison nodes.".format(len(created)),
        title="APIC Fluid"
    )


# -------------------------------------------------------
# Internal helpers
# -------------------------------------------------------

def _set_default_gravity(node):
    node.parmTuple("gravity").set((0, -9.8, 0))


def _init_progress_bar(node):
    pass  # progress is shown via the parameter pane label


def _cook_node(node):
    node.cook(force=True)


class _progress_context:
    def __init__(self, node, msg):
        self._node = node
        self._msg  = msg

    def __enter__(self):
        self._node.setComment(self._msg)
        return self

    def __exit__(self, *_):
        self._node.setComment("")
