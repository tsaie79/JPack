from qubitPack.workflow.tool_box import *

from fireworks import LaunchPad, Workflow

from atomate.vasp.database import VaspCalcDb
from atomate.vasp.workflows.jcustom.wf_full import get_wf_full_hse, get_wf_full_scan
from atomate.vasp.workflows.presets.core import wf_bandstructure
from atomate.vasp.fireworks.core import ScanOptimizeFW
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    add_modify_kpoints,
    set_queue_options,
    set_execution_options,
    clean_up_files
)

from pymatgen import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPScanRelaxSet, MPRelaxSet
from pymatgen.symmetry.bandstructure import HighSymmKpath

from pycdt.core.defectsmaker import ChargedDefectsStructures

from monty.serialization import loadfn
import os


def test_IR(cat="genWavecar"):
    lpad = LaunchPad.from_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/testIR/{}/my_launchpad.yaml".format(cat)))
    col = VaspCalcDb.from_db_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/symBaseBinaryQubit/scan_relax_pc/db.json")).collection

    input_st = Structure.from_dict(col.find_one({"chemsys":"Mo-S"})["output"]["structure"])

    wf = wf_bandstructure(input_st)
    for fw_name in ["static", "nscf"]:
        wf = add_modify_incar(wf, {"incar_update":{"LWAVE":True, "ISYM":2}}, fw_name_constraint=fw_name)

    for ispin in [1,2]:
        wf = add_modify_incar(wf, {"incar_update":{"ISPIN":ispin}})
        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=cat)
        wf = preserve_fworker(wf)
        wf = set_queue_options(wf, "01:00:00")
        if ispin == 1:
            wf.name = "MoS2_spin_{}_2D_k".format("off")
        elif ispin == 2:
            wf.name = "MoS2_spin_{}_2D_k".format("on")
        lpad.add_wf(wf)

