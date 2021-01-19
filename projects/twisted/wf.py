#%%
from atomate.vasp.fireworks.core import ScanOptimizeFW, StaticFW, NonSCFFW, OptimizeFW
from atomate.vasp.workflows.presets.core import wf_bandstructure
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    add_modify_kpoints,
    set_queue_options,
    set_execution_options,
    clean_up_files,
    modify_gzip_vasp
)
from atomate.vasp.database import VaspCalcDb

from pymatgen import Structure
from pymatgen.io.vasp.sets import MPRelaxSet

from fireworks import LaunchPad, Workflow

from monty.serialization import loadfn

import os

def ML_bs_wf(cat="Bi2Se3"):

    def bs_fws(structure):
        # opt = OptimizeFW(structure=structure)
        static_fw = StaticFW(structure=structure)
        line_fw = NonSCFFW(structure=structure,
                           mode="line",
                           parents=static_fw,
                           input_set_overrides={"other_params": {"two_d_kpoints": True}}
                           )

        wf = Workflow([static_fw, line_fw], name="{}:pbe_bs".format(structure.formula))

        updates = {
            # "ADDGRID": True,
            # "LASPH": True,
            # "LDAU": False,
            # "LMIXTAU": True,
            # "METAGGA": "SCAN",
            "NELM": 200,
            "EDIFF": 1E-5,
            "ISPIN": 1,
            "LAECHG": False,
        }

        wf = add_modify_incar(wf, {"incar_update": updates})
        # wf = add_modify_incar(wf, {"incar_update": {"LCHARG":False, "ISIF":2, "EDIFFG":-0.01, "EDIFF":1E-4}}, opt.name)
        wf = add_modify_incar(wf, {"incar_update": {"LCHARG":True, "LVHAR":True, "LWAVE":False, "NEDOS":9000,
                                                    "EMIN":-10, "EMAX":10}}, static_fw.name)
        wf = add_modify_incar(wf, {"incar_update": {"LWAVE":False, "LCHARG":False}}, line_fw.name)
        wf = clean_up_files(wf, files=["CHG*", "DOS*", "LOCPOT*"], fw_name_constraint=line_fw.name,
                            task_name_constraint="VaspToDb")
        return wf


    lpad = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/twisted/{}/my_launchpad.yaml".format(cat))))

    st = Structure.from_file("/home/tug03990/work/twisted/Bi2Se3/Bi2Se3_5_10.vasp")
    wf = bs_fws(st)

    encut = 1.3*max([potcar.enmax for potcar in MPRelaxSet(st).potcar])
    wf = add_modify_incar(wf, {"incar_update": {"ENCUT": encut}})

    wf = add_modify_incar(wf)
    wf = set_execution_options(wf, category=cat)
    wf = preserve_fworker(wf)
    wf = add_additional_fields_to_taskdocs(wf, {"twisted_angle":5.1})

    lpad.add_wf(wf)
    print(wf.name)

ML_bs_wf()