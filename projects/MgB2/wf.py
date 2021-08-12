#%% get dos wf
import glob

from atomate.vasp.fireworks import NonSCFFW
from atomate.vasp.powerups import *
import os
from fireworks import LaunchPad, Workflow

tk143 = "/home/tug03990/work/mgb2/hetero/block_2020-12-17-05-32-06-074371/launcher_2020-12-18-05-28-48-950427"

CATEGORY = "hetero"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/mgb2/{}/my_launchpad.yaml".format(CATEGORY)))

bs_fw = NonSCFFW(prev_calc_dir=tk143, mode="uniform", structure=Structure.from_file(os.path.join(tk143,"CONTCAR.gz")))

wf = Workflow([bs_fw])

wf = add_modify_incar(wf, {"incar_update": {"NEDOS":9000, "ENMAX":10, "ENMIN":-10}})

wf = add_additional_fields_to_taskdocs(wf, {"task_label":"bs_uniform", "terminate":"c_terminate"})
wf = set_execution_options(wf, category=CATEGORY)

LPAD.add_wf(wf)

#%%
#%% WF
import numpy as np
from pymatgen.io.vasp.sets import MPRelaxSet, MPMetalRelaxSet
from atomate.vasp.fireworks import OptimizeFW, StaticFW
from atomate.vasp.powerups import *
from fireworks import Workflow, LaunchPad, FileTransferTask
from atomate.vasp.fireworks import NonSCFFW
from atomate.vasp.powerups import *
import os
from fireworks import LaunchPad, Workflow
import os, glob

CATEGORY = "hetero"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/mgb2/{}/my_launchpad.yaml".format(CATEGORY)))


# c_term = Structure.from_file(os.path.join("models/heter", "two_MgO-C_term_SiC-O_bottom.vasp"))
# si_term = Structure.from_file(os.path.join("models/heter", "two_MgO-Si_term_SiC-O_bottom.vasp"))
# sic = Structure.from_file(os.path.join("models/unit_cell", "SiC_terminated.vasp"))
# mgo = Structure.from_file(os.path.join("models/unit_cell", "mgo.vasp"))

# sic = Structure.from_file(os.path.join("models/heter", "subs.vasp"))
mgo = Structure.from_file(os.path.join("models/heter/c_term", "mat2d.vasp"))


def dp(structure):
    weights = [s.species.weight for s in structure]
    center_of_mass = np.average(
        structure.frac_coords, weights=weights, axis=0
    )

    return {"IDIPOL":3, "LDIPOL": True, "DIPOL": center_of_mass}


config = {"subs": {-1: "si_term", 1: "c_term"}, "mat2d": {-1: "o_first", 1: "mg_first"}}

for subs in config["subs"]:
    for mat2d in config["mat2d"]:
        path = "models/heter/{}/{}".format(config["subs"][subs], config["mat2d"][mat2d])
        term = config["subs"][subs]
        first_layer = config["mat2d"][mat2d]
        for st in glob.glob(os.path.join(path, "*.vasp")):
            print(st)
            st = Structure.from_file(st)
            fws = []

            opt = OptimizeFW(st, vasp_input_set=MPRelaxSet(st, vdw="optPBE"))

            static = StaticFW(st, parents=opt,
                              vasptodb_kwargs={
                                  "parse_eigenvalues": False,
                                  "parse_dos": False,
                              })

            # bs = NonSCFFW(parents=static, mode="uniform", input_set_overrides=dict(reciprocal_density=250, nedos=9000), structure=st)

            fws.extend([static, opt])
            wf = Workflow(fws, name=term+":{}".format(st.formula))

            wf = add_additional_fields_to_taskdocs(wf, {
                "terminate": term, "first_layer": first_layer, "metal":False, "project": "MgO", "spin":1, "vdw":True
            })

            wf = add_tags(wf, [{
                "terminate": term, "first_layer": first_layer, "metal":False, "project": "MgO", "spin":1, "vdw":True
            }])

            uis = {"EDIFF":1E-4, "EDIFFG":-0.01, "ISIF":2}
            uis.update(dp(st))

            wf = add_modify_incar(wf, {"incar_update": uis}, fw_name_constraint=fws[0].name)
            uis.update({"EDIFF":1E-5})
            wf = add_modify_incar(wf, {"incar_update": uis}, fw_name_constraint=fws[1].name)
            # wf = add_modify_incar(wf, {"incar_update": {"ENMAX":-15, "ENMAX":15}},fw_name_constraint=fws[-1].name)
            wf = add_modify_incar(wf, {"incar_update": dict(ISPIN=1, MAGMOM=MPRelaxSet(st).incar.get("MAGMOM"))})
            wf = preserve_fworker(wf)
            wf = set_execution_options(wf, category=CATEGORY)

            print(wf)
            # LPAD.add_wf(wf)


