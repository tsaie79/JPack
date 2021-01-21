#%% get dos wf
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
import os

CATEGORY = "hetero"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/mgb2/{}/my_launchpad.yaml".format(CATEGORY)))

terminate = "si_terminate"
st_p = "/home/tug03990/work/mgb2/hetero"

heter = Structure.from_file(os.path.join(st_p, terminate, "0.vasp"))
mgb2 = Structure.from_file(os.path.join(st_p, terminate, "mat2d.vasp"))
sic = Structure.from_file(os.path.join(st_p, terminate, "subs.vasp"))

def dp(structure):
    weights = [s.species.weight for s in structure]
    center_of_mass = np.average(
        structure.frac_coords, weights=weights, axis=0
    )

    return {"IDIPOL":3, "LDIPOL": True, "DIPOL": center_of_mass}


uis = {"ISIF": 2, "EDIFFG":-0.01}
uis.update(dp(heter))

opt_vdw = OptimizeFW(heter, override_default_vasp_params={"user_incar_settings": uis})
uis = {}
uis.update({"EDIFF":1E-5})
static = StaticFW(heter, parents=opt_vdw, vasp_input_set_params={"user_incar_settings": uis})

hetero_fws = [opt_vdw, static]
wf = Workflow(hetero_fws, name=terminate)

add_additional_fields_to_taskdocs(wf, {"terminate": terminate})
wf = add_modify_incar(wf, {"incar_update": dict(ISPIN=1, MAGMOM=MPRelaxSet(heter).incar.get("MAGMOM"))})
wf = preserve_fworker(wf)
wf = set_execution_options(wf, category=CATEGORY)

print(wf)
# LPAD.add_wf(wf)


for st in [heter]:
    fws = []

    opt = OptimizeFW(st, override_default_vasp_params={"user_incar_settings":uis}, vasp_input_set=MPMetalRelaxSet(st))

    static = StaticFW(st, parents=opt, vasp_input_set_params={"user_incar_settings": uis},
                      vasptodb_kwargs={
                          "parse_eigenvalues": False,
                          "parse_dos": False
                      })

    bs = NonSCFFW(parents=static, mode="uniform", input_set_overrides=dict(reciprocal_density=250, nedos=9000), structure=st)

    fws.extend([static, opt, bs])
    wf = Workflow(fws, name=terminate+":{}".format(st.formula))

    wf = add_additional_fields_to_taskdocs(wf, {"terminate": terminate, "metal":True})

    uis = {"EDIFF":1E-5, "EDIFFG":-0.01, "ISIF":2}
    uis.update(dp(st))

    wf = add_modify_incar(wf, {"incar_update": uis}, fw_name_constraint=fws[0].name)
    uis.update({"EDIFF":1E-6})
    wf = add_modify_incar(wf, {"incar_update": uis}, fw_name_constraint=fws[1].name)
    wf = add_modify_incar(wf, {"incar_update": {"ENMAX":-15, "ENMAX":15}},fw_name_constraint=fws[-1].name)
    wf = add_modify_incar(wf, {"incar_update": dict(ISPIN=1, MAGMOM=MPRelaxSet(st).incar.get("MAGMOM"))})
    wf = preserve_fworker(wf)
    wf = set_execution_options(wf, category=CATEGORY)

    print(wf)
    LPAD.add_wf(wf)


