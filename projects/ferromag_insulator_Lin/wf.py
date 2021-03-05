pbe relax
hse static

#%%
import numpy as np
from pymatgen.io.vasp.sets import MPRelaxSet, MPMetalRelaxSet
from pymatgen.io.vasp.inputs import Kpoints
from atomate.vasp.fireworks.jcustom import JHSEStaticFW
from atomate.vasp.fireworks import OptimizeFW
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import cp_vdw_file
from fireworks import LaunchPad, Workflow
import os

CATEGORY = "slab"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/mag_insulator/{}/my_launchpad.yaml".format(CATEGORY)))

for f in ["3-cu.vasp"]:#, "copper-gr.vasp"]:
    fws = []

    st = Structure.from_file("/home/tug03990/scratch/copper_graphene/modeling/{}".format(f))
    encut = 1.3*max([potcar.enmax for potcar in MPRelaxSet(st).potcar])

    override_default_vasp_params={
        "user_incar_settings": {"EDIFF":1E-4, "EDIFFG":-0.01, "ISIF":3, "SIGMA":0.05, "ISPIN":1, "ENCUT":encut},
        "vdw": "df"}

    vis_set = MPMetalRelaxSet(st, force_gamma=True, **override_default_vasp_params)
    opt_vdw = OptimizeFW(st, vasp_input_set=vis_set, vasp_cmd=">>vasp_cmd<<")
    opt_vdw = cp_vdw_file(opt_vdw)


    hse_static = JHSEStaticFW(
        st, vasp_input_set_params={
            "user_incar_settings": {"EDIFF":1E-5, "SIGMA":0.05, "ENCUT":encut},
            "user_kpoints_settings": Kpoints.automatic_density_by_vol(st, 6**3)
        },
        parents=opt_vdw
    )

    fws.append(opt_vdw)
    fws.append(hse_static)

    wf = Workflow(fws, name="{}_{}".format(st.formula, f))
    wf = add_additional_fields_to_taskdocs(wf, {"source": f})
    wf = preserve_fworker(wf)
    wf = set_execution_options(wf, category=CATEGORY)
    wf = add_modify_incar(wf)

    LPAD.add_wf(wf)