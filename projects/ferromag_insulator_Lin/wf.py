#%%
import numpy as np
from pymatgen.io.vasp.sets import MPRelaxSet, MPMetalRelaxSet
from pymatgen.io.vasp.inputs import Kpoints
from atomate.vasp.fireworks.jcustom import JHSEStaticFW
from atomate.vasp.fireworks import OptimizeFW
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import cp_vdw_file
from fireworks import LaunchPad, Workflow
import os, glob

CATEGORY = "slab"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/mag_insulator/{}/my_launchpad.yaml".format(CATEGORY))))

f_path = "/home/tug03990/config/project/mag_insulator/slab/structures/tci_slab_hex"

for f in glob.glob(os.path.join(f_path, "icsd-157074_mp-169.vasp")):
    print(f.split("/")[-1])
    fws = []

    st = Structure.from_file(f)
    encut = 1.3*max([potcar.enmax for potcar in MPRelaxSet(st).potcar])

    override_default_vasp_params = {
        "user_incar_settings": {"EDIFF":1E-4, "EDIFFG":-0.01, "ISIF":3, "SIGMA":0.05, "ISPIN":1,
                                "ENCUT":encut, "LASPH":False, "ISMEAR":0,
                                "LCHARG":False},
        # "user_kpoints_settings": Kpoints.automatic_density_by_vol(st, 6**3),
        "vdw": None}

    vis_set = MPMetalRelaxSet(st, force_gamma=True, **override_default_vasp_params)
    opt_vdw = OptimizeFW(st, vasp_input_set=vis_set, vasp_cmd=">>vasp_cmd<<")

    hse_static = JHSEStaticFW(
        st, vasp_input_set_params={
            "user_incar_settings": {"EDIFF":1E-5, "ISMEAR":0, "SIGMA":0.05, "ENCUT":encut,
                                    "LASPH":False, "ISPIN":1,
                                    "LVHAR":True, "LCHARG": True},
            "user_kpoints_settings": Kpoints.automatic_density_by_vol(st, 6**3)
        },
        vasptodb_kwargs={"parse_chgcar":True},
        parents=opt_vdw
    )

    fws.append(opt_vdw)
    fws.append(hse_static)

    wf = Workflow(fws, name="{}_{}".format(st.formula, f.split("/")[-1]))
    wf = add_additional_fields_to_taskdocs(wf, {"source": f.split("/")[-1]})
    opt_vdw = cp_vdw_file(wf, fw_name_constraint=opt_vdw.name)

    wf = preserve_fworker(wf)
    wf = set_execution_options(wf, category=CATEGORY)
    wf = add_modify_incar(wf)

    wf = set_queue_options(wf, "32:00:00", fw_name_constraint=wf.fws[0].name)
    wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)

    LPAD.add_wf(wf)
