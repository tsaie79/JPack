#%% lattice relax
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
from qubitPack.tool_box import get_db

CATEGORY = "modeling"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/copper_graphene/{}/my_launchpad.yaml".format(CATEGORY)))

fws = []
for st in ["3-cu.vasp"]:#, "copper-gr.vasp"]:
    st = Structure.from_file("/home/tug03990/scratch/copper_graphene/modeling/{}".format(st))
    # db = get_db("copper_graphene", "modeling", port=12345)
    # st = Structure.from_dict(db.collection.find_one({"task_id":25})["output"]["structure"])
    # st.make_supercell([3,3,1])
    # print(st.num_sites)

    override_default_vasp_params={
        "user_incar_settings": {"EDIFF":1E-6, "EDIFFG":-0.001, "ISIF":3, "SIGMA":0.001, "ISPIN":1},
        "vdw": "df"}
    # override_default_vasp_params={
    #     "user_incar_settings": {"EDIFF":1E-5, "EDIFFG":-0.01, "ISIF":2, "ISPIN":1, "SIGMA":0.1},
    #     # "user_kpoints_settings": {"reciprocal_density":64},
    #     "vdw": "df"}

    vis_set = MPMetalRelaxSet(st, force_gamma=True, **override_default_vasp_params)
    # vis_set = MPRelaxSet(st, force_gamma=True, **override_default_vasp_params)

    opt_vdw = OptimizeFW(st, vasp_input_set=vis_set, vasp_cmd=">>vasp_cmd<<")

    fws.append(opt_vdw)

wf = Workflow(fws, name="{}_PBE_vdw_DF_relax".format(st.formula))
wf = add_additional_fields_to_taskdocs(wf, {"main":True})
wf = preserve_fworker(wf)
wf = set_execution_options(wf, category=CATEGORY)
wf = add_modify_incar(wf)


LPAD.add_wf(wf)

#%% make wf for supercell relaxation and then MD

from atomate.vasp.fireworks import MDFW, OptimizeFW
from atomate.vasp.powerups import *
from pymatgen.io.vasp.sets import MPMDSet, MPMetalRelaxSet, MPRelaxSet
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Kpoints
from qubitPack.tool_box import get_db
from fireworks import LaunchPad, Workflow
import os

def rnd_sampling(st):
    from random import sample
    rnd_site = sample(range(384), 126)
    for site in rnd_site:
        st.replace(site, "C")
    return st

CATEGORY = "modeling"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/copper_graphene/{}/my_launchpad.yaml".format(CATEGORY)))


db = get_db("copper_graphene", "modeling", port=12345)
st = Structure.from_dict(db.collection.find_one({"task_id":37})["output"]["structure"])


st.make_supercell([4, 4, 1])
print(st.num_sites)

override_default_vasp_params={
    "user_incar_settings": {"EDIFF":1E-5, "EDIFFG":-0.01, "ISIF":2, "SIGMA":0.04, "ISPIN":1},
    "vdw": "df"}

vis_set = MPMetalRelaxSet(st, force_gamma=True, **override_default_vasp_params)
opt_vdw = OptimizeFW(st, vasp_input_set=vis_set, vasp_cmd=GAMMA_VASP_CMD)
##MD
Q = 2
md_algo = 2
isif = 2
nsteps = 10000 #10ps

fws = []
fws.append(opt_vdw)
override_default_vasp_params={
    "user_incar_settings": {"MDALGO":md_algo, "ISIF":isif, "SMASS":Q, "SIGMA":0.04, "ISMEAR":1, "ISPIN":1, "POTIM":1}}
for T in [500, 800, 1200, 1400, 1600, 1800, 2000, 2200]:
    vis = MPMDSet(st, T, T, nsteps, True, **override_default_vasp_params)
    fw = MDFW(st, T, T, nsteps, vasp_input_set=vis, copy_vasp_outputs=False, vasp_cmd=GAMMA_VASP_CMD, parents=fws[0])
    fw.name += ":{}".format(T)
    fws.append(fw)

wf = Workflow(fws, name="{}".format(st.formula))
wf = add_modify_kpoints(wf, {"kpoints_update": {"kpts":[[1,1,1]]}}, fw_name_constraint="molecular dynamics")
wf = add_additional_fields_to_taskdocs(wf, {"main":True})
wf = remove_custodian(wf, fw_name_constraint=wf.fws[0].name)
wf = preserve_fworker(wf)
wf = set_execution_options(wf, category=CATEGORY)
wf = add_modify_incar(wf)

LPAD.add_wf(wf)

