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

# st = Structure.from_file("/home/tug03990/scratch/copper_graphene/modeling/cu-graphene-cu.vasp")
db = get_db("copper_graphene", "modeling", port=12345)
st = Structure.from_dict(db.collection.find_one()["output"]["structure"])
st.make_supercell([3,3,1])
print(st.num_sites)

# override_default_vasp_params={
#     "user_incar_settings": {"EDIFF":1E-6, "EDIFFG":-0.005, "ISIF":3},
#     "vdw": "df"}
override_default_vasp_params={
    "user_incar_settings": {"EDIFF":1E-4, "EDIFFG":-0.01, "ISIF":2, "ISPIN":1},
    "user_kpoints_settings": {"reciprocal_density":64},
    "vdw": "df"}

vis_set = MPMetalRelaxSet(st, force_gamma=True, **override_default_vasp_params)
# vis_set = MPRelaxSet(st, force_gamma=True, **override_default_vasp_params)

opt_vdw = OptimizeFW(st, vasp_input_set=vis_set, vasp_cmd=GAMMA_VASP_CMD)

fws = [opt_vdw]
wf = Workflow(fws, name="{}_PBE_vdw_DF_relax".format(st.formula))

wf = preserve_fworker(wf)
wf = set_execution_options(wf, category=CATEGORY)
wf = add_modify_incar(wf)


LPAD.add_wf(wf)

#%% MD
from atomate.vasp.fireworks import MDFW
from atomate.vasp.powerups import *
from pymatgen.io.vasp.sets import MPMDSet, MPMetalRelaxSet, MPRelaxSet
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Kpoints
from qubitPack.tool_box import get_db
from fireworks import LaunchPad, Workflow
import os

CATEGORY = "modeling"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/copper_graphene/{}/my_launchpad.yaml".format(CATEGORY)))

db = get_db("copper_graphene", "modeling", port=12345)
st = Structure.from_dict(db.collection.find_one({"task_id":10})["output"]["structure"])

Q = 2
md_algo = 2
isif = 2
nsteps = 5000 #10ps

override_default_vasp_params={
    "user_incar_settings": {"MDALGO":md_algo, "ISIF":isif, "SMASS":Q, "SIGMA":0.04, "ISMEAR":1, "ISPIN":1, "POTIM":2}}
for T in [800, 2200]:
    vis = MPMDSet(st, T, T, nsteps, True, **override_default_vasp_params)

    fw = MDFW(st, T, T, nsteps, vasp_input_set=vis, copy_vasp_outputs=False, vasp_cmd=GAMMA_VASP_CMD)

    fws = [fw]
    wf = Workflow(fws, name="{}:{}K".format(st.formula, T))

    wf = add_modify_kpoints(wf, {"kpoints_update": {"kpts":[[1,1,1]]}})

    wf = preserve_fworker(wf)
    wf = set_execution_options(wf, category=CATEGORY)
    wf = add_modify_incar(wf)

    LPAD.add_wf(wf)

#%% make wf for supercell relaxation and then MD
from pymatgen.io.vasp.sets import MPRelaxSet, MPMetalRelaxSet
from atomate.vasp.fireworks import OptimizeFW, StaticFW
from atomate.vasp.powerups import *
from fireworks import Workflow, LaunchPad, FileTransferTask
from atomate.vasp.fireworks import NonSCFFW
from atomate.vasp.powerups import *
import os
from fireworks import LaunchPad, Workflow
import os

CATEGORY = "modeling"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/copper_graphene/{}/my_launchpad.yaml".format(CATEGORY)))

db = get_db("copper_graphene", "modeling", port=12345)
st = Structure.from_dict(db.collection.find_one({"task_id":10})["output"]["structure"])
# st.make_supercell([4,4,1])

override_default_vasp_params={
    "user_incar_settings": {"EDIFF":1E-4, "EDIFFG":-0.01, "ISIF":2},
    "vdw": "df"}

vis_set = MPMetalRelaxSet(st, force_gamma=True, **override_default_vasp_params)

opt_vdw = OptimizeFW(st, vasp_input_set=vis_set)

fws = [opt_vdw]

Q = 2
md_algo=2
isif = 2
T = 500
nsteps = 10000 #10ps

override_default_vasp_params={
    "user_incar_settings": {"MDALGO":md_algo, "ISIF":isif, "SMASS":Q, "POTIM":1, "SIGMA":0.04}}
vis = MPMDSet(st, T, T, nsteps, True, **override_default_vasp_params)

md = MDFW(st, T, T, nsteps, vasp_input_set=vis, copy_vasp_outputs=False, parents=fws[0])

fws.append(md)
wf = Workflow(fws, name="{}:{}K".format(st.formula, T))

wf = preserve_fworker(wf)
wf = set_execution_options(wf, category=CATEGORY)
wf = add_modify_incar(wf)

# LPAD.add_wf(wf)