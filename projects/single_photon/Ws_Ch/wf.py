#%% standard cdft
from projects.antisiteQubit.wf_defect import ZPLWF
from fireworks import LaunchPad
import os
from atomate.vasp.powerups import *

CATEGORY = "cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))
"""
NBVN:
    ENCUT = 336
    FERWE = 329*1 1*0 1*1 144*0
    FERDO = 328*1 147*0
    NBANDS = 475
    A="/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869/
    launcher_2021-01-13-02-58-44-383950"
"""

anti_triplet = ZPLWF("/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869/"
                     "launcher_2021-01-13-02-58-44-383950", None)

wf = anti_triplet.wf(
    "B-C-D", 0, up_occupation="329*1 1*0 1*1 144*0",
    down_occupation="328*1 147*0", nbands=475, gamma_only=True, selective_dyn=None, specific_structure=None
)

wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

wf = add_modify_incar(wf)

wf = set_execution_options(wf, category=CATEGORY)
wf = preserve_fworker(wf)
LPAD.add_wf(wf)
#%% make interpolate

from qubitPack.tool_box import get_db
from qubitPack.tool_box import get_interpolate_sts
from pymatgen import Structure
import numpy as np

cdft = get_db("single_photon_emitter", "cdft")

b_st = Structure.from_dict(cdft.collection.find_one({"task_id":158})["input"]["structure"])
c_st = Structure.from_dict(cdft.collection.find_one({"task_id":160})["output"]["structure"])

dx = np.linspace(0, 2,11)
get_interpolate_sts(b_st,
                    c_st,
                    disp_range=dx,
                    output_dir='/Users/jeng-yuantsai/Research/project/single_photon/calculations/cdft/NBVN/'
                               'interpolate_st/excited_state/')

#%% run calcs for parabola

import glob
from fireworks import LaunchPad
import os
from projects.antisiteQubit.wf_defect import ZPLWF
from monty.serialization import loadfn
from pymatgen.io.vasp.inputs import Poscar
from pymatgen import Structure
from atomate.vasp.powerups import *

A = "/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869/launcher_2021-01-20-22-05-19-931342"

CATEGORY = "cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))
"""

"""
def antistie_triplet_ZPL(sts_path, state="excited_state"): #6.5
    #efrc

    anti_triplet = ZPLWF(A, None)

    # sd_sites = anti_triplet.set_selective_sites(25, 5)
    # sd_sites.sort()
    # print(sd_sites)
    process = None
    if state == "ground_state":
        process = "D"
    elif state == "excited_state":
        process = "B"
    else:
        return

    info = loadfn(os.path.join(sts_path, "{}/info.json".format(state)))

    for poscar in glob.glob(os.path.join(sts_path, "{}/*.vasp".format(state))):

        st = Structure.from_file(poscar)

        wf = anti_triplet.wf(
            process, 0, up_occupation="98*1 1*0 1*1 50*0",
            down_occupation="98*1 52*0", nbands=150, gamma_only=True, selective_dyn=None,
            specific_structure=Poscar.from_file(poscar)
        )

        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

        wf = add_modify_incar(wf)
        wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")
        wf = add_additional_fields_to_taskdocs(
            wf,
            {
                "cdft_info": "{}_NBVN".format(state),
                "interpolate_st_idx":poscar.split("/")[-1], "sts_info": info,
                "poscar_idx": float(poscar.split("/")[-1].split("_")[-1].split(".")[0])/10
            })
        wf.name = wf.name+":{}".format(poscar.split("/")[-1])
        wf = set_execution_options(wf, category=CATEGORY)
        wf = preserve_fworker(wf)
        LPAD.add_wf(wf)

antistie_triplet_ZPL(sts_path="/home/tug03990/work/single_photon_emitter/cdft/NBVN/interpolate_st", state="ground_state")

#%% read energy
from qubitPack.tool_box import get_db
import pandas as pd

state = "ground_state"

es = get_db("single_photon_emitter", "cdft", port=1234).collection.aggregate([
    {"$match": {"cdft_info": "{}_NBVN".format(state), "poscar_idx": {"$exists":True}}},
    {"$project": {"energy": "$output.energy", "poscar_idx":1, "task_id":1, "_id":0}}
])

df = pd.DataFrame(es)
df = df.sort_values("poscar_idx")
df.to_clipboard()
