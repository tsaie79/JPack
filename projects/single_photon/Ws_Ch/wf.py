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


#%% soc standard defect
import glob
from fireworks import LaunchPad
import os
from projects.antisiteQubit.wf_defect import ZPLWF
from qubitPack.tool_box import get_db
from monty.serialization import loadfn
from pymatgen.io.vasp.inputs import Poscar
from pymatgen import Structure
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from pymatgen.io.vasp.sets import MPRelaxSet
from atomate.vasp.workflows.jcustom.wf_full import get_wf_full_hse

d_name = "single_photon_emitter"
col_name = "standard_defect"

CATEGORY = "soc_standard_defect"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))

std_defect = get_db(d_name, col_name, port=12345)
for tk_id in [264, 265, 266]:
    entry = std_defect.collection.find_one({"task_id": tk_id})

    scf_dir = entry["calcs_reversed"][0]["dir_name"].split("/")[-1]
    scf_dir = os.path.join("/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869", scf_dir)

    st = Structure.from_dict(entry["output"]["structure"])

    magmom = [[0,0,mag_z] for mag_z in MPRelaxSet(st).incar.get("MAGMOM", None)]

    wf = get_wf_full_hse(
        structure=st,
        task_arg={"prev_calc_dir": scf_dir},
        charge_states=[0],
        gamma_only=False,
        gamma_mesh=True,
        scf_dos=True,
        nupdowns=[-1],
        task="hse_soc",
        vasptodb={
            "category": CATEGORY,
            "NN": entry["NN"],
            "defect_entry": entry["defect_entry"],
            "nonsoc_from": "{}/{}/{}/{}".format(d_name, col_name, entry["task_label"], entry["task_id"])
        },
        wf_addition_name="{}:{}".format(entry["perturbed"], entry["task_id"]),
        category=CATEGORY,
    )

    wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)

    wf = add_modify_incar(wf, {"incar_update":{"LCHARG":False, "LWAVE":True, "LAECHG":False}})

    wf = add_modify_incar(wf)

    wf = scp_files(
        wf,
        "/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(CATEGORY),
        fw_name_constraint=wf.fws[-1].name,
    )

    # wf = clear_to_db(wf, fw_name_constraint=wf.fws[0].name)

    # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")

    wf = set_execution_options(wf, category=CATEGORY)
    wf = preserve_fworker(wf)
    LPAD.add_wf(wf)
#%%
d_name = "single_photon_emitter"
col_name = "c3v_to_ch_test"
std_defect = get_db(d_name, col_name, port=1234)

for tk_id in [241]:
    entry = std_defect.collection.find_one({"task_id": tk_id})

    st = Structure.from_dict(entry["output"]["structure"])
    print("=="*20)
    print(entry["chemsys"])
    for i in entry["NN"]:

        print("{:d}: {:.3f}".format(i, st.get_distance(entry["NN"][-1], i)))

#%% soc cdft
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from fireworks import LaunchPad
from pymatgen.io.vasp import Poscar
import os, glob
from monty.serialization import loadfn
from projects.antisiteQubit.wf_defect import ZPLWF
from qubitPack.tool_box import get_db


CATEGORY = "soc_cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))

"""
WSe2: 
    FERWE = 657*1 1*0 1*1 301*0
    NBANDS = 960
    ENCUT = 336
    NCORE = 10 /5 nodes 100 cores
    soc=/home/tug03990/work/single_photon_emitter/soc_standard_defect/block_2021-01-17-04-46-01-108825/launcher_2021-01-27-00-13-26-988037"
    
WS2:
    FERWE = 657*1 1*0 1*1 301*0
    NBANDS = 960
    ENCUT = 336
    NCORE = 10 /5 nodes 100 cores
    soc=/home/tug03990/work/single_photon_emitter/soc_standard_defect/block_2021-01-17-04-46-01-108825/launcher_2021-01-27-03-41-44-622947
    nonsoc = '/home/tug03990/work/single_photon_emitter/soc_standard_defect/block_2021-01-17-04-46-01-108825/launcher_2021-01-26-23-58-04-516424'
"""
tk_id = 268# 267te, 268s, 269se
nbands = 960

db_a = get_db("single_photon_emitter", "soc_standard_defect", port=12345)
entry_a = db_a.collection.find_one({"task_id":tk_id})
a_path = entry_a["dir_name"].split("/")[-1]

nonsoc_a_entry = entry_a["nonsoc_from"].split("/")
db_nonsoc_a = get_db(nonsoc_a_entry[0], nonsoc_a_entry[1], port=12345)
nonsoc_a_path = db_nonsoc_a.collection.find_one({"task_id":int(nonsoc_a_entry[-1])})["dir_name"].split("/")[-1]


anti_triplet = ZPLWF(os.path.join("/home/tug03990/work/single_photon_emitter/soc_standard_defect/"
                     "block_2021-01-17-04-46-01-108825/", a_path), None)

wf = anti_triplet.wf(
    "B-C-D", 0, up_occupation=anti_triplet.get_lowest_unocc_band_idx(tk_id, db_a, nbands),
    down_occupation=None, nbands=nbands, gamma_only=True, selective_dyn=None,
    nonsoc_prev_dir=os.path.join("/home/tug03990/work/single_photon_emitter/standard_defect/"
                                 "block_2021-01-12-17-28-08-335869", nonsoc_a_path)
)

wf = jmodify_to_soc(wf, nbands=nbands, structure=anti_triplet.structure)

wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

wf = add_modify_incar(wf)

wf = set_execution_options(wf, category=CATEGORY)
wf = preserve_fworker(wf)
LPAD.add_wf(wf)


