#%% soc cdft B parabola
import glob
from fireworks import LaunchPad
import os
from projects.antisiteQubit.wf_defect import ZPLWF
from monty.serialization import loadfn
from pymatgen.io.vasp.inputs import Poscar
from pymatgen import Structure
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from qubitPack.tool_box import get_db, get_interpolate_sts
from pymatgen import Structure
import numpy as np

CATEGORY = "soc_cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))


tk_id = 268# 267te, 268s, 269se
nbands = 960
soc_p = "/home/tug03990/scratch/single_photon_emitter/soc_standard_defect/block_2021-01-22-20-56-50-305107"
nonsoc_p = "/home/tug03990/scratch/single_photon_emitter/standard_defect/block_2021-01-11-17-39-35-861559"


db_a = get_db("single_photon_emitter", "soc_standard_defect", port=12345)
entry_a = db_a.collection.find_one({"task_id":tk_id})
a_path = entry_a["dir_name"].split("/")[-1]

nonsoc_a_entry = entry_a["nonsoc_from"].split("/")
db_nonsoc_a = get_db(nonsoc_a_entry[0], nonsoc_a_entry[1], port=12345)
nonsoc_a_path = db_nonsoc_a.collection.find_one({"task_id":int(nonsoc_a_entry[-1])})["dir_name"].split("/")[-1]


soc_cdft = get_db("single_photon_emitter", "soc_cdft", port=12345)
b_st = Structure.from_dict(soc_cdft.collection.find_one(
    {"chemsys":entry_a["chemsys"], "task_label": {"$regex": "B-HSE_scf_soc"}})["input"]["structure"])
c_st = Structure.from_dict(soc_cdft.collection.find_one(
    {"chemsys":entry_a["chemsys"], "task_label": {"$regex": "C-HSE_relax_soc"}})["output"]["structure"])

dx = np.linspace(0, 2,11)
sts, info = get_interpolate_sts(
    b_st,
    c_st,
    disp_range=dx,
    output_dir=None
)


anti_triplet = ZPLWF(os.path.join(soc_p, a_path), None)

for frac, st in sts.items():
    wf = anti_triplet.wf(
        "B", 0, up_occupation=anti_triplet.get_lowest_unocc_band_idx(tk_id, db_a, nbands),
        down_occupation=None, nbands=nbands, gamma_only=True, selective_dyn=None, specific_structure=Poscar(st),
        nonsoc_prev_dir=os.path.join(nonsoc_p, nonsoc_a_path)
    )


    wf = jmodify_to_soc(wf, nbands=nbands, structure=anti_triplet.structure)

    wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)

    wf = add_modify_incar(wf)

    wf.name += wf.name+"{}".format(frac)
    wf = add_additional_fields_to_taskdocs(wf, {"info":info, "dx": frac})
    wf = set_execution_options(wf, category=CATEGORY)
    wf = preserve_fworker(wf)
    LPAD.add_wf(wf)

#%% soc cdft D parabola
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


CATEGORY = "soc_cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))


tk_id = 268# 267te, 268s, 269se
nbands = 960
soc_p = "/home/tug03990/scratch/single_photon_emitter/soc_standard_defect/block_2021-01-22-20-56-50-305107"
nonsoc_p = "/home/tug03990/scratch/single_photon_emitter/standard_defect/block_2021-01-11-17-39-35-861559"


db_a = get_db("single_photon_emitter", "soc_standard_defect", port=12345)
entry_a = db_a.collection.find_one({"task_id":tk_id})
a_path = entry_a["dir_name"].split("/")[-1]

nonsoc_a_entry = entry_a["nonsoc_from"].split("/")
db_nonsoc_a = get_db(nonsoc_a_entry[0], nonsoc_a_entry[1], port=12345)
nonsoc_a_path = db_nonsoc_a.collection.find_one({"task_id":int(nonsoc_a_entry[-1])})["dir_name"].split("/")[-1]


soc_cdft = get_db("single_photon_emitter", "soc_cdft", port=12345)
b_st = Structure.from_dict(soc_cdft.collection.find_one(
    {"chemsys":entry_a["chemsys"], "task_label": {"$regex": "B-HSE_scf_soc"}})["input"]["structure"])
c_st = Structure.from_dict(soc_cdft.collection.find_one(
    {"chemsys":entry_a["chemsys"], "task_label": {"$regex": "C-HSE_relax_soc"}})["output"]["structure"])

dx = np.linspace(-1, 1,11)
sts, info = get_interpolate_sts(
    b_st,
    c_st,
    disp_range=dx,
    output_dir=None
)

for frac, st in sts.items():

    wf = get_wf_full_hse(
        structure=st,
        task_arg={"prev_calc_dir": os.path.join(nonsoc_p, nonsoc_a_path)},
        charge_states=[0],
        gamma_only=False,
        gamma_mesh=True,
        scf_dos=False,
        nupdowns=[-1],
        task="hse_soc",
        vasptodb={
            "task_label": "{}-CDFT-D-HSE_scf_soc".format(entry_a["formula_pretty"])
        },
        wf_addition_name="D-{}".format(frac),
        category=CATEGORY,
    )

    wf = write_PMGObjects(wf, dict(poscar=Poscar(st)))

    wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)

    wf = add_modify_incar(wf, {"incar_update":{"LCHARG":False, "LWAVE":False, "LAECHG":False}})

    wf = add_modify_incar(wf)

    wf = add_additional_fields_to_taskdocs(wf, {"info":info, "dx": frac})

    wf = set_execution_options(wf, category=CATEGORY)
    wf = preserve_fworker(wf)
    # LPAD.add_wf(wf)