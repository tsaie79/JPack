#%% make standard defect calc
from fireworks import LaunchPad
from qubitPack.tool_box import *
from projects.antisiteQubit.wf_defect import get_wf_full_hse
from atomate.vasp.jpowerups import *
from atomate.vasp.powerups import *


CATEGORY = "standard_defect"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))

def ws2_anion_antisite(defect_type="substitutions", dtort=0.01): #1e-4
    db_name, col_name = "single_photon_emitter", "pc"
    col = get_db(db_name, col_name, port=12345).collection

    mx2s = col.find({"task_id":{"$in":[188, 187]}})
    geo_spec = {5* 5 * 3: [20]}
    for mx2 in mx2s:
        pc = Structure.from_dict(mx2["output"]["structure"])
        defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
        cation, anion = find_cation_anion(pc)

        for sub in range(len(defect[defect_type])):
            print(cation, anion)
            # cation vacancy
            if "as_1_{}_on_{}".format(cation,anion) not in defect[defect_type][sub]["name"]:
                continue
            for na, thicks in geo_spec.items():
                for thick in thicks:
                    se_antisite = GenDefect(
                        orig_st=pc,
                        defect_type=(defect_type, sub),
                        natom=na,
                        vacuum_thickness=thick,
                        distort=dtort,
                        sub_on_side=None
                    )
                    wf = get_wf_full_hse(
                        structure=se_antisite.defect_st,
                        task_arg=None,
                        charge_states=[0],
                        gamma_only=False,
                        gamma_mesh=True,
                        scf_dos=False,
                        nupdowns=[-1],
                        task="opt-hse_relax-hse_scf",
                        vasptodb={
                            "category": CATEGORY, "NN": se_antisite.NN,
                            "defect_entry": se_antisite.defect_entry,
                            "lattice_constant": "HSE",
                            "perturbed": se_antisite.distort,
                            "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
                                  },
                        wf_addition_name="{}:{}".format(na, thick),
                        category=CATEGORY,
                    )

                    # wf = scp_files(
                    #     wf,
                    #     "/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(CATEGORY),
                    #     fw_name_constraint=wf.fws[-1].name,
                    # )
                    # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb",
                    #                     fw_name_constraint="HSE_scf")

                    wf = clear_to_db(wf, fw_name_constraint=wf.fws[0].name)
                    wf = clear_to_db(wf, fw_name_constraint=wf.fws[1].name)

                    wf = set_queue_options(wf, "24:00:00", fw_name_constraint="PBE_relax")
                    wf = set_queue_options(wf, "32:00:00", fw_name_constraint="HSE_relax")
                    wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")
                    # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_soc")
                    wf = add_modify_incar(wf)
                    # related to directory
                    wf = preserve_fworker(wf)
                    wf.name = wf.name+":dx[{}]".format(se_antisite.distort)
                    print(wf.name)
                    LPAD.add_wf(wf)
ws2_anion_antisite()

#%% SOC STANDARD
from projects.antisiteQubit.wf_defect import get_wf_full_hse
from fireworks import LaunchPad
import os
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from pymatgen.io.vasp.sets import MPRelaxSet
from qubitPack.tool_box import get_db

CATEGORY = "soc_standard_defect"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))

tk_id = 22
src = get_db("single_photon_emitter", "standard_defect", port=12345).collection.find_one({"task_id":tk_id})

scf_dir = src["calcs_reversed"][0]["dir_name"]

st = Structure.from_dict(src["input"]["structure"])

magmom = [[0,0,mag_z] for mag_z in MPRelaxSet(st).incar.get("MAGMOM", None)]

wf = get_wf_full_hse(
    structure=st,
    task_arg={"prev_calc_dir": scf_dir},
    charge_states=[0],
    gamma_only=False,
    gamma_mesh=True,
    scf_dos=False,
    nupdowns=[-1],
    task="hse_scf-hse_soc",
    vasptodb={"category": CATEGORY, "perturbed": src["perturbed"],
              "pc_from": "single_photon_emitter/standard_defect/{}".format(tk_id)},
    wf_addition_name="soc_Ch:pert_{}".format(src["perturbed"]),
    category=CATEGORY,
)

wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)

wf = add_modify_incar(wf)

wf = clear_to_db(wf, fw_name_constraint="HSE_scf")

wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb", fw_name_constraint=wf.fws[-1].name)

wf = set_execution_options(wf, category=CATEGORY)
wf = preserve_fworker(wf)


LPAD.add_wf(wf)

#%% make interpolated structures
from qubitPack.tool_box import get_db, get_interpolate_sts
from pymatgen import Structure
import numpy as np
from projects.antisiteQubit.wf_defect import get_wf_full_hse
from fireworks import LaunchPad
import os
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from pymatgen.io.vasp.sets import MPRelaxSet

db = get_db("single_photon_emitter", "standard_defect", port=12345)

init_task = 22
final_task = 20
init_st = Structure.from_dict(db.collection.find_one({"task_id":init_task})["output"]["structure"])
final_st = Structure.from_dict(db.collection.find_one({"task_id":final_task})["output"]["structure"])

a = get_interpolate_sts(init_st, final_st, disp_range=np.linspace(-2, 2, 21), output_dir=None)
print(list(a[0].keys()))

CATEGORY = "c3v_to_ch_test"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))

for dr, st in a[0].items():
    wf = get_wf_full_hse(
        structure=st,
        task_arg=None,
        charge_states=[0],
        gamma_only=False,
        gamma_mesh=True,
        scf_dos=False,
        nupdowns=[-1],
        task="hse_scf",
        vasptodb={"category": CATEGORY, "perturbed": dr,
                  "init_st": "single_photon_emitter/standard_defect/{}".format(init_task),
                  "final_st": "single_photon_emitter/standard_defect/{}".format(final_task),
                  "dr_info": a[1],
                  "ps": "strict_ch"
                  },
        wf_addition_name="Ch:pert_{}".format(dr),
        category=CATEGORY,
    )

    wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)

    wf = add_modify_incar(wf)

    wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb", fw_name_constraint=wf.fws[-1].name)

    wf = set_execution_options(wf, category=CATEGORY)
    wf = preserve_fworker(wf)


    LPAD.add_wf(wf)

#%% check energy for different pertubation
from pymatgen import Structure
from qubitPack.tool_box import get_db
import pandas as pd
from matplotlib import pyplot as plt

db = get_db("single_photon_emitter", "c3v_to_ch_test")

se = db.collection.aggregate(
    [
        {"$match": {"chemsys":"Se-W"}},
        {"$project":{"ch": "$chemsys","pert":"$perturbed", "energy":"$output.energy", "_id":0, "task_id":1, "r": "$dr_info.Delta_R"}}
    ]
)
se = list(se)

s = db.collection.aggregate(
    [
        {"$match": {"chemsys":"S-W", "ps": {"$exists": False}}},
        {"$project":{"ch": "$chemsys","pert":"$perturbed", "energy":"$output.energy", "_id":0, "task_id":1, "r": "$dr_info.Delta_R"}}
    ]
)
s = list(s)

E1 = db.collection.find_one({"chemsys":"Se-W", "perturbed":0})["output"]["energy"]
E0 = db.collection.find_one({"chemsys":"S-W", "perturbed":0, "ps": {"$exists": False}})["output"]["energy"]
dE = E1-E0

for i in se:
    print(i["task_id"])
    Structure.from_dict(db.collection.find_one({"task_id":i["task_id"]})["output"]["structure"]).to("poscar", "/Users/jeng-yuantsai/Research/project/"
                                                                                                              "single_photon/calculations/c3v_to_ch_test/structures/WSe2/{}.vasp".format(i["pert"]))

for i in s:
    print(i["task_id"])
    Structure.from_dict(db.collection.find_one({"task_id":i["task_id"]})["output"]["structure"]).to("poscar", "/Users/jeng-yuantsai/Research/project/"
                                                                                                              "single_photon/calculations/c3v_to_ch_test/structures/WS2/{}.vasp".format(i["pert"]))

figure, ax = plt.subplots()

for e in [se, s]:
    e = pd.DataFrame(e)
    print(e[e.energy == e.energy.min()])
    if e["ch"][0] == "Se-W":
        e["energy"] -= dE
    e["energy"] -= E0
    e["nr"] = e["pert"] * e["r"] * 0.529
    e.sort_values("pert", inplace=True)
    ax.plot(e["nr"], e["energy"], label=e["ch"][0],)
    ax.legend()
    ax.set_xlabel("displacement (Å)")
    ax.set_ylabel("relative energy (eV)")
    ax.set_title("C3v to Ch")

figure.show()

#%% check energy for different pertubation
from pymatgen import Structure
from qubitPack.tool_box import get_db
import pandas as pd
from matplotlib import pyplot as plt

db = get_db("single_photon_emitter", "c3v_to_ch_test")

se = db.collection.aggregate(
    [
        {"$match": {"chemsys":"S-W", "ps": {"$exists": False}}},
        {"$project":{"ps":"easy","ch": "$chemsys","pert":"$perturbed", "energy":"$output.energy", "_id":0, "task_id":1, "r": "$dr_info.Delta_R"}}
    ]
)
se = list(se)

s = db.collection.aggregate(
    [
        {"$match": {"chemsys":"S-W", "ps": "strict_ch"}},
        {"$project":{"ps": "strict_ch", "ch": "$chemsys","pert":"$perturbed", "energy":"$output.energy", "_id":0, "task_id":1, "r": "$dr_info.Delta_R"}}
    ]
)
s = list(s)

for i in se:
    print(i["task_id"])
    Structure.from_dict(db.collection.find_one({"task_id":i["task_id"]})["output"]["structure"]).to("poscar", "/Users/jeng-yuantsai/Research/project/"
                                                                                                              "single_photon/calculations/c3v_to_ch_test/structures/WS2/ch/{}.vasp".format("_".join(str(i["pert"]).split("."))))

for i in s:
    print(i["task_id"])
    Structure.from_dict(db.collection.find_one({"task_id":i["task_id"]})["output"]["structure"]).to("poscar", "/Users/jeng-yuantsai/Research/project/"
                                                                                                              "single_photon/calculations/c3v_to_ch_test/structures/WS2/strict_ch/{}.vasp".format("_".join(str(i["pert"]).split("."))))

figure, ax = plt.subplots()

for e in [se, s]:
    e = pd.DataFrame(e)
    print(e[e.energy == e.energy.min()])
    if e["ps"][0] == "strict_ch":
        e["nr"] = e["pert"] * e["r"] * 0.529
    else:
        e["nr"] = e["pert"] * e["r"] * 0.529
    e.sort_values("pert", inplace=True)
    ax.plot(e["nr"], e["energy"], label=e["ps"][0])
    ax.legend()
    ax.set_xlabel("displacement (Å)")
    ax.set_ylabel("relative energy (eV)")
    ax.set_title("C3v to Ch")

figure.show()

#%%
from pymatgen import Structure
from qubitPack.tool_box import get_db
import pandas as pd
from matplotlib import pyplot as plt
from qubitPack.qc_searching.analysis.main import get_defect_state

db = get_db("single_photon_emitter", "c3v_to_ch_test")

e = db.collection.find_one({"chemsys": "S-W", "task_id":241, "task_label":"HSE_scf"})
Structure.from_dict(e["output"]["structure"]).to("poscar", "/Users/jeng-yuantsai/Research/project/"
                                                           "single_photon/calculations/c3v_to_ch_test/structures/WS2/241.vasp")
print(e["task_id"])

# tot, proj, d_df = get_defect_state(
#     db,
#     {"task_id": 189},
#     1,-3,
#     None,
#     False,
#     "dist",
#     None,
#     0.1
# )