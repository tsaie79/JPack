#%% owls/mx2_antisite_basic_aexx0.25_cdft
"""
ZPL of MoS2
"""
from atomate.vasp.database import VaspCalcDb
import pandas as pd

db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/project/"
                             "qubit/calculations/mx2_antisite_basic_aexx0.25_cdft/db.json")

r = db.collection.aggregate(
    [
        {"$match": {"chemsys":"Mo-S"}},
        {"$project": {
                "_id":0,
                "A": "$source.total_energy",
                "src_dir": "$source.prev_path",
                "tk": "$task_label",
                "en": "$output.energy",
                "tid": "$task_id"
            }
         },
        {"$group":
             {"_id": "$tk",
              "tid": {"$push": "$tid"},
              # "src_dir": {"$push": "$src_dir"},
              "en": {"$push": "$en"},
              "A": {"$push": "$A"},

              }
         }
    ]
)

df = pd.DataFrame(list(r)).set_index("_id")

data = {}
idx = -1
data["A-B"] = [df.loc["CDFT-B-HSE_scf"]["en"][idx] - df.loc["CDFT-B-HSE_scf"]["A"][idx]]
data["B-C"] = [df.loc["CDFT-B-HSE_scf"]["en"][idx] - df.loc["CDFT-C-HSE_relax"]["en"][idx]]
data["C-D"] = [df.loc["CDFT-C-HSE_relax"]["en"][idx] - df.loc["CDFT-D-HSE_scf"]["en"][idx]]
data["D-A"] = [df.loc["CDFT-D-HSE_scf"]["en"][idx] - df.loc["CDFT-B-HSE_scf"]["A"][idx]]

result = pd.DataFrame(data)

#%% mx2_antisite_basic_aexx0.25_sigma_test
"""
Run CDFT on MoS2
"""
from atomate.vasp.powerups import *
from fireworks import LaunchPad
from projects.antisiteQubit.wf_defect import ZPLWF
import os

CATEGORY = "mx2_antisite_basic_aexx0.25_cdft"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/category/{}/my_launchpad.yaml".format(CATEGORY)))


def antistie_triplet_ZPL(): #6.5
    # MoS2
    anti_triplet = ZPLWF("/home/tug03990/work/mx2_antisite_basic_aexx0.25_sigma_test/"
                       "block_2020-11-16-00-55-05-740824/launcher_2020-11-23-05-01-59-204100", "triplet")

    moving_sites = anti_triplet.set_selective_sites(25, 5)
    print(len(moving_sites))

    wf = anti_triplet.wf(
        "C", 0, up_occupation="302*1.0 1*0.5 1*0.5 1*1.0 175*0", #"303*1.0 1*0 1*1.0 175*0"
        down_occupation="302*1.0 178*0.0", nbands=480, gamma_only=True, selective_dyn=None
    )
    wf = add_modify_incar(wf, {"incar_update": {"POTIM":0.125, "IBRION":1, "EDIFFG":-0.02, "EDIFF":1E-4}},
                          fw_name_constraint=wf.fws[0].name)
    # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
    wf = set_queue_options(wf, "48:00:00", fw_name_constraint=wf.fws[0].name)
    # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)
    wf = set_execution_options(wf, category=CATEGORY)
    wf = add_modify_incar(wf)
    wf = preserve_fworker(wf)
    wf = clean_up_files(wf, files=["WAVECAR*"], task_name_constraint="VaspToDb")
    wf = add_additional_fields_to_taskdocs(wf, {"defect_type": "TIM0.125:IB1"})
    wf.name = wf.name + ":delta{:.2f}".format(len(moving_sites) / len(anti_triplet.structure.sites))
    wf.name = wf.name + ":TIM0.125:IB1"
    LPAD.add_wf(wf)
    print(wf.name)

antistie_triplet_ZPL()
