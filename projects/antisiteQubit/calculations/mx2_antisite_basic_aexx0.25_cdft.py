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
idx = -2
data["A-B"] = [df.loc["CDFT-B-HSE_scf"]["en"][idx] - df.loc["CDFT-B-HSE_scf"]["A"][idx]]
data["B-C"] = [df.loc["CDFT-B-HSE_scf"]["en"][idx] - df.loc["CDFT-C-HSE_relax"]["en"][idx]]
data["C-D"] = [df.loc["CDFT-C-HSE_relax"]["en"][idx] - df.loc["CDFT-D-HSE_scf"]["en"][idx]]
data["D-A"] = [df.loc["CDFT-D-HSE_scf"]["en"][idx] - df.loc["CDFT-B-HSE_scf"]["A"][idx]]

result = pd.DataFrame(data)

#%% mx2_antisite_basic_aexx0.25_sigma_test
"""
Run CDFT on MoS2 w/o perturbation
"""
from atomate.vasp.powerups import *
from fireworks import LaunchPad
from projects.antisiteQubit.wf_defect import ZPLWF
import os

CATEGORY = "mos2_cdft_div"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/antisiteQubit/{}/my_launchpad.yaml".format(CATEGORY)))


def antistie_triplet_ZPL(): #6.5
    # MoS2
    #owls
    # anti_triplet = ZPLWF("/home/tug03990/work/mx2_antisite_basic_aexx0.25_cdft/"
    #                      "block_2020-11-16-00-55-05-740824/launcher_2020-11-23-05-01-59-204100", "triplet")
    #efrc
    anti_triplet = ZPLWF("/home/tug03990/work/mx2_antisite_basic_aexx0.25_sigma_test/"
                       "block_2020-11-16-00-55-05-740824/launcher_2020-11-23-05-01-59-204100", "triplet")

    # anti_triplet = ZPLWF('/home/tug03990/work/mx2_antisite_basic_aexx0.25_cdft/'
    #               'block_2020-07-27-19-52-24-555470/launcher_2020-11-29-20-24-21-823909', "triplet")

    sd_sites = anti_triplet.set_selective_sites(25, 3)
    # print(len(moving_sites))

    def wf(potim, ibrion, ediffg, ediff, iopt,maxmove, moving_sites):
        wf = anti_triplet.wf(
            "C", 0, up_occupation="302*1.0 1*0.5 1*0.5 1*1.0 175*0",#"303*1.0 1*0 1*1.0 175*0",
            down_occupation="302*1.0 178*0.0", nbands=480, gamma_only=True, selective_dyn=moving_sites
        )


        wf = add_modify_incar(wf, {"incar_update": {"ENCUT":320, "LASPH":True}})

        wf = add_modify_incar(wf, {"incar_update": {
            "IOPT": iopt, "MAXMOVE": maxmove,
            "POTIM": potim, "IBRION": ibrion, "EDIFFG": ediffg, "EDIFF": ediff}}, fw_name_constraint=wf.fws[0].name)

        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

        wf = add_modify_incar(wf)
        wf = clean_up_files(wf, files=["WAVECAR*"], task_name_constraint="VaspToDb")
        wf = add_additional_fields_to_taskdocs(wf, {"defect_type": "TIM{}:IB{}:IOPT{}".format(potim,ibrion,iopt),
                                                    "selective_dynamics":moving_sites})

        wf = set_execution_options(wf, category=CATEGORY)
        wf = preserve_fworker(wf)

        if moving_sites:
            wf.name = wf.name + ":delta{:.2f}".format(len(moving_sites) / len(anti_triplet.structure.sites))
        wf.name = wf.name + ":TIM{}:IB{}:IOPT{}:{}".format(potim, ibrion,iopt, maxmove)
        LPAD.add_wf(wf)
        print(wf.name)

    # for po in [0.015625]:
    #     for ib in [1]:
    #         for sd in [[]]:
    #             wf(po,ib,-0.01, 1E-4, sd)
    for po in [0]:
        for ib in [3]:
            for sd in [sd_sites]:
                for iopt in [7]:
                    for maxmv in [1E-4]:
                        wf(po,ib,-0.02, 1E-4, iopt, maxmv, sd)


antistie_triplet_ZPL()
#%% antisiteQubit/perturbed
"""
Run CDFT on MoS2 w/ perturbation
"""
from atomate.vasp.powerups import *
from fireworks import LaunchPad
from projects.antisiteQubit.wf_defect import ZPLWF
import os

CATEGORY = "Mos_Ch_CDFT"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/antisiteQubit/{}/my_launchpad.yaml".format(CATEGORY)))


def antistie_triplet_ZPL(): #6.5
    #efrc
    anti_triplet = ZPLWF("/home/tug03990/work/antisiteQubit/perturbed/"
                         "block_2020-11-21-04-37-28-746929/launcher_2020-11-23-22-19-31-039641", "triplet")

    sd_sites = anti_triplet.set_selective_sites(25, 5)
    sd_sites.sort()
    print(sd_sites)

    def wf(potim, ibrion, ediffg, ediff, moving_sites):
        wf = anti_triplet.wf(
            "B-C-D", 0, up_occupation="302*1.0 1*0.5 1*0.5 1*1.0 155*0.0",#"303*1.0 1*0 1*1.0 175*0",
            down_occupation="302*1.0 158*0.0", nbands=460, gamma_only=True, selective_dyn=moving_sites
        )

        wf = add_modify_incar(wf, {"incar_update": {"ENCUT": 320, "LASPH": True}})

        wf = add_modify_incar(wf, {"incar_update": {
            "IOPT": 7, "MAXMOVE": 5E-4,
            "POTIM": potim, "IBRION": ibrion, "EDIFFG": ediffg, "EDIFF": ediff}}, fw_name_constraint=wf.fws[1].name)

        wf = add_modify_incar(wf, {"incar_update": {"EDIFF":1E-5}}, fw_name_constraint=wf.fws[0].name)

        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

        wf = add_modify_incar(wf)
        wf = clean_up_files(wf, files=["WAVECAR*"], task_name_constraint="VaspToDb")
        wf = add_additional_fields_to_taskdocs(wf, {"defect_type": "TIM{}:IB{}".format(potim,ibrion),
                                                    "selective_dynamics":moving_sites, "perturbed":True})

        wf = set_execution_options(wf, category=CATEGORY)
        wf = preserve_fworker(wf)

        if moving_sites:
            wf.name = wf.name + ":delta{:.2f}".format(len(moving_sites) / len(anti_triplet.structure.sites))
        wf.name = wf.name + ":TIM{}:IB{}".format(potim, ibrion)
        # LPAD.add_wf(wf)
        print(wf.name)

    for po in [0]:
        for ib in [3]:
            for sd in [sd_sites]:
                wf(po,ib,-0.02, 1E-4, sd)


antistie_triplet_ZPL()

#%%
from projects.antisiteQubit.wf_defect import ZPLWF

anti_triplet = ZPLWF("/home/tug03990/work/antisiteQubit/perturbed/"
                     "block_2020-11-21-04-37-28-746929/launcher_2020-11-23-22-19-31-039641", "triplet")

sd_sites = anti_triplet.set_selective_sites(25, 3)
print(sd_sites)