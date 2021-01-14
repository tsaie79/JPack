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
from pymatgen.io.vasp.inputs import Poscar
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
    # source with encut=336
    # anti_triplet = ZPLWF('/gpfs/work/tug03990/antisiteQubit/mx2s_up_encut/'
    #                      'block_2020-12-01-03-21-58-237215/launcher_2020-12-02-13-25-41-348623', 'triplet')
    # occ = ["302*1.0 1*0.5 1*0.5 1*1.0 115*0.0", "302*1.0 118*0.0"]


    #efrc
    anti_triplet = ZPLWF("/home/tug03990/work/mx2_antisite_basic_aexx0.25_sigma_test/"
                       "block_2020-11-16-00-55-05-740824/launcher_2020-11-23-05-01-59-204100", "triplet")
    occ = ["302*1.0 1*0.5 1*0.5 1*1.0 215*0", "302*1.0 218*0.0"]

    # anti_triplet = ZPLWF('/home/tug03990/work/mx2_antisite_basic_aexx0.25_cdft/'
    #               'block_2020-07-27-19-52-24-555470/launcher_2020-11-29-20-24-21-823909', "triplet")
    #occ = ["303*1.0 1*0 1*1.0 175*0","302*1.0 178*0.0"]

    # this is the source with sigma=0.05
    # anti_triplet = ZPLWF("/home/tug03990/work/mx2_antisite_basic_aexx0.25_final/"
    #                      "block_2020-07-03-20-26-36-261691/launcher_2020-07-04-02-43-56-914816", "triplet")
    # occ = ["302*1.0 1*0.5 1*0.5 1*1.0 175*0.0", "302*1.0 178*0.0"]

    sd_sites = anti_triplet.set_selective_sites(25, 5)
    # print(len(moving_sites))

    def wf(potim, ibrion, ediffg, ediff, iopt, maxmove, moving_sites):
        wf = anti_triplet.wf(
            "C", 0, up_occupation=occ[0], down_occupation=occ[1],
            nbands=520, gamma_only=True, selective_dyn=moving_sites,
            specific_structure=Poscar.from_file(
                '/home/tug03990/work/antisiteQubit/mos2_cdft_div/block_2020-11-30-17-05-10-860374/'
                'launcher_2020-12-08-05-21-14-092530/mv_110.vasp'
            )
        )

        wf = add_modify_incar(wf, {"incar_update": {"ENCUT":320, "LASPH":True, "SIGMA":0.002}})

        wf = add_modify_incar(wf, {"incar_update": {
            # "IOPT": iopt, "MAXMOVE": maxmove,
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

    for po in [0.01562]:
        for ib in [1]:
            for sd in [[]]:
                wf(po,ib,-0.02, 1E-4, 7, 0.5, sd)

    # for po in [0]:
    #     for ib in [3]:
    #         for sd in [[]]:
    #             for iopt in [7]:
    #                 for maxmv in [1E-3]:
    #                     wf(po,ib,-0.02, 1E-4, iopt, maxmv, sd)


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
    # anti_triplet = ZPLWF("/home/tug03990/work/antisiteQubit/perturbed/"
    #                      "block_2020-11-21-04-37-28-746929/launcher_2020-11-23-22-19-31-039641", "triplet")
    anti_triplet = ZPLWF("/home/jengyuantsai/Research/projects/single_photon_emitter/standard_defect/"
                         "launcher_2021-01-13-04-09-57-563577", "triplet")

    sd_sites = anti_triplet.set_selective_sites(25, 5)
    sd_sites.sort()
    print(sd_sites)

    def wf(potim, ibrion, ediffg, ediff, moving_sites):
        wf = anti_triplet.wf(
            "B", 0, up_occupation="302*1.0 1*0.5 1*0.5 1*1.0 155*0.0",#"303*1.0 1*0 1*1.0 175*0",
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
        print(wf.name)
        LPAD.add_wf(wf)

    for po in [0]:
        for ib in [3]:
            for sd in [sd_sites]:
                wf(po,ib,-0.02, 1E-4, sd)


antistie_triplet_ZPL()

#%%
from atomate.vasp.powerups import *
from fireworks import LaunchPad
from projects.antisiteQubit.wf_defect import ZPLWF
import os
from glob import glob
from pymatgen.io.vasp.inputs import Incar
import pandas as pd

CATEGORY = "mos2_cdft_div"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/antisiteQubit/{}/my_launchpad.yaml".format(CATEGORY)))

fw_ids = LPAD.get_fw_ids({"spec._category":"mos2_cdft_div", "name":{"$regex":"CDFT-C-HSE"},
                          "state":{"$in":["RUNNING","COMPLETED", "FIZZLED"]}})
data = []
for fw_id in fw_ids:
    fw = LPAD.get_fw_by_id(fw_id)
    fworker = fw.launches[-1].fworker.name
    state = fw.launches[-1].state
    cat = fw.spec["_category"]
    src = fw.spec["_tasks"][1].as_dict()["calc_dir"]
    if fworker == "efrc":
        dir_name = LPAD.get_launchdir(fw_id)
        if dir_name:
            e = {}

            e["state"] = state
            e["cat"] = cat
            e["fw_id"] = fw_id

            incar_path = glob(os.path.join(dir_name, "INCAR*"))[0]
            incar = Incar.from_file(incar_path).as_dict()

            e["IBRION"] = incar.get("IBRION")
            e["POTIM"] = incar.get("POTIM")
            e["IOPT"] = incar.get("IOPT", None)
            e["MAXMOVE"] = incar.get("MAXMOVE", None)

            e["ENCUT"] = incar.get("ENCUT")
            e["SIGMA"] = incar.get("SIGMA")
            e["EDIFF"] = incar.get("EDIFF")
            e["EDIFFG"] = incar.get("EDIFFG")
            e["LASPH"] = incar.get("LASPH")
            e["ALGO"] = incar.get("ALGO")
            e["NBANDS"] = incar.get("NBANDS")

            e["src"] = src
            e["hope"] = None
            data.append(e)
        else:
            continue

df = pd.DataFrame(data)
df = df.set_index("fw_id")
for idx in [1102, 1101, 1076]:
    df.loc[idx, "hope"] = True
df.loc[1222, "hope"] = "ST USED"
df.loc[1173, "hope"] = "FINAL"
# df.to_csv("~/efrc.csv")

#%%
"""
get displacements of antisites
"""
from pymatgen.io.vasp.outputs import Vasprun
from fireworks import LaunchPad
import os
import pandas as pd
from matplotlib import pyplot as plt


CATEGORY = "mos2_cdft_div"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/antisiteQubit/{}/my_launchpad.yaml".format(CATEGORY)))

fw_id = 1222
p = LPAD.get_launchdir(fw_id)

vr = Vasprun(os.path.join(p,"vasprun.xml.gz"))
sts = vr.ionic_steps

NN = [5,6,0,25]

x0, y0, z0 = sts[0]["structure"].get_distance(25, 5), sts[0]["structure"].get_distance(25, 6), \
             sts[0]["structure"].get_distance(25, 0)
print(x0, y0, z0)
data = []
for st in sts:
    e = {}
    en = st["e_fr_energy"]
    force = st["forces"][25]
    pos = []
    for nn in [5,6,0]:
        d = st["structure"].get_distance(25, nn)
        pos.append(d)
    print(pos)
    x, y, z = pos[0], pos[1], pos[2]
    e["x"] = x - x0
    e["y"] = y - y0
    e["z"] = z - z0
    print(e)
    e["Fx"] = force[0]
    e["Fy"] = force[1]
    e["Fz"] = force[2]
    e["E0"] = en
    data.append(e)
    print(data)
df = pd.DataFrame(data)

fig, axs = plt.subplots(5, sharex=True)
axs[0].plot(df["x"], label="d5")
axs[1].plot(df["y"], label="d6")
axs[2].plot(df["z"], label="d0")


# axs[3].plot(df["Fx"], label="Fx")
# axs[4].plot(df["Fy"], label)
axs[3].plot(df["Fz"], label="Fz")

axs[4].plot(df["E0"], label="E0")

for i in range(5):
    axs[i].legend()

plt.show()