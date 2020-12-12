#%% antisiteQubit/perturbed
"""
Read defect state from antisiteQubit/perturbed
"""
import os
from qubitPack.qc_searching.analysis.main import get_defect_state

proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/perturbed"
db_json = os.path.join(proj_path, "db.json")

host_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_pc"
db_host_json = os.path.join(host_path, "db.json")

perturb = 0.005
chem = "S-W"
tot, proj, d_df = get_defect_state(
    db_json,
    {"task_label":"HSE_scf", "chemsys":chem, "perturbed":perturb},
    5,-5,
    None,#os.path.join(proj_path, chem),
    True,
    "dist",
    db_host_json,
    0.05
)
tot.to_clipboard()

#%% antisiteQubit/perturbed
"""
Look at database 
"""
from pymatgen import Structure
from atomate.vasp.database import VaspCalcDb
import pandas as pd
import os

proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/perturbed"
db_json = os.path.join(proj_path, "db.json")

db = VaspCalcDb.from_db_file(db_json)


es = db.collection.aggregate(
    [
        {"$match": {"task_label":"HSE_scf", "chemsys":"Mo-S"}},
        {"$group": {"_id":"$output.spacegroup.point_group",
                    "en": {"$push": "$output.energy"},
                    "tid": {"$push": "$task_id"},
                    "defect": {"$push": "$defect_entry.name"},
                    "chem": {"$push": "$chemsys"},
                    "mag": {"$push": "$calcs_reversed.output.outcar.total_magnetization"}
                    }
         },
    ]
)

df = pd.DataFrame(es)
df = df.transpose()

e = db.collection.find_one({"task_id":481})
init = Structure.from_dict(e["input"]["structure"])
final = Structure.from_dict(e["output"]["structure"])
init.to("poscar", proj_path+"/{}_{}_init.vasp".format(e["formula_pretty"], e["task_id"]))
final.to("poscar", proj_path+"/{}_{}_final.vasp".format(e["formula_pretty"], e["task_id"]))
for idx in e["NN"][:-1]:
    print(idx, init.get_distance(25, idx), final.get_distance(25,idx))
print(e["output"]["energy"])

#%%
from atomate.vasp.powerups import *
from fireworks import LaunchPad
from projects.antisiteQubit.wf_defect import ZPLWF
import os
from glob import glob
from pymatgen.io.vasp.inputs import Incar
import pandas as pd

CATEGORY = "perturbed"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/antisiteQubit/{}/my_launchpad.yaml".format(CATEGORY)))

fw_ids = LPAD.get_fw_ids({"spec._category":CATEGORY, "state":{"$in":["RUNNING"]}})
data = []
for fw_id in fw_ids:
    fw = LPAD.get_fw_by_id(fw_id)
    fworker = fw.launches[-1].fworker.name
    state = fw.launches[-1].state
    cat = fw.spec["_category"]
    # src = fw.spec["_tasks"][1].as_dict()["calc_dir"]
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

            data.append(e)
        else:
            continue

df = pd.DataFrame(data)
for idx in [1102, 1101, 1076]:
    df.loc[df["fw_id"] == idx, ["hope"]] = True
df.loc[df["fw_id"] == 1173, ["hope"]] = "FINAL"