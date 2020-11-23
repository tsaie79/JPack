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

proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/perturbed"
db_json = os.path.join(proj_path, "db.json")

db = VaspCalcDb.from_db_file(db_json)


es = db.collection.aggregate(
    [
        {"$match": {"task_label":"HSE_scf"}},
        {"$group": {"_id":"$output.spacegroup.point_group",
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