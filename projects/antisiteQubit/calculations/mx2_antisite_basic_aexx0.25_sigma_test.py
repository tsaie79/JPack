import os
from qubitPack.qc_searching.analysis.main import get_defect_state

proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_sigma_test"
db_json = os.path.join(proj_path, "db.json")

host_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_pc"
db_host_json = os.path.join(host_path, "db.json")

perturb = None
chem = "Mo-S"
tot, proj, d_df = get_defect_state(
    db_json,
    {"task_label":"HSE_scf", "chemsys":chem, "perturbed":perturb},
    5,-5,
    os.path.join(proj_path),
    True,
    "dist",
    db_host_json,
    0.05
)
tot.to_clipboard()

#%% mx2_antisite_aexx0.25_sigma_test vs perturbed
"""
compare energy
"""
import os
import pandas as pd
from atomate.vasp.database import VaspCalcDb

p1 = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_sigma_test"
db1 = VaspCalcDb.from_db_file(os.path.join(p1, "db.json"))

p2 = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/perturbed"
db2 = VaspCalcDb.from_db_file(os.path.join(p2, "db.json"))

p3 = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_final"
db3 = VaspCalcDb.from_db_file(os.path.join(p3, "db.json"))

e1 = db1.collection.aggregate(
    [
        {"$match": {"task_label":"HSE_scf"}},
        {"$project":{
            "_id":0,
            "energy": "$output.energy",
            "chem": "$chemsys",
            "tid": "$task_id"
        },
        }
    ]
)
e2 = db2.collection.aggregate(
    [
        {"$match": {"task_label":"HSE_scf"}},
        {"$project":{
            "_id":0,
            "energy_perturb": "$output.energy",
            "chem": "$chemsys",
        },
        }
    ]
)
e3 = db3.collection.aggregate(
    [
        {"$match": {"task_label":"HSE_scf", "nupdown_set":2}},
        {"$project":{
            "_id":0,
            "energy_no_sigma": "$output.energy",
            "chem": "$chemsys",
        },
        }
    ]
)


e1 = pd.DataFrame(list(e1)).set_index("chem")
e2 = pd.DataFrame(list(e2)).set_index("chem")
e3 = pd.DataFrame(list(e3)).set_index("chem")

e = pd.concat([e1,e2], axis=1).round(3)
e["dE"] = e["energy"] - e["energy_perturb"]
# e1 = db1.find_one({"task_id": 4256})["output"]["energy"]
# e2 = db2.find_one({"task_id": 481})["output"]["energy"]
"""
MoS2 has little splitting on x2-y2 and xy, and z2 is occupied. I don't think this is a ground state based on 
large energetic difference between C3v and Ch.
"""