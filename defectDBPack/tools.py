#%% search_triplet_from_defect_db
from atomate.vasp.database import VaspCalcDb

from pymatgen import Element

import pandas as pd

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'search_triplet_from_defect_db/db.json')

pc_db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc/db.json")

d = []

filter = {
    "$or":[{"calcs_reversed.output.outcar.total_magnetization":{"$gte":1.5}},
           {"calcs_reversed.output.outcar.total_magnetization":{"$lte":-1.5}}],
    "nsites":75,
    "task_label":"SCAN_scf"
}

for e in db.collection.find(filter):
    entry = {}
    entry["task_id"] = e["task_id"]
    entry["formula_pretty"] = e["formula_pretty"]
    entry["mag"] = abs(e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"])
    entry["sym_data"] = pc_db.collection.find_one({"task_id":int(e["pc_from"].split("/")[-1])})["sym_data"]["site_symmetry_symbols"]
    entry["species"] = pc_db.collection.find_one({"task_id":int(e["pc_from"].split("/")[-1])})["sym_data"]["species"]
    species = list(dict.fromkeys(entry["species"]))
    entry["elec_configs"] = (Element(species[0]).electronic_structure, Element(species[1]).electronic_structure)
    entry["group"] = (Element(species[0]).group, Element(species[1]).group)
    entry["NN"] = e["NN"]
    entry["bandgap"] = pc_db.collection.find_one({"task_id":int(e["pc_from"].split("/")[-1])})["output"]["bandgap"]
    entry["pc_from"] = e["pc_from"]
    entry["a"] = e["input"]["structure"]["lattice"]["a"]
    entry["b"] = e["input"]["structure"]["lattice"]["b"]



    d.append(entry)

df = pd.DataFrame(d).sort_values(["sym_data"])# "group", "species"]) #"mag", "pc_from", "bandgap"])
df = df.set_index("task_id")
df.to_clipboard()

#%%