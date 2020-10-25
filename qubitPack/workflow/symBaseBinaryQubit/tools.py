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

#%% search for triplet
from atomate.vasp.database import VaspCalcDb

from pymatgen import Element

import pandas as pd

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'search_triplet_from_defect_db/db.json')

es = db.collection.aggregate(
    [
        {"$match": {"task_label":"SCAN_scf"}},
        {"$addFields":{"mag": {"$abs": {"$arrayElemAt": ["$calcs_reversed.output.outcar.total_magnetization",0]}}}},
        {"$addFields":{"defect_name":"$defect_entry.name", "point_group":"$output.spacegroup.point_group"}},
        {"$project": {"_id":0, "formula_pretty":1, "nsites":1, "elements":1,
                      "mag":1, "charge_state":1, "chemsys":1, "defect_name":1,
                      "point_group":1
                      }},
        {"$match":{"mag":{"$gte":1.9, "$lte":2.1}}},
        {"$sort":{"mag":-1}},
        # {"$group": {"_id":"$chemsys", "q": {"$push":"$charge_state"},"m": {"$push":"$formula_pretty"}, "nsites":{"$push":"$nsites"}}}
    ]
)

ess = []
for e in es:
    els = [Element(e["elements"][0]).group, Element(e["elements"][1]).group]
    els.sort()
    e.update({"group": "{}-{}".format(els[0], els[1])})
    # e.pop("elements")
    ess.append(e)
df = pd.DataFrame(ess).sort_values(["group","defect_name", "mag"], inplace=False).set_index("formula_pretty")

#%% Build the alignment of cbm and vbm of defect and host materials
from atomate.vasp.database import VaspCalcDb

import numpy as np

from matplotlib import pyplot as plt


defect = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                                 'search_triplet_from_defect_db/db.json')


host = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                               'scan_relax_pc/db.json')

test_host = host.collection.find_one({"task_id": 1})
vac_host = max(test_host["calcs_reversed"][2]["output"]["locpot"]["2"])

vac_vbm = test_host["output"]["vbm"]
vac_cbm = test_host["output"]["cbm"]

print(vac_host-vac_cbm, vac_host-vac_vbm)


