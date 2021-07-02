from pymatgen import Structure
from atomate.vasp.database import VaspCalcDb
from pycdt.core.defectsmaker import ChargedDefectsStructures
import pandas as pd
from matplotlib import pyplot as plt




#%% host bandgap
from atomate.vasp.database import VaspCalcDb
from pymatgen import Element
import pandas as pd
from matplotlib import pyplot as plt
import os
from collections import Counter
import numpy as np
from qubitPack.qc_searching.analysis.read_eigen import DetermineDefectState

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'scan_relax_pc/db.json')


es = db.collection.aggregate(
    [
        {"$match": {"task_label":"SCAN_relax"}},
        {"$project": {
            "_id":0,
            "SCAN_gap": "$output.bandgap",
            "HSE_nosoc_gap": "$c2db_info.gap_hse_nosoc"}
        },
    ]
)

df = pd.DataFrame(es)
ax = df.plot.hist(bins=15, alpha=0.5)
plt.show()

#%% host space group
from atomate.vasp.database import VaspCalcDb
import pandas as pd
from matplotlib import pyplot as plt


db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'scan_relax_pc/db.json')


es = db.collection.aggregate(
    [
        {"$match": {"task_label":"SCAN_relax"}},
        {"$group":{"_id": "$output.spacegroup.symbol", "count": {"$sum":1}}},
        {"$project": {
            "_id":0,
            "spg": "$_id",
            "count": "$count"
            # "space_group": "$output.spacegroup.symbol"
        }},
        {"$sort": {"spg":1}}
    ]
)

df1 = pd.DataFrame(es)
df1.plot(x="spg", y="count", kind="bar", rot=20)
plt.show()



#%% host space group
from pymatgen import Structure
from atomate.vasp.database import VaspCalcDb
from pycdt.core.defectsmaker import ChargedDefectsStructures
import pandas as pd
from matplotlib import pyplot as plt


db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'search_triplet_from_defect_db/db.json')


es = db.collection.aggregate(
    [
        {"$match": {"task_label":"SCAN_scf",
                    "$or": [
                        {"calcs_reversed.output.outcar.total_magnetization": {"$gte":1.9, "$lte":2.1}},
                        {"calcs_reversed.output.outcar.total_magnetization": {"$gte":-2.1, "$lte":-1.9}}
                    ]}},
        {"$group": {"_id": {
                            "pc": "$pc_from",
                            # "pg": "$output.spacegroup.point_group",
                            # "q": "$charge_state",
                            # "defect": "$defect_entry.defect_type"
        },
                    "tid": {"$push":"$task_id"},
                    # "pc": {"$push": "$pc_from"},
                    "count": {"$sum":1},
                    }},
        {"$project": {
            "pc": "$_id.pc",
            "_id":0,
            # "pg": "$_id.pg",
            # "q": "$_id.q",
            # "d": "$_id.defect",
            "count": "$count",
            # "tid": "$tid"
        }},
        {"$sort": {"pc":1}}
    ]
)

df = pd.DataFrame(es)
df.to_clipboard("/t")
#%% site sym
from atomate.vasp.database import VaspCalcDb
import pandas as pd
from matplotlib import pyplot as plt


db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'scan_relax_pc/db.json')


es = db.collection.aggregate(
    [
        {"$match": {"task_label":"SCAN_relax"}},
        {"$group": {
            "_id": {"site_sym": "$sym_data.combo.syms", "spg": "$sym_data.international"},
            "count": {"$sum":1},
            "tid": {"$push": "$task_id"}
        }},
        {"$project": {"_id": 0,
                      "spg": "$_id.spg",
                      "site_sym": "$_id.site_sym",
                      "counts": "$count",
                      "tid": "$tid"
                      }},
        {"$sort": {"spg":-1}}
    ]
)

df = pd.DataFrame(es)
# df.to_json("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/site_sym_01012021.json", indent=4)
# df.to_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/site_sym_01012021.xlsx")

df.plot(x="site_sym", y="counts", kind="barh", rot=0, figsize=(10,8))

plt.show()
#%% Categorize defects by its site symmetry
from monty.serialization import loadfn, dumpfn

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'search_triplet_from_defect_db/db.json')

site_sym_cat = loadfn("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/site_sym_01012021.json")

site_sym_cat.update({"triplet_pc_tids": {}})
site_sym_cat.update({"triplet__pc_counts": {}})

pc_cat_tid = site_sym_cat["tid"]

for pc in range(12):
    pc = "{}".format(pc)

    es = db.collection.aggregate(
        [
            {"$match": {"pc_from_id": {"$in": pc_cat_tid[pc]},
                        "task_label": "SCAN_scf",
                        "$or": [
                            {"calcs_reversed.output.outcar.total_magnetization": {"$gte":1.9, "$lte":2.1}},
                            {"calcs_reversed.output.outcar.total_magnetization": {"$gte":-2.1, "$lte":-1.9}}
                        ]}},
            {"$group": {"_id": {
                "pc_id": "$pc_from_id",
                # "defect_name": "$defect_entry.name",
                # "pg": "$output.spacegroup.point_group",
                # "q": "$charge_state",
                # "tids": "$task_id",
                # "mag": "$calcs_reversed.output.outcar.total_magnetization"
            },
                "chemsys": {"$push": "$chemsys"},
                "defect_name": {"$push": "$defect_entry.name"},
                "pg": {"$push": "$output.spacegroup.point_group"},
                "q": {"$push": "$charge_state"},
                "class": {"$push": "$class"},
                "tids": {"$push": "$task_id"},
                "count": {"$sum":1},
            }},
            {"$project": {
                "_id":0,
                "pc_id": "$_id.pc_id",
                "chemsys": "$chemsys",
                "q": "$q",
                "d_n": "$defect_name",
                "tids": "$tids",
                "pg": "$pg",
                "class": "$class",
                # "mag": "$_id.mag",
                "count": "$count",
            }},
            {"$sort": {"pc_id":1}}
        ]
    )

    df = pd.DataFrame(es)
    entry = df.to_dict()

    if not entry == {}:
        triplet_tids = list(entry["pc_id"].values())
        site_sym_cat["triplet_pc_tids"].update({pc: triplet_tids})
        site_sym_cat["triplet_pc_counts"].update({pc: len(triplet_tids)})
        df.to_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/results/{}.xlsx".format(pc))
    else:
        site_sym_cat["triplet_pc_tids"].update({pc: []})
        site_sym_cat["triplet_pc_counts"].update({pc: 0})

pd.DataFrame(site_sym_cat).to_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/"
                                    "triplet_site_sym_01012021.xlsx")

#%% look into detail for each site-sym cat of defect
from monty.serialization import loadfn, dumpfn
import json

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'search_triplet_from_defect_db/db.json')

triplet_site_sym_cat = pd.read_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/"
                                     "xlsx/site_sym_cat/triplet_site_sym_01012021.xlsx", index_col=0)

triplet_cat_tids = triplet_site_sym_cat["triplet_pc_tids"][4]

es = db.collection.aggregate(
    [
        {"$match": {"pc_from_id": {"$in": json.loads(triplet_cat_tids)},
                    "task_label": "SCAN_scf",
                    "$or": [
                        {"calcs_reversed.output.outcar.total_magnetization": {"$gte":1.9, "$lte":2.1}},
                        {"calcs_reversed.output.outcar.total_magnetization": {"$gte":-2.1, "$lte":-1.9}}
                    ]}},
        {"$group": {
            "_id": {
            # "pc_id": "$pc_from_id",
            # "defect_name": "$defect_entry.name",
            "defect_type": "$defect_entry.defect_type",
            # "pg": "$output.spacegroup.point_group",
            "q": "$charge_state",
            # "tids": "$task_id",
            # "mag": "$calcs_reversed.output.outcar.total_magnetization"
        },
            "chemsys": {"$push": "$chemsys"},
            "pg": {"$push": "$output.spacegroup.point_group"},
            # "q": {"$push": "$charge_state"},
            "class": {"$push": "$class"},
            "tids": {"$push": "$task_id"},
            "defect_type": {"$push": "defect_entry.defect_type"},
            "defect": {"$push":"$defect_entry.name"},
            "count": {"$sum":1},
        }},
        {"$project": {
            "_id":0,
            "defect_type": "$_id.defect_type",
            "defect": "$defect",
            "chemsys": "$chemsys",
            "q": "$_id.q",
            "tids": "$tids",
            "pg": "$pg",
            "class": "$class",
            # "mag": "$_id.mag",
            "count": "$count",
        }},
        {"$sort": {"defect":1}}
    ]
)

df = pd.DataFrame(es)

#%%
d = pd.read_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/"
                  "xlsx/site_sym_cat/temp.xlsx")
d.fillna(0, inplace=True)

d["ant"] = d[2] / d[1]
d["vac"] = d[3] / d[1]
d["q1_ant"] = d[4] / d[1]
d["q0_ant"] = d[6] / d[1]
d["q-1_ant"] = d[8] / d[1]
d["q1_vac"] = d[5] / d[1]
d["q0_vac"] = d[7] / d[1]
d["q-1_vac"] = d[9] / d[1]

d.fillna(0, inplace=True)
d = d[["q1_ant", "q0_ant", "q-1_ant", "q1_vac", "q0_vac", "q-1_vac"]]

d.to_clipboard(excel=True)

#%% plausible qubit
from qubitPack.qc_searching.analysis.main import get_defect_state
import os

proj_path = '/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/search_triplet_from_defect_db'
db_json = os.path.join(proj_path, "db.json")

host_path = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc"
db_host_json = os.path.join(host_path, "db.json")
for i in [1,2,3,4,5,7,8,9]:
    if i != 3:
        continue
    df = pd.read_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/results/{}.xlsx".format(i)
                       , index_col=0)
    data = {}
    data.update({"qubit_check": []})
    data.update({"error": []})
    print(df["tids"])
    for n, tids in enumerate(df["tids"]):
        ers = []
        checks = []
        print(n)
        from json import loads
        tids = loads(tids)
        print(tids)
        for tid in tids:
            try:
                tot, proj, d_df = get_defect_state(
                    db_json,
                    {"task_id": tid},
                    2,-2,
                    None,
                    None,
                    "dist",
                    db_host_json,
                    0.2
                )
                ers.append(None)
                checks.append(True)

            except Exception as er:
                print(er)
                ers.append(er)
                checks.append(False)

        print(data)
        data["error"].append(ers)
        data["qubit_check"].append(checks)
        add_info = pd.DataFrame(data)
    print(add_info)
    result = pd.concat([df, add_info], axis=1)
    result.to_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/results/{}.xlsx".format(i))
#%%
a = []
for i in [1,2,3,4,5,7,8,9]:
    a.append(pd.read_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/results/{}.xlsx".format(i), index_col=0))
a = pd.concat(a, axis=0, ignore_index=True)
a.to_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/results/combine.xlsx", index=False, sheet_name="combine")
#%% Pick pc that only has good output of read_defect
df = pd.read_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/results/combine.xlsx", sheet_name="combine")

a = df.to_dict()
a.update({"good": {}})
for col, checks in a["qubit_check"].items():
    if "True" in checks:
        a["good"].update({col: True})
    else:
        a["good"].update({col: False})
a = pd.DataFrame(a)
good = a[a.good==True]
good.drop("good", axis=1,inplace=True)


bad = a[a.good==False]
bad.drop("good", axis=1,inplace=True)

with pd.ExcelWriter("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/results/temp.xlsx") as f:
    good.to_excel(f, index=False, sheet_name="good")
    bad.to_excel(f, index=False, sheet_name="bad")

#%% Pick plausible qubit
# 1. electronic structure of defect center
# 1. defect state from vbm > 0.5 eV
# 2. at least one unoccupied state in at least one of spin channels

"""
triplet_from: dn
    dn [0, 1, 1] up [0, 0, 0]
    dn [1, 1] up []
     
    
"""


from qubitPack.qc_searching.analysis.main import get_defect_state
import os

proj_path = '/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/search_triplet_from_defect_db'
db_json = os.path.join(proj_path, "db.json")

host_path = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc"
db_host_json = os.path.join(host_path, "db.json")

df = pd.read_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/results/combine.xlsx"
                   ,index_col=0, sheet_name="good")
data = {}
data.update({"qubit_check": []})
data.update({"error": []})
print(df["tids"])
for n, tids in enumerate(df["tids"]):
    ers = []
    checks = []
    print(n)
    from json import loads
    tids = loads(tids)
    print(tids)
    for tid in tids:
        try:
            tot, proj, d_df = get_defect_state(
                db_json,
                {"task_id": tid},
                2,-2,
                None,
                None,
                "dist",
                db_host_json,
                0.1
            )

            d_df = d_df.to_dict()[0]


            ers.append(None)
            checks.append(True)

        except Exception as er:
            print(er)
            ers.append(er)
            checks.append(False)

    print(data)
    data["error"].append(ers)
    data["qubit_check"].append(checks)
    add_info = pd.DataFrame(data)
print(add_info)
result = pd.concat([df, add_info], axis=1)
result.to_excel("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/xlsx/site_sym_cat/results/{}.xlsx".format(i))