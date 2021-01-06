from pymatgen import Structure
from atomate.vasp.database import VaspCalcDb
from pycdt.core.defectsmaker import ChargedDefectsStructures
import pandas as pd
from matplotlib import pyplot as plt


#%%
db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/code/JPack/qubitPack/database_profile/db_qimin_db.json")
# db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
#                              'scan_relax_pc/db.json')
# db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/code/JPack/qubitPack/database_profile/db_dk_local.json")

es = db.collection.aggregate(
    [
        {"$match": {
            "functional_type": "SCAN",
            "task_type":"static",
            "calcs_reversed.output.outcar.total_magnetization": {"$gt":0}
        }},
        {"$group": {"_id": {"gap": "$output.bandgap"}}},
        {"$project": {
            "_id":0,
            "c2db_PBE_gap": "$_id.gap",
        }}
    ]
)
# print(len(list(es)))

df = pd.DataFrame(es)
ax = df.plot.hist(bins=50, alpha=1)
plt.show()

#%%
db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/code/JPack/qubitPack/database_profile/db_qimin_db.json")

es = db.collection.aggregate(
    [   {"$unwind": "$calcs_reversed"},
        {"$match": {"task_type":"static", "functional_type": "SCAN", "output.bandgap":{"$gte":0.1},
                    "calcs_reversed.output.outcar.total_magnetization": {"$gte":0.1}
                    }},
        {"$group": {"_id": {"mag": "$calcs_reversed.output.outcar.total_magnetization"},
                    "tid": {"$push":"$task_id"},
                    "chemsys": {"$push": "$chemsys"},
                    "formula": {"$push": "$formula_pretty"},
                    "cbm": {"$push": "$output.cbm"},
                    "vbm": {"$push": "$output.vbm"},
                    "bandgap": {"$push": "$output.bandgap"},
                    "locpot": {"$push": "$vacuum_locpot.value"},
                    "space_group": {"$push": "$output.spacegroup.symbol"},
                    "a": {"$push": "$input.structure.lattice.a"},
                    "b": {"$push": "$input.structure.lattice.b"},
                    "c": {"$push": "$input.structure.lattice.c"},
                    }},
        {"$project": {"formula": "$formula","chemsys": "$chemsys", "mag": "$_id.mag","bandgap": "$bandgap", "cbm": "$cbm", "vbm": "$vbm",
                      "locpot": "$locpot", "space_group": "$space_group", "a": "$a", "b":"$b", "c": "$c", "tids": "$tid", "_id":0}},

        {"$sort": {"bandgap":1, "mag":1}}
    ]
)


df = pd.DataFrame(es)
d = df.to_json(indent=4)
from json import loads
d = loads(d)

for k, v in d.items():
    print(k, v)
    if k != "mag":
        d[k].update(dict(zip(v.keys(), [i[0] for i in v.values()])))
    else:
        continue
df = pd.DataFrame(d)
df["workfunction"] = df["locpot"]-df["vbm"]

df.to_clipboard(excel=True)
# ax = df.plot.hist(bins=15, alpha=1)
# plt.show()

#%%
db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/code/JPack/qubitPack/database_profile/db_dk_local.json")


es = db.collection.aggregate(
    [   {"$match": {"gap_hse_nosoc":{"$gte":0.1},
                    "magstate": {"$in": ["FM", "AFM"]},
                    }},
        {"$project": {"_id":0, "formula": "$formula", "uid": "$uid", "gap_hse_nosoc": "$gap_hse_nosoc",
                      "vbm_hse_nosoc": "$vbm_hse_nosoc",
                      "cbm_hse_nosoc": "$cbm_hse_nosoc",
                      "evacmean": "$evacmean",
                      "evac_diff": "$evacdiff",
                      "spacegroup": "$spacegroup",
                      "a": "$structure.lattice.a",
                      "b": "$structure.lattice.b",
                      "c": "$structure.lattice.c"
                      }
         },

        {"$sort": {"gap_hse_nosoc":1}}
        ]
)


df = pd.DataFrame(es)
df["workfunction"] = df["evacmean"] - df["vbm_hse_nosoc"]
# d = df.to_json(indent=4)
# from json import loads
# d = loads(d)
#
# for k, v in d.items():
#     print(k, v)
#     if k != "mag":
#         d[k].update(dict(zip(v.keys(), [i[0] for i in v.values()])))
#     else:
#         continue
# df = pd.DataFrame(d)
df.to_clipboard(excel=True)
ax = df.plot.hist(bins=15, alpha=1)
plt.show()