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
                      "compounds": "$count",
                      "tid": "$tid"
                      }},
        {"$sort": {"spg":-1}}
    ]
)

df = pd.DataFrame(es)

df.plot(x="spg", y="compounds", kind="barh", rot=0, figsize=(10,8))

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
