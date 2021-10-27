#%% 1
import json
import os.path


class Pc:
    @classmethod
    def hist_nsites(cls):
        import pandas as pd
        from matplotlib import pyplot as plt
        from qubitPack.tool_box import get_db

        db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan")
        filter = {
            "magstate": "NM",
            "gap_hse_nosoc": {"$gt":0},
            "nkinds": 3
        }

        es = db.collection.aggregate(
            [
                {"$match": filter},
                {"$group": {
                    "_id": {"gap": "$gap_hse_nosoc"},
                    "count": {"$sum":1},
                    "uid": {"$push": "$uid"}
                }},
                {"$project":
                    {
                    "_id": 0,
                    "gap": "$_id.gap",
                    "count": "$count",
                    }
                },
                {"$sort": {"gap":1}}
            ]
        )

        df = pd.DataFrame(es)
        print(df)

        df.plot(x="gap", y="count", kind="bar", rot=0, figsize=(10,8), title="Hist nsites under filter")
        plt.ylabel("counts")

        plt.show()

#%% 2 read defect states
from qubitPack.qc_searching.analysis.main import get_defect_state
from qubitPack.tool_box import get_db

from pymatgen import Structure
import os

defect_db = get_db("Scan2dDefect", "calc_data", port=1236)
# defect_db = get_db("antisiteQubit", "perturbed", port=1234)
host_db = get_db("Scan2dMat", "calc_data", port=1236)


tk_id = 2153
pc_from_id = defect_db.collection.find_one({"task_id": tk_id})["pc_from_id"]
# c2db_uid = host_db.collection.find_one({"task_id": pc_from_id})["c2db_info"]["uid"]

tot, proj, d_df = get_defect_state(
    defect_db,
    {"task_id": tk_id},
    1, -5,
    None,
    True,
    "proj",
    (host_db, pc_from_id, 0.5), #(get_db("antisiteQubit", "W_S_Ef"), 312, 0.),
    0.01,
    locpot_c2db=None #(c2db, c2db_uid, 0)
)

# tgt_path = "/Users/jeng-yuantsai/Research/project/Scan2dDefect/calculations/calc_data/defect_st"
# e = defect_db.collection.find_one({"task_id": tk_id})
# st = Structure.from_dict(e["output"]["structure"])
# st.to("poscar", os.path.join(tgt_path, "{}-{}.vasp".format(e["host_info"]["c2db_info"]["uid"], tk_id)))

#%% 3 site sym
from atomate.vasp.database import VaspCalcDb
import pandas as pd
from matplotlib import pyplot as plt
from qubitPack.tool_box import get_db


db = get_db("Scan2dMat", "calc_data", port=1236)


es = db.collection.aggregate(
    [
        {"$match": {"task_label":"SCAN_nscf line", "output.bandgap": {"$gte": 1}}},
        {"$group": {
            "_id": {"site_sym": "$sym_data.good_ir_info.site_sym",
                    "pg": "$sym_data.pmg_pg",
                    "spg": "$sym_data.pmg_spg",
                    "spg_number": "$sym_data.pmg_spg_number"
                    },
            "count": {"$sum":1},
            "tid": {"$push": "$task_id"},
            "uid": {"$push": "$c2db_info.uid"},
        }},
        {"$project": {"_id": 0,
                      "spg": "$_id.spg",
                      "spg_number": "$_id.spg_number",
                      "pg": "$_id.pg",
                      "site_sym": "$_id.site_sym",
                      "counts": "$count",
                      "uid": "$uid",
                      "tid": "$tid",
                      }},
        {"$sort": {"site_sym": -1, "spg":1}}
    ]
)

df = pd.DataFrame(es)
# df.to_json("/Users/jeng-yuantsai/Research/project/Scan2dDefect/xlsx/site_sym_cat/site_sym_07272021.json", indent=4)
df.to_clipboard()

df.plot(x="site_sym", y="counts", kind="barh", rot=0, figsize=(10,8))

#%% 4 update entries by including info of site sym
#%% 5 generate "cat_host_bg>1_and_homo" json file
import pandas as pd
df = pd.read_excel(
    "/Users/jeng-yuantsai/Research/project/Scan2dDefect/xlsx/site_sym_cat/site_sym_07272021.xlsx",
    sheet_name="cat_host_bg>=1_and_homo"
)
df.to_json("/Users/jeng-yuantsai/Research/project/Scan2dDefect/xlsx/site_sym_cat/bg_gte_1_and_homo_07272021", indent=4, orient="records")

#%% 6
from atomate.vasp.database import VaspCalcDb
import pandas as pd
from matplotlib import pyplot as plt
from qubitPack.tool_box import get_db
from monty.serialization import loadfn
from pymatgen import Structure
from pymatgen.electronic_structure.plotter import BSDOSPlotter
import json, os

db = get_db("Scan2dMat", "calc_data", port=1236)
tgt_path = "/Users/jeng-yuantsai/Research/project/Scan2dDefect/calculations/calc_data/host_st"

groups = loadfn("/Users/jeng-yuantsai/Research/project/Scan2dDefect/xlsx/site_sym_cat/bg_lge_1_and_homo_07062021")

for group in groups:
    tids = group["tid"]
    for tid in json.loads(tids)[:1]:
        e = db.collection.find_one({"task_id": tid})
        uid = e["c2db_info"]["uid"]
        st = Structure.from_dict(e["output"]["structure"])
        st.to("poscar", os.path.join(tgt_path, "{}-{}.vasp".format(uid, tid)))

        bs = db.get_band_structure(tid)
        dos = db.get_dos(db.collection.find_one({"task_label": "SCAN_nscf uniform", "c2db_info.uid": uid})["task_id"])
        plot = BSDOSPlotter()
        plot.get_plot(bs, dos)
        plt.suptitle("Homo-{}-{}-{}".format(group["group_id"], uid, tid))
        plt.show()

