#%% 1
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

defect_db = get_db("Scan2dDefect", "calc_data", port=1236)
host_db = get_db("Scan2dMat", "calc_data", port=1236)


tk_id = 37
pc_from_id = defect_db.collection.find_one({"task_id": tk_id})["pc_from_id"]
# c2db_uid = host_db.collection.find_one({"task_id": pc_from_id})["c2db_info"]["uid"]

tot, proj, d_df = get_defect_state(
    defect_db,
    {"task_id": tk_id},
    1, -5,
    None,
    True,
    "proj",
    (host_db, pc_from_id, 0), #(get_db("antisiteQubit", "W_S_Ef"), 312, 0.),
    0.1,
    locpot_c2db=None #(c2db, c2db_uid, 0)
)

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
            "_id": {"site_sym": "$sym_data.high_dim_ir_info.syms", "spg": "$sym_data.international"},
            "count": {"$sum":1},
            "tid": {"$push": "$task_id"},
            "uid": {"$push": "$c2db_info.uid"},
        }},
        {"$project": {"_id": 0,
                      "spg": "$_id.spg",
                      "site_sym": "$_id.site_sym",
                      "counts": "$count",
                      "uid": "$uid",
                      "tid": "$tid",
                      }},
        {"$sort": {"site_sym": -1, "spg":1}}
    ]
)

df = pd.DataFrame(es)
df.to_json("/Users/jeng-yuantsai/Research/project/defectDB/symBaseBinaryQubit/xlsx/site_sym_cat/site_sym_07062021.json", indent=4)
df.to_clipboard()

df.plot(x="site_sym", y="counts", kind="barh", rot=0, figsize=(10,8))

#%% 4 update entries by including info of site sym
from qubitPack.tool_box import get_db, get_good_ir_sites
from monty.json import jsanitize
from pymatgen.io.vasp.inputs import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

db = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=1236)

for e in db.collection.find():
    st = Structure.from_dict(e["output"]["structure"])
    sym = SpacegroupAnalyzer(st, symprec=1e-1).get_symmetry_dataset()
    d = jsanitize(sym)
    print(d)
    d.update({"symprec":1e-1})
    species = [str(sp) for sp in st.species]
    combine = dict(zip(d["wyckoffs"], zip(species, d["site_symmetry_symbols"])))
    species, syms = get_good_ir_sites(dict(combine.values()).keys(), dict(combine.values()).values())
    combine.update({"species": species, "syms": syms})

    d.update({"high_dim_ir_info": combine})
    print(d)
    db.collection.update_one({"task_id": e["task_id"]}, {"$set":{"sym_data":d}})

#%% 5 generate "cat_host_bg>1_and_homo" json file

df = pd.read_excel(
    "/Users/jeng-yuantsai/Research/project/defectDB/symBaseBinaryQubit/xlsx/site_sym_cat/host_scan_nscf_line_07062021.xlsx",
    sheet_name="cat_host_bg>1_and_inhomo"
)
df.to_json("/Users/jeng-yuantsai/Research/project/defectDB/symBaseBinaryQubit/xlsx/site_sym_cat/bg_lge_1_and_inhomo_07062021", indent=4, orient="records")

#%% 6
from atomate.vasp.database import VaspCalcDb
import pandas as pd
from matplotlib import pyplot as plt
from qubitPack.tool_box import get_db


db = get_db("Scan2dMat", "calc_data", port=1236)


es = db.collection.aggregate(
    [
        {"$match": {"task_label":"SCAN_nscf line"}},
        {"$group": {
            "_id": {"nsites": "$nsites"},
            "count": {"$sum":1},
            "AB": {"$push": "$formula_anonymous"}
        }},
        {"$project": {"_id": 0,
                      "nsites": "$_id.nsites",
                      "counts": "$count",
                      "AB": "$AB"
                      }},
        {"$sort": {"nsites": -1}}
    ]
)

df = pd.DataFrame(es)
df.to_clipboard()