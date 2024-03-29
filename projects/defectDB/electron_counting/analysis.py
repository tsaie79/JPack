from mpinterfaces.utils import *

import pandas as pd

from pymatgen import Structure

from qubitPack.qc_searching.analysis.main import get_defect_state
from qubitPack.tool_box import *

from matplotlib import pyplot as plt
import os


C2DB = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan", port=1234)
SCAN2dMat = get_db("Scan2dMat", "calc_data",  user="Jeng_ro", password="qimin", port=1236)
SCAN2dDefect = get_db("Scan2dDefect", "calc_data",  user="Jeng_ro", password="qimin", port=1236)
C2DB_IR = get_db("C2DB_IR", "calc_data",  user="Jeng_ro", password="qimin", port=1234)
#%%
class PureQuery:
    @classmethod
    def query(cls):
        db = SCAN2dDefect
        col = db.collection
        es = col.find({"task_id": {"$in": [24, 28]}})


    @classmethod
    def aggregate_defect_site_sym(cls):
        """
        refer to defect_site-sym_07292021.xlsx
        :return:
        """
        db = SCAN2dDefect
        filter = {
            "task_label":"SCAN_nscf uniform",
            "defect_name": {"$regex": "as"},
            "site_symmetry_uniform": True
        }

        es = db.collection.aggregate(
            [
                {"$match": filter},
                {"$group": {
                    "_id": {"pc_from_id": "$pc_from_id", "gp_id": "$group_id"},
                    "count": {"$sum":1},
                    "id": {"$push": "$task_id"},
                    "nsites": {"$push": "$nsites"},
                    "chemsys": {"$push": "$chemsys"},
                    "defect_name": {"$push": "$defect_name"},
                    "charge_state": {"$push": "$charge_state"},
                    "mag": {"$push":{"$round": [{"$arrayElemAt": ["$calcs_reversed.output.outcar.total_magnetization", 0]}, 3]}},
                    "pg": {"$push": "$output.spacegroup.point_group"},
                    "host_vbm_cbm_same": {"$push": "$host_info.scan_bs.band_edges.is_vbm_cbm_from_same_element"},
                    "scan_bg": {"$push": {"$round": ["$host_info.scan_bs.bandgap",3]}},
                    "is_gap_direct": {"$push": "$host_info.scan_bs.is_gap_direct"},
                    "cbm_up": {"$push": "$host_info.scan_bs.band_edges.cbm.up.max_element"},
                    "cbm_dn": {"$push": "$host_info.scan_bs.band_edges.cbm.down.max_element"},
                    "vbm_up": {"$push": "$host_info.scan_bs.band_edges.vbm.up.max_element"},
                    "vbm_dn": {"$push": "$host_info.scan_bs.band_edges.vbm.down.max_element"}
                }},
                {"$project":
                    {
                        "_id": 0,
                        "group_id": "$_id.gp_id",
                        "pc_from_id": "$_id.pc_from_id",
                        "chemsys": {"$arrayElemAt": ["$chemsys", 0]},
                        "task_id": "$id",
                        "pg": "$pg",
                        "defect_name": "$defect_name",
                        "charge_state": "$charge_state",
                        "mag": "$mag",
                        "scan_bg": {"$arrayElemAt": ["$scan_bg", 0]},
                        "vbm_el": [{"$arrayElemAt": ["$vbm_up", 0]}, {"$arrayElemAt": ["$vbm_dn", 0]}],
                        "cbm_el": [{"$arrayElemAt": ["$cbm_up", 0]}, {"$arrayElemAt": ["$cbm_dn", 0]}],
                        "vbm_eq_cbm": {"$arrayElemAt": ["$host_vbm_cbm_same", 0]},

                        "is_bg_direct": {"$arrayElemAt": ["$is_gap_direct", 0]},
                        "nsites": {"$arrayElemAt": ["$nsites", 0]},
                        "count": "$count",
                    }
                },
                {"$sort": {"group_id":1, "vbm_eq_cbm":-1, "scan_bg":-1}}
            ]
        )

        df = pd.DataFrame(es)
        return df

    @classmethod
    def aggregate_defect_triplet_and_sole_bandedges(cls):
        db = SCAN2dDefect
        filter = {
            "task_label":"SCAN_nscf uniform",
            "defect_name": {"$regex": "vac"},
            "mag": {"$nin": [2, -2]},
            # "vbm_eq_cbm": True
            # "scan_bg": {"$gte":0}
            # "site_symmetry_uniform": True
        }

        es = db.collection.aggregate(
            [
                {"$project": {
                    "_id": 0,
                    "group_id":1,
                    "task_label": 1,
                    "mag": {"$round": [{"$arrayElemAt": ["$calcs_reversed.output.outcar.total_magnetization", 0]}, 3]},
                    "group_id": 1,
                    "pc_from_id": 1,
                    "chemsys": 1,
                    "task_id": 1,
                    "pg": "$output.spacegroup.point_group",
                    "defect_name": 1,
                    "charge_state": 1,
                    "site_symmetry_uniform": 1,
                    "nsites": 1,

                    "scan_bg": {"$ifNull": [{"$round": ["$host_info.scan_bs.bandgap",3]}, None]},
                    "vbm_el": {"$ifNull": [["$host_info.scan_bs.band_edges.vbm.up.max_element","$host_info.scan_bs.band_edges.vbm.down.max_element"],None]},
                    "cbm_el": {"$ifNull": [["$host_info.scan_bs.band_edges.cbm.up.max_element","$host_info.scan_bs.band_edges.cbm.down.max_element"], None]},
                    "vbm_eq_cbm": {"$ifNull": ["$host_info.scan_bs.band_edges.is_vbm_cbm_from_same_element", None]},
                    "is_bg_direct": {"$ifNull": ["$host_info.scan_bs.is_gap_direct", None]},

                }},
                {
                    "$match": filter
                },
                {
                    "$group": {
                        "_id": {"pc_from_id": "$pc_from_id", "gp_id": "$group_id"},
                        "id": {"$push": "$task_id"},
                        "nsites": {"$push": "$nsites"},
                        "chemsys": {"$push": "$chemsys"},
                        "defect_name": {"$push": "$defect_name"},
                        "charge_state": {"$push": "$charge_state"},
                        "mag": {"$push":"$mag"},
                        "scan_bg": {"$push": "$scan_bg"},
                        "vbm_eq_cbm": {"$push": "$vbm_eq_cbm"},
                        "vbm_el": {"$push": "$vbm_el"},
                        "cbm_el": {"$push": "$cbm_el"},
                        "pg": {"$push": "$pg"},
                        "is_gap_direct": {"$push": "$is_bg_direct"},
                        "count": {"$sum":1},
                    }
                },
                {
                    "$project": {
                        "_id": 0,
                        'id': 1,
                        "gp_id": "$_id.gp_id",
                        "pc_from_id": "$_id.pc_from_id",
                        "nsites": 1,
                        "chemsys": {"$arrayElemAt": ["$chemsys", 0]},
                        "defect_name": 1,
                        "charge_state": 1,
                        "mag": 1,
                        "scan_bg": {"$arrayElemAt": ["$scan_bg", 0]},
                        "vbm_eq_cbm": {"$arrayElemAt": ["$vbm_eq_cbm", 0]},
                        "vbm_el": {"$arrayElemAt": ["$vbm_el", 0]},
                        "cbm_el": {"$arrayElemAt": ["$cbm_el", 0]},
                        "pg": 1,
                        "is_gap_direct": {"$arrayElemAt": ["$is_gap_direct", 0]},
                        "count": 1,
                    }
                },
                {"$sort": { "gp_id":1, "scan_bg":-1, "vbm_eq_cbm":-1}}

            ]
        )

        df = pd.DataFrame(es)
        return df








df = PureQuery.aggregate_defect_triplet_and_sole_bandedges()
df.to_clipboard()
#%%
class ExtractStructure:
    @classmethod
    def defect_structure(cls, task_id):
        os.chdir("/Users/jeng-yuantsai/anaconda3/envs/workflow/lib/python3.7/site-packages/JPack/projects/defectDB/electron_counting")
        db = SCAN2dDefect
        col = db.collection
        es = col.find({"task_id": {"$in": [task_id]}})

        for e in es:
            st = Structure.from_dict(e["output"]["structure"])
            os.makedirs("defect_structures/group_id-{}".format(e["group_id"]), exist_ok=True)
            st.to("POSCAR","defect_structures/group_id-{}/gp_{}-tk_{}.vasp".format(e["group_id"], e["group_id"], e["task_id"]))

ExtractStructure.defect_structure(973)
#%%
class ExtractDefectES:
    @classmethod
    def defect_levels(cls, task_id):
        defect_db = get_db("Scan2dDefect", "calc_data", port=1236)
        ir_db = get_db("Scan2dDefect", "ir_data", port=1236)
        # defect_db = get_db("antisiteQubit", "perturbed", port=1234)
        host_db = get_db("Scan2dMat", "calc_data", port=1236)


        tk_id = task_id
        defect = defect_db.collection.find_one({"task_id": tk_id})
        pc_from_id = defect["pc_from_id"]
        defect_name = defect["defect_name"]
        charge_state = defect["charge_state"]

        tot, proj, d_df = get_defect_state(
            defect_db,
            {"task_id": tk_id},
            5, -5,
            None,
            True,
            "proj",
            (host_db, pc_from_id, 0, 0, 0),
            0.001,
            is_vacuum_aligment_on_plot=True,
            locpot_c2db=None,#(c2db, c2db_uid, 0)
            ir_db=ir_db,
            ir_entry_filter={"pc_from_id": pc_from_id, "defect_name": defect_name, "charge_state": charge_state},
            top_texts=None
        )


        return tot, proj, d_df

for j in [3298]:
    tot, proj, d_df = ExtractDefectES.defect_levels(j)

#%%
class Potential:
    @classmethod
    def plot_pot(cls, tk_id):
        defect = SCAN2dDefect.collection.find_one({"task_id": tk_id})
        pc_from_id = defect["pc_from_id"]
        host = SCAN2dMat.collection.find_one({"task_id": pc_from_id})
        defect_name = defect["defect_name"]
        charge_state = defect["charge_state"]

        defect_pot = defect["calcs_reversed"][0]["output"]["locpot"]["2"]
        host_pot = host["calcs_reversed"][0]["output"]["locpot"]["2"]
        print("defct: {:.4}, host: {:.4}, delta: {:.4}".format(min(defect_pot), min(host_pot), min(defect_pot) - min(host_pot)))
        plt.plot(defect_pot)
        plt.plot(host_pot)
        plt.show()
        return plt1, plt2
plt1, plt2 = Potential.plot_pot(802)


#%%
def main():
    ExtractStructure.defect_structure()
    # dd = PureQuery.aggregate()
    # return dd
if __name__ == '__main__':
    main()
