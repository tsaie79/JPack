import numpy as np
from mpinterfaces.utils import *

import pandas as pd

from pymatgen import Structure

from qubitPack.qc_searching.analysis.main import get_defect_state_v1, get_defect_state_v2, get_defect_state_v3
from qubitPack.tool_box import *

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import os

C2DB = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan", port=12345)
SCAN2dMat = get_db("Scan2dMat", "calc_data",  user="Jeng_ro", password="qimin", port=12347)
SCAN2dDefect = get_db("Scan2dDefect", "calc_data",  user="Jeng_ro", password="qimin", port=12347)
C2DB_IR = get_db("C2DB_IR", "calc_data",  user="Jeng_ro", password="qimin", port=12345)

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







#
# df = PureQuery.aggregate_defect_triplet_and_sole_bandedges()
# df.to_clipboard()
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

# ExtractStructure.defect_structure(973)
class ExtractDefectES:
    @classmethod
    def defect_levels(cls, task_id):
        defect_db = get_db("Scan2dDefect", "calc_data", port=12347)
        ir_db = get_db("Scan2dDefect", "ir_data", port=12347)
        host_db = get_db("Scan2dMat", "calc_data", port=12347)

        tk_id = task_id
        defect = defect_db.collection.find_one({"task_id": tk_id})
        pc_from_id = defect["pc_from_id"]
        defect_name = defect["defect_name"]
        charge_state = defect["charge_state"]

        tot, proj, d_df, levels, in_gap_levels = get_defect_state_v3(
            defect_db,
            {"task_id": tk_id},
            -10, 10,
            None,
            "all",
            None,
            None, #(host_db, pc_from_id, 0, 0, 0),
            0.2,
            edge_tol=(.25, .25),
            is_vacuum_aligment_on_plot=True,
            locpot_c2db=None, #(c2db, c2db_uid, 0)
            ir_db=ir_db,
            ir_entry_filter={"pc_from_id": pc_from_id, "defect_name": defect_name, "charge_state": charge_state},
        )
        return tot, proj, d_df, levels, in_gap_levels

class Potential:
    def __init__(self, tk_id, defect_db, host_db):
        self.defect = defect_db.collection.find_one({"task_id": tk_id})
        self.pc_from_id = self.defect["pc_from_id"]
        self.host = host_db.collection.find_one({"task_id": self.pc_from_id})
        self.defect_pot = self.defect["calcs_reversed"][0]["output"]["locpot"]["2"]
        self.defect_pot_max, self.d_pot_max_index = max(self.defect_pot), self.defect_pot.index(max(self.defect_pot))

        self.host_pot = self.host["calcs_reversed"][0]["output"]["locpot"]["2"]
        self.host_pot_max, self.h_pot_max_index = max(self.host_pot), self.host_pot.index(max(self.host_pot))

    def model_f(self, x, a, b):
        return a*x+b

    def plot_pot(self):
        fig = plt.figure(figsize=(10,8), dpi=300)
        ax = fig.add_subplot(1, 3, 1)
        defect_name = self.defect["defect_name"]
        charge_state = self.defect["charge_state"]
        print("defct: {:.4}, host: {:.4}, delta: {:.4}".format(min(self.defect_pot), min(self.host_pot),
        min(self.defect_pot) - min(self.host_pot)))

        ax.plot(self.defect_pot, color="r", label="{}-{}-q{}".format(self.defect["task_id"], defect_name,
                                                                     charge_state))
        ax.plot(self.host_pot, color="blue", label="{}-{}".format(self.host["task_id"], self.host["c2db_info"]["uid"]))
        ax.text(self.d_pot_max_index, self.defect_pot_max, round(self.defect_pot_max, 3))
        ax.text(self.h_pot_max_index+15, self.host_pot_max, round(self.host_pot_max, 3))
        ax.legend(loc="center right", fancybox=False, edgecolor="black")

        #host electric field fit
        host_popt, defect_popt = self.linear_fit_potential()
        ax = fig.add_subplot(1, 3, 2)
        ax.plot(self.host_pot, color="blue", label="{}-{}".format(self.host["task_id"], self.host["c2db_info"]["uid"]))
        a_opt, b_opt = host_popt
        x_model = range(len(self.host_pot[:50]))
        y_model = self.model_f(x_model, a_opt, b_opt)
        ax.plot(x_model, y_model, color="orange", label="{:.1e}x+{:.2f}".format(a_opt, b_opt))
        ax.legend(loc="center right", fancybox=False, edgecolor="black")

        # defect electric field fig
        ax = fig.add_subplot(1, 3, 3)
        ax.plot(self.defect_pot, color="r", label="{}-{}-q{}".format(self.defect["task_id"], defect_name, charge_state))
        a_opt, b_opt = defect_popt
        x_model = range(len(self.defect_pot[:50]))
        y_model = self.model_f(x_model, a_opt, b_opt)
        ax.plot(x_model, y_model, color="orange", label="{:.1e}x+{:.2f}".format(a_opt, b_opt))
        ax.legend(loc="center right", fancybox=False, edgecolor="black")

        fig.show()

    def linear_fit_potential(self):
        host_popt, host_pcov = curve_fit(self.model_f, range(len(self.host_pot[:50])), self.host_pot[:50], p0=[0.1,
                                                                                                     self.host_pot_max])
        defect_popt, defect_pcov =\
            curve_fit(self.model_f, range(len(self.defect_pot[:50])), self.defect_pot[:50], p0=[0.1,
                                                                                                self.defect_pot_max])
        return host_popt, defect_popt

def main():
    for j in [3284]:
        tot, proj, d_df, levels, in_gap_levels = ExtractDefectES.defect_levels(j)    # dd = PureQuery.aggregate()
        return tot, proj, d_df, levels, in_gap_levels

if __name__ == '__main__':
    tot, proj, d_df, levels, in_gap_levels = main()
