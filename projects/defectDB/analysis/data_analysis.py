import os, datetime

# from analysis.deep_defect import Potential

from qubitPack.tool_box import get_db, get_band_edges_characters, get_unique_sites_from_wy, get_good_ir_sites
from qubitPack.tool_box import IOTools

import pandas as pd
import numpy as np
# pd.set_option('display.max_rows', 500)

from pymatgen import Spin, Element
from pymatgen.io.vasp.inputs import Structure

from monty.serialization import loadfn

from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator

SCAN2dMat = get_db("Scan2dMat", "calc_data",  user="Jeng_ro", password="qimin", port=12347)
SCAN2dDefect = get_db("Scan2dDefect", "calc_data",  user="Jeng_ro", password="qimin", port=12347)
SCAN2dIR = get_db("Scan2dDefect", "ir_data", port=12347)

HSEQubitDefect = get_db("HSE_triplets_from_Scan2dDefect", "calc_data-pbe_pc", port=12347)
HSEQubitIR = get_db("HSE_triplets_from_Scan2dDefect", "ir_data-pbe_pc", port=12347)
HSEQubitCDFT = get_db("HSE_triplets_from_Scan2dDefect", "cdft-pbe_pc", port=12347, user="Jeng")

input_path = "analysis/input"
save_xlsx_path = "analysis/output/xlsx"
save_plt_path = "analysis/output/plt"


class Tools:
    @classmethod
    def extract_defect_levels(cls, defect_taskid):
        from qubitPack.qc_searching.analysis.main import get_defect_state_v3
        from qubitPack.tool_box import get_db

        from pymatgen import Structure
        import os


        defect_db =SCAN2dDefect
        host_db = SCAN2dMat
        ir_db = SCAN2dIR

        defect = defect_db.collection.find_one({"task_id": defect_taskid})
        pc_from_id = defect["pc_from_id"]
        defect_name = defect["defect_name"]
        charge_state = defect["charge_state"]
        level_info, levels = None, None
        try:
            state = get_defect_state_v3(
                defect_db,
                {"task_id": defect_taskid},
                -10, 10,
                None,
                False,
                None,
                None,  #(host_db, host_taskid, 0, vbm_dx, cbm_dx),
                0.2,  #0.2
                locpot_c2db=None,  #(c2db, c2db_uid, 0)
                is_vacuum_aligment_on_plot=True,
                ir_db=ir_db,
                ir_entry_filter={"pc_from_id": pc_from_id, "defect_name": defect_name, "charge_state": charge_state},
            )
            tot, proj, d_df, levels = state
            level_info = d_df.T.to_dict("records")[0]
        except Exception as er:
            print(er)
            level_info = {}
            levels = {}
        return level_info, levels

    @classmethod
    def extract_defect_levels_v2(cls, defect_taskid):
        from qubitPack.qc_searching.analysis.main import get_defect_state_v3
        from qubitPack.tool_box import get_db

        from pymatgen import Structure
        import os


        defect_db =SCAN2dDefect
        host_db = SCAN2dMat
        ir_db = SCAN2dIR

        defect = defect_db.collection.find_one({"task_id": defect_taskid})
        pc_from_id = defect["pc_from_id"]
        defect_name = defect["defect_name"]
        charge_state = defect["charge_state"]
        level_info, levels, defect_levels = None, None, None
        try:
            state = get_defect_state_v3(
                defect_db,
                {"task_id": defect_taskid},
                -10, 10,
                None,
                None,
                None,
                None,  #(host_db, host_taskid, 0, vbm_dx, cbm_dx),
                0.2,  #0.2
                locpot_c2db=None,  #(c2db, c2db_uid, 0)
                is_vacuum_aligment_on_plot=True,
                edge_tol=(0.25, 0.25),
                ir_db=ir_db,
                ir_entry_filter={"pc_from_id": pc_from_id, "defect_name": defect_name, "charge_state": charge_state},
            )
            tot, proj, d_df, levels, defect_levels = state
            level_info = d_df.to_dict("records")[0]

        except Exception as er:
            print(er)
            level_info = {}
            levels = {}
            defect_levels = {}
        return level_info, levels, defect_levels

    @classmethod
    def extract_defect_levels_v2_hse(cls, defect_taskid):
        from qubitPack.qc_searching.analysis.main import get_defect_state_v3
        from qubitPack.tool_box import get_db

        from pymatgen import Structure
        import os


        defect_db =HSEQubitDefect
        ir_db = HSEQubitIR

        defect = defect_db.collection.find_one({"task_id": defect_taskid})
        level_info, levels, defect_levels = None, None, None
        try:
            state = get_defect_state_v3(
                defect_db,
                {"task_id": defect_taskid},
                -10, 10,
                None,
                None,
                None,
                None,  #(host_db, host_taskid, 0, vbm_dx, cbm_dx),
                0.2,  #0.2
                locpot_c2db=None,  #(c2db, c2db_uid, 0)
                is_vacuum_aligment_on_plot=True,
                edge_tol=(-0.025, -0.025), # defect state will be picked only if it's above vbm by 0.025 eV and below
                # cbm by 0.025 eV
                ir_db=ir_db,
                ir_entry_filter={"prev_fw_taskid": defect_taskid},
            )
            tot, proj, d_df, levels, defect_levels = state
            level_info = d_df.to_dict("records")[0]
        except Exception as er:
            print(er)
            level_info = {}
            levels = {}
            defect_levels = {}
        return level_info, levels, defect_levels

class Host:
    def __init__(self, host_xlsx="host_2021-10-25"):
        if host_xlsx:
            self.host_df = IOTools(cwd=input_path, excel_file=host_xlsx).read_excel()

    def get_host_df(self):
        col = SCAN2dMat.collection
        filter = {
            "output.bandgap":{"$gte":1},
            #             "c2db_info.gap_hse_nosoc":{"$gte":1.5},
            "task_label": "SCAN_nscf line"
        }
        show = {
            "formula":1,
            "class":1,
            "ehull":1,
            "spacegroup":1,
            "gap_hse_nosoc":1
        }

        hosts = []
        host_candidates = list(col.find(filter))
        for e in host_candidates[:]:
            print(e["task_id"])
            print(e["task_id"], e["c2db_info"]["uid"])
            sym_data = e["sym_data"]
            unique_sites_from_wy = sym_data["unique_wyckoff"]
            good_ir = sym_data["good_ir_info"]

            print(e["c2db_info"]["uid"])
            d = {
                "host_taskid": e["task_id"],
                "species": tuple(unique_sites_from_wy["species"]),
                "reduced_species": tuple(good_ir["species"]),
                "site_sym": tuple(unique_sites_from_wy["site_sym"]),
                "reduced_site_sym": tuple(good_ir["site_sym"]),
                "wyckoffs": tuple(unique_sites_from_wy["wyckoffs"]),
                "reduced_wy": tuple(good_ir["wyckoffs"]),
                "site_idx": tuple(unique_sites_from_wy["site_idx"]),
                "reduced_site_idx": tuple(good_ir["site_idx"]),
                "prototype": e["c2db_info"]["prototype"],
                "formula": e["c2db_info"]["formula"],
                "c2db_spacegroup": e["c2db_info"]["spacegroup"],
                "spg": sym_data["pmg_spg"],
                "spg_number": sym_data["pmg_spg_number"],
                "pg": sym_data["pmg_pg"],
                "gap_scan": e["output"]["bandgap"],
                "c2db_gap_hse": e["c2db_info"]["gap_hse"],
                "c2db_gap_hse_nosoc": e["c2db_info"]["gap_hse_nosoc"],
                "c2db_soc_intensity": e["c2db_info"]["gap_hse_nosoc"] - e["c2db_info"]["gap_hse"],
                "c2db_ehull": e["c2db_info"]["ehull"],
                "c2db_uid":e["c2db_info"]["uid"],
            }
            d.update(e["band_edges"])
            hosts.append(d)
        # print(bad_pg)
        self.host_df = pd.DataFrame(hosts)
        self.host_df.fillna("None", inplace=True)
        print(self.host_df)
        IOTools(cwd=save_xlsx_path, pandas_df=self.host_df).to_excel("host")

    def get_host_gp(self):
        # all_triplet_df = defect_df.loc[~defect_df["reduced_site_sym"].isin([("-6m2", "-6m2"), ("3m", "3m"), ("4mm", "4mm")]), :]

        grouping = ["prototype", "spg", "reduced_site_sym",
                    "is_vbm_cbm_from_same_element", "reduced_species",
                    "vbm_up_max_el",
                    "vbm_up_max_proj_orbital",
                    "vbm_down_max_el",
                    "vbm_down_max_proj_orbital", "cbm_up_max_el", "cbm_up_max_proj_orbital",
                    "cbm_down_max_el", "cbm_down_max_proj_orbital"]
        agg_func = {"c2db_uid": ["count", "unique"], "host_taskid": ["unique"]}
        # agg_func.update(level_info)
        host_gp = self.host_df.groupby(grouping).agg(agg_func)
        IOTools(cwd=save_xlsx_path, pandas_df=host_gp).to_excel("host_gp", index=True)

class Defect:
    def __init__(self, defect_xlsx=None):
        if defect_xlsx:
            self.defect_df = IOTools(cwd=input_path, excel_file=defect_xlsx).read_excel()
    def get_defect_df(self):
        data = []
        col = SCAN2dDefect.collection
        for e in list(col.find({"task_label": "SCAN_scf", "task_id": 751}))[:]:
            print(e["task_id"])
            host_c2db_info = e["host_info"]["c2db_info"]
            for field in ["spacegroup", "pmg_point_gp", "irreps", "formula"]:
                if host_c2db_info.get(field):
                    host_c2db_info.pop(field)

            host_band_edges = e["host_info"]["scan_bs"]["band_edges"]
            for field in ["vbm_up_proj_on_el", "vbm_up_orbital_proj_on_el", "vbm_down_proj_on_el",
                          "vbm_down_orbital_proj_on_el", "cbm_up_proj_on_el", "cbm_up_orbital_proj_on_el"]:
                if host_band_edges.get(field):
                    host_band_edges.pop(field)

            host_sym_data = e["host_info"]["sym_data"]
            host_sym_data.update({"reduced_site_sym": tuple(e["host_info"]["sym_data"]["good_ir_info"]["site_sym"]),
                                  "reduced_site_specie": tuple(e["host_info"]["sym_data"]["good_ir_info"]["species"])
                                 })
            for field in ["unique_wyckoff", "good_ir_info"]:
                host_sym_data.pop(field)

            is_nelect_even = None
            if e["input"]["incar"]["NELECT"] % 2 == 0:
                is_nelect_even = True
            else:
                is_nelect_even = False

            host_pot_a, defect_pot_a = Potential(e["task_id"], SCAN2dDefect, SCAN2dMat).linear_fit_potential()

            info = {
                "task_id": e["task_id"],
                "host_taskid": e["pc_from_id"],
                "gap_scan": e["host_info"]["scan_bs"]["bandgap"],
                "defect_name": e["defect_name"],
                "defect_type": e["defect_entry"]["defect_type"],
                "charge": e["charge_state"],
                "mag": e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"],
                "chemsys": e["chemsys"],
                "is_nelect_even": is_nelect_even,
                "is_host_pot_steep": abs(host_pot_a[0]) > 5e-4,
                "is_defect_pot_steep": abs(defect_pot_a[0]) > 5e-4,
                "host_pot_a": host_pot_a[0],
                "defect_pot_a": defect_pot_a[0]
            }
            for host_info in [host_c2db_info, host_band_edges, host_sym_data]:
                info.update(host_info)

            d_df, levels = Tools.extract_defect_levels(e["task_id"])
            info.update(d_df)
            info.update(levels)
            data.append(info)


        self.defect_df = pd.DataFrame(data)
        self.defect_df.fillna("None", inplace=True)

        print(self.defect_df)
        IOTools(cwd=save_xlsx_path, pandas_df=self.defect_df).to_excel("defect")

    def get_defect_df_v2(self):
        data = []
        col = SCAN2dDefect.collection
        condition = {"task_label": "SCAN_scf", "task_id": {"$lte": 6084}}

        for e in list(col.find(condition))[:]:
            print(e["task_id"])
            host_c2db_info = e["host_info"]["c2db_info"]
            for field in ["spacegroup", "pmg_point_gp", "irreps", "formula"]:
                if host_c2db_info.get(field):
                    host_c2db_info.pop(field)

            host_band_edges = e["host_info"]["scan_bs"]["band_edges"]
            for field in ["vbm_up_proj_on_el", "vbm_up_orbital_proj_on_el", "vbm_down_proj_on_el",
                          "vbm_down_orbital_proj_on_el", "cbm_up_proj_on_el", "cbm_up_orbital_proj_on_el"]:
                if host_band_edges.get(field):
                    host_band_edges.pop(field)

            host_sym_data = e["host_info"]["sym_data"]
            host_sym_data.update({"reduced_site_sym": tuple(e["host_info"]["sym_data"]["good_ir_info"]["site_sym"]),
                                  "reduced_site_specie": tuple(e["host_info"]["sym_data"]["good_ir_info"]["species"])
                                 })
            for field in ["unique_wyckoff", "good_ir_info"]:
                host_sym_data.pop(field)

            is_nelect_even = None
            if e["input"]["incar"]["NELECT"] % 2 == 0:
                is_nelect_even = True
            else:
                is_nelect_even = False

            # host_pot_a, defect_pot_a = Potential(e["task_id"], SCAN2dDefect, SCAN2dMat).linear_fit_potential()

            info = {
                "task_id": e["task_id"],
                "host_taskid": e["pc_from_id"],
                "gap_scan": e["host_info"]["scan_bs"]["bandgap"],
                "defect_name": e["defect_name"],
                "defect_type": e["defect_entry"]["defect_type"],
                "charge": e["charge_state"],
                "mag": e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"],
                "chemsys": e["chemsys"],
                "is_nelect_even": is_nelect_even,
                # "is_host_pot_steep": abs(host_pot_a[0]) > 5e-4,
                # "is_defect_pot_steep": abs(defect_pot_a[0]) > 5e-4,
                # "host_pot_a": host_pot_a[0],
                # "defect_pot_a": defect_pot_a[0],
                "site_oxi_state": tuple([tuple(i) for i in e["host_info"]["scan_bs"]["site_oxi_state"]]),
                "number_NN": len(e["NN"]),
                "nbands": e["input"]["parameters"]["NBANDS"]

            }
            for host_info in [host_c2db_info, host_band_edges, host_sym_data]:
                info.update(host_info)

            d_df, levels, in_gpa_levels = Tools.extract_defect_levels_v2(e["task_id"])
            info.update(d_df)
            info.update(levels)
            info.update(in_gpa_levels)
            data.append(info)

        self.defect_df = pd.DataFrame(data)
        self.defect_df.fillna("None", inplace=True)
        self.defect_df.replace({"up_tran_en": "None"}, 0, inplace=True)
        self.defect_df.replace({"dn_tran_en": "None"}, 0, inplace=True)

        print(self.defect_df)
        print(os.getcwd())
        IOTools(cwd=save_xlsx_path, pandas_df=self.defect_df).to_excel("defect")

    def get_defect_df_v2_hse(self):
        data = []
        col = HSEQubitDefect.collection
        condition = {"task_label": "HSE_scf", "nupdown_set": 2}
        for e in list(col.find(condition))[:]:
            print(e["task_id"])
            e_from_scan = SCAN2dDefect.collection.find_one({"task_id": e["taskid_Scan2dDefect"]})
            host_c2db_info = e_from_scan["host_info"]["c2db_info"]
            for field in ["spacegroup", "pmg_point_gp", "irreps", "formula"]:
                if host_c2db_info.get(field):
                    host_c2db_info.pop(field)

            host_band_edges = e_from_scan["host_info"]["scan_bs"]["band_edges"]
            for field in ["vbm_up_proj_on_el", "vbm_up_orbital_proj_on_el", "vbm_down_proj_on_el",
                          "vbm_down_orbital_proj_on_el", "cbm_up_proj_on_el", "cbm_up_orbital_proj_on_el"]:
                if host_band_edges.get(field):
                    host_band_edges.pop(field)

            host_sym_data = e_from_scan["host_info"]["sym_data"]
            host_sym_data.update({"reduced_site_sym": tuple(e_from_scan["host_info"]["sym_data"]["good_ir_info"]["site_sym"]),
                                  "reduced_site_specie": tuple(e_from_scan["host_info"]["sym_data"]["good_ir_info"]["species"])
                                 })
            for field in ["unique_wyckoff", "good_ir_info"]:
                host_sym_data.pop(field)

            is_nelect_even = None
            if e["input"]["incar"]["NELECT"] % 2 == 0:
                is_nelect_even = True
            else:
                is_nelect_even = False


            info = {
                "task_id": e["task_id"],
                "host_taskid": e["pc_from"].split("/")[-1],
                "gap_scan": e_from_scan["host_info"]["scan_bs"]["bandgap"],
                "defect_name": e["defect_entry"]["name"],
                "defect_type": e["defect_entry"]["defect_type"],
                "charge": e["charge_state"],
                "mag": round(e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]),
                "chemsys": e["chemsys"],
                "is_nelect_even": is_nelect_even,
                "site_oxi_state": tuple([tuple(i) for i in e_from_scan["host_info"]["scan_bs"]["site_oxi_state"]]),
                "number_NN": len(e["NN"]),
                "nbands": e["input"]["parameters"]["NBANDS"]
            }
            for host_info in [host_c2db_info, host_band_edges, host_sym_data]:
                info.update(host_info)

            d_df, levels, in_gpa_levels = Tools.extract_defect_levels_v2_hse(e["task_id"])

            info.update(d_df)
            info.update(levels)
            info.update(in_gpa_levels)
            data.append(info)

        self.defect_df = pd.DataFrame(data)
        self.defect_df.fillna("None", inplace=True)
        self.defect_df.replace({"up_tran_en": "None"}, 0, inplace=True)
        self.defect_df.replace({"dn_tran_en": "None"}, 0, inplace=True)

        print(self.defect_df)
        print(os.getcwd())
        IOTools(cwd=save_xlsx_path, pandas_df=self.defect_df).to_excel("hse_qubit")

    def get_screened_defect_df(self):
        defect_df = self.defect_df.copy()
        triplet_plt_df = defect_df.loc[(defect_df["mag"] == 2) & (defect_df["level_cbm"] != "None")
                                       & (defect_df["level_vbm"] != "None"), :]
        def find_certain_df():
            screened_df = defect_df.loc[defect_df["fworker"] == "nersc", :]
            return screened_df

        def triplet_plt():
            screened_df = defect_df.loc[(defect_df["mag"] == 2)]
            return screened_df

        def qubit_candidate():
            # triplet_df.replace({"up_tran_en": "None"}, 0, inplace=True)
            # triplet_df.replace({"is_vbm_cbm_from_same_element": "True"}, 0, inplace=True)
            # screened_df = triplet_plt_df.loc[triplet_plt_df["level_cat"] == "4", :]
            condition = (defect_df["level_cbm"] != "None") & (defect_df["level_vbm"] != "None")
            screened_df = triplet_plt().loc[condition]
            screened_df = screened_df.loc[
                          (screened_df["up_tran_en"] >= 0.5) |
                          (screened_df["dn_tran_en"] >= 0.5), :]
            return screened_df

        def qubit_and_follow_level_hypth():
            screened_df = triplet_plt_df.loc[
                          (triplet_plt_df["up_tran_en"] >= 0.5) |
                          (triplet_plt_df["dn_tran_en"] >= 0.5), :]
            screened_df = screened_df.loc[screened_df["is_follow_level_hypth"] != True, :]
            return screened_df

        def vacancy():
            screened_df = triplet_plt_df.loc[(triplet_plt_df["defect_type"] == "vacancy"), :]
            return screened_df

        def band_edges():
            screened_df = triplet_plt_df.loc[(triplet_plt_df["is_vbm_cbm_from_same_element"] == False), :]
            return screened_df

        return triplet_plt()

    def get_defect_df_gp(self, df, output_gp_name):
        # all_triplet_df = defect_df.loc[~defect_df["reduced_site_sym"].isin([("-6m2", "-6m2"), ("3m", "3m"), ("4mm", "4mm")]), :]

        grouping = ["prototype", "pmg_spg", "reduced_site_sym", "reduced_site_specie", "uid", "host_taskid",
                    "site_oxi_state"]
        # grouping.extend([
        #     "vbm_up_max_el", "vbm_up_max_proj_orbital", "vbm_down_max_el",
        #     "vbm_down_max_proj_orbital", "cbm_up_max_el", "cbm_up_max_proj_orbital",
        #     "cbm_down_max_el", "cbm_down_max_proj_orbital", "is_vbm_cbm_from_same_element"
        # ])
        grouping.extend([
            "vbm_max_el", "vbm_max_proj_orbital", "cbm_max_el", "cbm_max_proj_orbital", "is_vbm_cbm_from_same_element"
        ])
        grouping.extend(["defect_type", "charge", "defect_name", "mag"])
        grouping.extend(["level_source_specie", "tran_top", "level_gap", "level_from_edge",
                         "is_tran_top_near_vbm", "is_follow_level_hypth"])
        agg_func = {"up_tran_en": ["unique"], "dn_tran_en": ["unique"], "up_tran_bottom": "unique",
                    "dn_tran_bottom": "unique", "task_id": ["count", "unique"], }
        # agg_func.update(level_info)
        df = df.groupby(grouping).agg(agg_func)
        IOTools(cwd=save_xlsx_path, pandas_df=df).to_excel(output_gp_name, index=True)
        return df

    @classmethod
    def bar_plot_defect_levels(cls, defect_df):
        all_protos = ["BN"]#"AgBr3", "BN", "BiI3", "CH", "CdI2", "FeSe", "GaS", "GaSe", "GeS2", "GeSe", "MoS2", "TiCl3"]
        # defect_df_copy = defect_df.loc[defect_df["charge"] == 0, :]
        for set_prototype in all_protos:
            defect_df_copy = defect_df.copy()
            print(set_prototype)
            conditions = (defect_df_copy["prototype"] == set_prototype)
            defect_df_copy = defect_df_copy.loc[conditions, :]
            if defect_df_copy.empty:
                print("defect_df_copy is empty!")
                continue
            host_taskids = defect_df_copy.loc[:, "host_taskid"].unique()
            uids = defect_df_copy.loc[:, "uid"].unique()
            print(host_taskids)
            for host_taskid, uid in zip(host_taskids,uids):
                plot_df = defect_df_copy.loc[(defect_df_copy["host_taskid"] == host_taskid) &
                                             (defect_df_copy["level_cbm"] != "None") &
                                             (defect_df_copy["level_vbm"] != "None"), :]
                if plot_df.empty:
                    continue

                defect_levels = []
                taskid_labels = []
                mags =[]
                defect_name_labels = []
                charge_labels = []
                for taskid in plot_df["task_id"].tolist():
                    defect = plot_df.loc[(plot_df["task_id"] == taskid) , :]
                    mag = plot_df["mag"].iloc[0]
                    print("defect_taskid:{}, host_taskid:{}".format(defect["task_id"].iloc[0],
                                                                    defect["host_taskid"].iloc[0]))
                    # levels = extract_defect_levels(defect["d_taskid"])
                    defect_levels.append(
                        {
                            "cbm": defect["level_cbm"].iloc[0],
                            "vbm": defect["level_vbm"].iloc[0],
                            "up": defect["1"].iloc[0],
                            "dn": defect["-1"].iloc[0],
                            "up_deg": defect["level_up_deg"].iloc[0] ,
                            "dn_deg": defect["level_dn_deg"].iloc[0],
                            "up_ir": defect["level_up_ir"].iloc[0],
                            "dn_ir": defect["level_dn_ir"].iloc[0]
                        }
                    )
                    taskid_labels.append(defect["task_id"].iloc[0])
                    mags.append(defect["mag"].iloc[0])
                    defect_name_labels.append(defect["defect_name"].iloc[0])
                    charge_labels.append(defect["charge"].iloc[0])




                fig, ax = plt.subplots(figsize=(12,11), dpi=300)
                print(np.arange(0, len(defect_levels)*2, 2))
                print(defect_levels)
                for defect, x in zip(defect_levels, np.arange(0, len(defect_levels)*5, 5)):
                    cbm, vbm, up_deg, dn_deg, up_ir, dn_ir = defect["cbm"], defect["vbm"], defect["up_deg"], \
                                                             defect["dn_deg"], defect["up_ir"], defect["dn_ir"]
                    print(defect["up"].keys(), defect["dn"].keys())
                    if not up_ir:
                        up_ir = [None for i in defect["up"].keys()]
                    if not up_deg:
                        up_deg = [None for i in defect["up"].keys()]

                    if not dn_ir:
                        dn_ir = [None for i in defect["dn"].keys()]
                    if not dn_deg:
                        dn_deg = [None for i in defect["dn"].keys()]

                    for edge in [cbm, vbm]:
                        ax.bar(x, vbm--15, 4, -15, color="deepskyblue")
                        ax.bar(x, 5-cbm, 4, cbm, color="orange")

                    dx = 0.35
                    for level, occupied, deg, ir in zip(defect["up"].keys(), defect["up"].values(), up_deg, up_ir):
                        # dx += 0.2
                        print(ir)
                        ax.text(x-2, level, ir)
                        color = None
                        if deg == 1:
                            color = "red"
                        elif deg == 2:
                            color = "blue"
                        elif deg == 3:
                            color = "green"
                        elif deg == 4:
                            color = "yellow"
                        else:
                            color = "black"

                        if occupied:
                            ax.hlines(level, x-2, x-0.5, colors=color)
                        else:
                            ax.hlines(level, x-2, x-0.5, colors=color)
                            ax.text(x-0.5, level, "*")

                    dx = 1
                    for level, occupied, deg, ir in zip(defect["dn"].keys(), defect["dn"].values(), dn_deg, dn_ir):
                        ax.text(x+2, level, ir)
                        color = None
                        if deg == 1:
                            color = "red"
                        elif deg == 2:
                            color = "blue"
                        elif deg == 3:
                            color = "green"
                        elif deg == 4:
                            color = "yellow"
                        else:
                            color = "black"

                        if occupied:
                            ax.hlines(level, x+0.5, x+2, colors=color)
                            print(level)
                        else:
                            ax.hlines(level, x+0.5, x+2, colors=color)
                            ax.text(x+0.5, level, "*")

                ax.set_ylim(-15, 5)
                ax.yaxis.set_minor_locator(AutoMinorLocator(5))
                for tick in ax.get_yticklabels():
                    tick.set_fontsize(15)


                ax.tick_params(axis="x", bottom=False, labelbottom=True)
                ax.set_xticks(np.arange(0, len(defect_levels)*5, 5))
                print(taskid_labels, mags)
                ax.set_xticklabels(["{}\nmag:{}".format(taskid_labels[i], mags[i]) for i in range(len(mags))],
                                   fontdict={"fontsize": 10}, rotation=0)

                for x, name, q in zip(np.arange(0, len(defect_levels)*5, 5), defect_name_labels, charge_labels):
                    ax.text(x-2, 5+0.1, "{}-{}".format(q, name), rotation=0, size=10)

                ax.set_ylabel("Energy relative to vacuum (eV)", fontsize=15)
                ax.text(-0.1, 1.05, "{}-{}".format(uid, host_taskid), size=30, transform=ax.transAxes)

                # fig.savefig(
                #     os.path.join(
                #         save_plt_path,
                #         "triplets_levels(incorrected)",
                #         "triplets_polarized_hosts",
                #         "neutral_triplets_unpolarized_hosts",
                #         "2021-10-27",
                #         "proto_{}_triplets.png".format(set_prototype)
                #     )
                # )
                os.makedirs(os.path.join(save_plt_path, set_prototype), exist_ok=True)
                fig.savefig(
                    os.path.join(
                        save_plt_path,
                        set_prototype,
                        "host_taskid-{}.png".format(host_taskid)
                    )
                )

    @classmethod
    def bar_plot_defect_levels_v2(cls, defect_df):
        plt.style.use(["science", "notebook", "grid"])
        all_protos = ["AgBr3", "BN", "BiI3", "CH", "CdI2", "FeSe", "GaS", "GaSe", "GeS2", "GeSe", "MoS2",
                      "TiCl3"]
        # defect_df_copy = defect_df.loc[defect_df["charge"] == 0, :]
        for set_prototype in defect_df.prototype.unique()[:]:
            defect_df_copy = defect_df.copy()
            print(set_prototype)
            conditions = (defect_df_copy["prototype"] == set_prototype)
            defect_df_copy = defect_df_copy.loc[conditions, :]
            if defect_df_copy.empty:
                print("defect_df_copy is empty!")
                continue
            host_taskids = defect_df_copy.loc[:, "host_taskid"].unique()
            uids = defect_df_copy.loc[:, "uid"].unique()
            print(host_taskids)
            for host_taskid, uid in zip(host_taskids,uids):
                plot_df = defect_df_copy.loc[(defect_df_copy["host_taskid"] == host_taskid) &
                                             (defect_df_copy["level_cbm"] != "None") &
                                             (defect_df_copy["level_vbm"] != "None"), :]
                if plot_df.empty:
                    continue

                defects = []
                taskid_labels = []
                mags =[]
                defect_name_labels = []
                charge_labels = []
                for taskid in plot_df["task_id"].tolist():
                    defect = plot_df.loc[(plot_df["task_id"] == taskid), :]
                    mag = plot_df["mag"].iloc[0]
                    print("=="*20)
                    print("defect_taskid:{}, host_taskid:{}".format(defect["task_id"].iloc[0],
                                                                    defect["host_taskid"].iloc[0]))

                    defects.append(
                        {
                            "cbm": defect["level_cbm"].iloc[0],
                            "vbm": defect["level_vbm"].iloc[0],
                            "level_edge_ir": defect["level_edge_ir"].iloc[0],
                            "up": defect["up_in_gap_level"].iloc[0],
                            "dn": defect["dn_in_gap_level"].iloc[0],
                            "up_deg": defect["up_in_gap_deg"].iloc[0] ,
                            "dn_deg": defect["dn_in_gap_deg"].iloc[0],
                            "up_ir": defect["up_in_gap_ir"].iloc[0],
                            "dn_ir": defect["dn_in_gap_ir"].iloc[0],
                            "up_occ": defect["up_in_gap_occ"].iloc[0],
                            "dn_occ": defect["dn_in_gap_occ"].iloc[0],
                            "up_id": defect["up_in_gap_band_id"].iloc[0],
                            "dn_id": defect["dn_in_gap_band_id"].iloc[0],
                            "up_band": defect["up_in_gap_band_index"].iloc[0],
                            "dn_band": defect["dn_in_gap_band_index"].iloc[0],

                            "vbm_max_el": defect["vbm_max_el"].iloc[0],
                            "vbm_max_proj_orbital": defect["vbm_max_proj_orbital"].iloc[0],
                            "cbm_max_el": defect["cbm_max_el"].iloc[0],
                            "cbm_max_proj_orbital": defect["cbm_max_proj_orbital"].iloc[0],

                            "up_tran_bottom": defect["up_tran_bottom"].iloc[0] + defect["level_vbm"].iloc[0] if
                            defect["up_tran_bottom"].iloc[0] != "None" else None,
                            "up_tran_en": defect["up_tran_en"].iloc[0] if defect["up_tran_en"].iloc[0] else None,
                            "dn_tran_bottom": defect["dn_tran_bottom"].iloc[0] + defect["level_vbm"].iloc[0] if
                            defect["dn_tran_bottom"].iloc[0] != "None" else None,
                            "dn_tran_en": defect["dn_tran_en"].iloc[0] if defect["dn_tran_en"].iloc[0] else None,

                        }
                    )

                    taskid_labels.append(defect["task_id"].iloc[0])
                    mags.append(defect["mag"].iloc[0])
                    defect_name_labels.append(defect["defect_name"].iloc[0])
                    charge_labels.append(defect["charge"].iloc[0])

                vbm_lim, cbm_lim = plot_df["level_vbm"].median(), plot_df["level_cbm"].median()

                fig, ax = plt.subplots(figsize=(12,11), dpi=300)
                fig_y_height = 5, 5
                ax.set_ylim(vbm_lim-fig_y_height[0], cbm_lim+fig_y_height[1])

                print(np.arange(0, len(defects)*2, 2))
                print(defects)
                for defect, x in zip(defects, np.arange(0, len(defects)*5, 5)):
                    vbm = defect["vbm"]
                    cbm = defect["cbm"]
                    level_edge_ir = defect["level_edge_ir"],
                    up = defect["up"]
                    dn = defect["dn"]
                    up_deg = defect["up_deg"]
                    dn_deg = defect["dn_deg"]
                    up_ir = defect["up_ir"]
                    dn_ir = defect["dn_ir"]
                    up_occ = defect["up_occ"]
                    dn_occ = defect["dn_occ"]
                    up_id = defect["up_id"]
                    dn_id = defect["dn_id"]
                    up_band = defect["up_band"]
                    dn_band = defect["dn_band"]
                    vbm_max_el = defect["vbm_max_el"]
                    vbm_max_proj_orbital = defect["vbm_max_proj_orbital"]
                    cbm_max_el = defect["cbm_max_el"]
                    cbm_max_proj_orbital = defect["cbm_max_proj_orbital"]
                    up_tran_bottom = defect["up_tran_bottom"]
                    up_tran_en = defect["up_tran_en"]
                    dn_tran_bottom = defect["dn_tran_bottom"]
                    dn_tran_en = defect["dn_tran_en"]

                    if not up_ir:
                        up_ir = ["" for i in up]
                    if not up_deg:
                        up_deg = ["" for i in up]
                    if not up_occ:
                        up_occ = ["" for i in up]

                    if not dn_ir:
                        dn_ir = ["" for i in dn]
                    if not dn_deg:
                        dn_deg = ["" for i in dn]
                    if not dn_occ:
                        dn_occ = ["" for i in dn]

                    ax.bar(x, vbm-(vbm_lim-fig_y_height[0]), 4, (vbm_lim-fig_y_height[0]), color="deepskyblue")
                    ax.bar(x, (cbm_lim+fig_y_height[1])-cbm, 4, cbm, color="orange")

                    dx = 0.35
                    for level, occ, deg, id in zip(up, up_occ, up_deg, up_id):
                        print(occ)
                        # dx += 0.2
                        ax.text(x-2, level, id)
                        color = None
                        if deg == 1:
                            color = "red"
                        elif deg == 2:
                            color = "blue"
                        elif deg == 3:
                            color = "green"
                        elif deg == 4:
                            color = "yellow"
                        else:
                            color = "black"

                        ax.hlines(level, x-2, x-0.5, colors=color)
                        if occ == 0:
                            ax.text(x-0.5, level, "*")
                        elif round(occ, 1) == 0.5:
                            ax.text(x-0.5, level, "%")

                    up_info = '\n'.join(
                        [
                            "{}/ {}/ {}/ {}/ {}/ {}".format(level, occ, deg, ir, band_id, band_index) for
                            level, occ, deg, ir, band_id, band_index in zip(up, up_occ, up_deg, up_ir, up_id, up_band)
                        ]
                    )
                    ax.text(x-1.5, cbm+0.45, up_info, bbox=dict(facecolor='white', edgecolor='black'), size=8)
                    ax.text(x, cbm+3, "{}-{}-{}".format(cbm_max_el, cbm_max_proj_orbital, level_edge_ir[0][1]),
                            bbox=dict(
                        facecolor='white', edgecolor='black'), size=12)
                    if up_tran_en and up_tran_bottom:
                        ax.arrow(x=x-1.25, y=up_tran_bottom, dx=0, dy=up_tran_en, width=0.02,
                                 length_includes_head=True, facecolor="gray")

                    dx = 1
                    for level, occ, deg, id in zip(dn, dn_occ, dn_deg, dn_id):
                        ax.text(x+2, level, id)
                        color = None
                        if deg == 1:
                            color = "red"
                        elif deg == 2:
                            color = "blue"
                        elif deg == 3:
                            color = "green"
                        elif deg == 4:
                            color = "yellow"
                        else:
                            color = "black"

                        ax.hlines(level, x+0.5, x+2, colors=color)
                        if occ == 0:
                            ax.text(x+0.5, level, "*")
                        elif round(occ, 1) == 0.5:
                            ax.text(x+0.5, level, "%")

                    dn_info = '\n'.join(
                        [
                            "{}/ {}/ {}/ {}/ {}/ {}".format(level, occupied, deg, ir, band_id, band_index) for
                            level, occupied, deg, ir, band_id, band_index in zip(dn, dn_occ, dn_deg, dn_ir, dn_id,
                                                                                 dn_band)
                        ]
                    )
                    ax.text(x-1.5, vbm_lim-fig_y_height[0]+0.1, dn_info, bbox=dict(facecolor='white',
                                                                                   edgecolor='black'),
                            size=8)

                    ax.text(x, vbm_lim-fig_y_height[0]+4, "{}-{}-{}".format(vbm_max_el, vbm_max_proj_orbital,
                                                                            level_edge_ir[0][0]),
                            bbox=dict(facecolor='white', edgecolor='black'), size=12)
                    if dn_tran_en and dn_tran_bottom:
                        ax.arrow(x=x+1.25, y=dn_tran_bottom, dx=0, dy=dn_tran_en, width=0.02,
                                 length_includes_head=True, facecolor="gray")

                ax.yaxis.set_minor_locator(AutoMinorLocator(5))
                for tick in ax.get_yticklabels():
                    tick.set_fontsize(15)

                ax.tick_params(axis="x", bottom=False, labelbottom=True)
                ax.tick_params(axis="y", right=False, which="both")

                ax.set_xticks(np.arange(0, len(defects)*5, 5))
                print(taskid_labels, mags)
                ax.set_xticklabels(["{}\nmag:{}".format(taskid_labels[i], mags[i]) for i in range(len(mags))],
                                   fontdict={"fontsize": 10}, rotation=0)

                for x, name, q in zip(np.arange(0, len(defects)*5, 5), defect_name_labels, charge_labels):
                    ax.text(x-2, cbm_lim+fig_y_height[1]+0.1, "{}-{}".format(q, name), rotation=0, size=10)

                ax.set_ylabel("Energy relative to vacuum (eV)", fontsize=15)
                ax.text(-0.1, 1.05, "{}-{}".format(uid, host_taskid), size=30, transform=ax.transAxes)

                os.makedirs(os.path.join(save_plt_path, set_prototype), exist_ok=True)
                fig.savefig(
                    os.path.join(
                        save_plt_path,
                        set_prototype,
                        "host_taskid-{}.png".format(host_taskid)
                    )
                )
                plt.close()


    def statistics(self, df):
        df = df.loc[df["defect_type"] == "antisite", :]
        df = df.groupby("charge")
        print(df.agg({"task_id": ["count"]}))
        # print(df["task_id"].count())

class BackProcess:
    def __init__(self, input_df):
        self.input_df = input_df

    def df_to_excel(self, save_path=save_xlsx_path, excel_name="test"):
        IOTools(cwd=save_path, pandas_df=self.input_df).to_excel(excel_name)

    def add_number_occ_electrons(self):
        """
        Regarding electron counting
        Q: charge state for triplets
        Vac_e: vacancy-type number of electrons
        Anti-e: antisite-type number of electrons (0 if no antisite)
        Occ-e: Vac_e + Anti-e
        Triplet-e: triplet-configuration number of electrons (from calculation)
        Q = Occ-e - Triplet-e
        :return:
        """
        def fun(x):
            oxi = x["site_oxi_state"]
            defect_name = x["defect_name"]
            NN = x["number_NN"]
            occ_e, vac_e, anti_e = None, None, None
            if "vac" in defect_name:
                vac_specie = defect_name.split("_")[-1]
                for site_oxi in oxi:
                    if site_oxi[0] == vac_specie:
                        continue
                    NN_oxi = site_oxi[1]
                    vac_e = abs(NN_oxi)*NN
                    anti_e = 0
                    occ_e = vac_e + anti_e
                    break

            else:
                ant_specie = defect_name.split("_")[2]
                orig_specie = defect_name.split("_")[-1]
                for site_oxi in oxi:
                    if site_oxi[0] != ant_specie:
                        continue
                    NN_oxi = site_oxi[1]
                    vac_e = abs(NN_oxi)*(NN-1)
                    ant_val, orig_val = None, None
                    try:
                        ant_val = Element(ant_specie).valence[-1]
                    except ValueError:
                        print("valence not found:{} or {}".format(ant_specie, orig_specie))
                        group, row = Element(ant_specie).group, Element(ant_specie).row
                        try:
                            ant_val = Element.from_row_and_group(row+1, group).valence[-1]
                        except:
                            ant_val = Element.from_row_and_group(row-1, group).valence[-1]
                    try:
                        orig_val = Element(orig_specie).valence[-1]
                    except ValueError:
                        group, row = Element(orig_specie).group, Element(orig_specie).row
                        try:
                            orig_val = Element.from_row_and_group(row+1, group).valence[-1]
                        except:
                            orig_val = Element.from_row_and_group(row-1, group).valence[-1]
                    anti_e = ant_val - orig_val
                    occ_e = vac_e + anti_e
                    break
            return (vac_e, anti_e, occ_e)
        defect_df = self.input_df
        defect_df["occ_e"] = defect_df.apply(lambda x: fun(x), axis=1)
        defect_df["triplet_e"] = defect_df.apply(lambda x: x["occ_e"][-1] - x["charge"], axis=1)
        defect_df.fillna("None", inplace=True)

    def add_band_edges_and_defects(self):
        def fun(x):
            is_vbm_cbm = x["is_vbm_cbm_from_same_element"]
            vbm_up_max_el = x["vbm_up_max_el"] if x["vbm_up_max_el"] != "None" else None
            vbm_up_max_proj_orbital = x["vbm_up_max_proj_orbital"] if x["vbm_up_max_proj_orbital"] != "None" else \
                None

            vbm_down_max_el = x["vbm_down_max_el"] if x["vbm_down_max_el"] != "None" else None
            vbm_down_max_proj_orbital = x["vbm_down_max_proj_orbital"] if x["vbm_down_max_proj_orbital"] != \
                                                                          "None" else None

            vbm_max_el = vbm_up_max_el or vbm_down_max_el
            vbm_max_proj_orbital= vbm_up_max_proj_orbital or vbm_down_max_proj_orbital

            cbm_up_max_el = x["cbm_up_max_el"] if x["cbm_up_max_el"] != "None" else None
            cbm_up_max_proj_orbital = x["cbm_up_max_proj_orbital"] if x["cbm_up_max_proj_orbital"] != "None" else None
            cbm_down_max_el = x["cbm_down_max_el"] if x["cbm_down_max_el"] != "None" else None
            cbm_down_max_proj_orbital = x["cbm_down_max_proj_orbital"] if x["cbm_down_max_proj_orbital"] !="None" else None
            cbm_max_el = cbm_up_max_el or cbm_down_max_el
            cbm_max_proj_orbital= cbm_up_max_proj_orbital or cbm_down_max_proj_orbital



            up_tran_top = x["up_tran_top"] if x["up_tran_top"] != "None" else 0
            dn_tran_top = x["dn_tran_top"] if x["dn_tran_top"] != "None" else 0
            tran_top = max(up_tran_top, dn_tran_top)
            level_gap = None
            try:
                level_gap = x["level_cbm"] - x["level_vbm"]
            except Exception:
                level_gap = 0

            level_source_specie = None
            if x["defect_type"] == "antisite":
                level_source_specie = x["defect_name"].split("_")[2]
            elif x["defect_type"] == "vacancy":
                for el in x["chemsys"].split("-"):
                    if x["defect_name"].split("_")[2] != el:
                        level_source_specie = el

            is_tran_top_near_vbm = None
            if level_gap == 0:
                is_tran_top_near_vbm = "None"
            elif tran_top >= 0.5*level_gap:
                is_tran_top_near_vbm = False
            else:
                is_tran_top_near_vbm = True


            level_from_edge = None
            if vbm_max_el == cbm_max_el:
                if vbm_max_el == level_source_specie:
                    level_from_edge = "both"
                else:
                    level_from_edge = "None"
            else:
                if vbm_max_el == level_source_specie:
                    level_from_edge = "vbm"
                elif cbm_max_el == level_source_specie:
                    level_from_edge = "cbm"
                else:
                    level_from_edge = "None"

            is_follow_level_hypth = None
            if level_from_edge == "both":
                if is_tran_top_near_vbm is True:
                    is_follow_level_hypth = False
                elif is_tran_top_near_vbm is False:
                    is_follow_level_hypth = True
                else:
                    is_follow_level_hypth = "None"

            elif level_from_edge == "cbm":
                if is_tran_top_near_vbm is False:
                    is_follow_level_hypth = True
                elif is_tran_top_near_vbm is True:
                    is_follow_level_hypth = False
                else:
                    is_follow_level_hypth = "None"

            elif level_from_edge == "vbm":
                if is_tran_top_near_vbm is True:
                    is_follow_level_hypth = True
                elif is_tran_top_near_vbm is False:
                    is_follow_level_hypth = False
                else:
                    is_follow_level_hypth = "None"
            else:
                is_follow_level_hypth = "None"

            return vbm_max_el, cbm_max_el, tran_top, level_gap, level_from_edge, is_tran_top_near_vbm, \
                   is_follow_level_hypth, level_source_specie, vbm_max_proj_orbital, cbm_max_proj_orbital

        defect_df = self.input_df
        print(defect_df)
        defect_df["vbm_max_el"] = defect_df.apply(lambda x: fun(x)[0], axis=1)
        defect_df["cbm_max_el"] = defect_df.apply(lambda x: fun(x)[1], axis=1)
        defect_df["tran_top"] = defect_df.apply(lambda x: fun(x)[2], axis=1)
        defect_df["level_gap"] = defect_df.apply(lambda x: fun(x)[3], axis=1)
        defect_df["level_from_edge"] = defect_df.apply(lambda x: fun(x)[4], axis=1)
        defect_df["is_tran_top_near_vbm"] = defect_df.apply(lambda x: fun(x)[5], axis=1)
        defect_df["is_follow_level_hypth"] = defect_df.apply(lambda x: fun(x)[6], axis=1)
        defect_df["level_source_specie"] = defect_df.apply(lambda x: fun(x)[7], axis=1)
        defect_df["vbm_max_proj_orbital"] = defect_df.apply(lambda x: fun(x)[8], axis=1)
        defect_df["cbm_max_proj_orbital"] = defect_df.apply(lambda x: fun(x)[9], axis=1)
        defect_df.fillna("None", inplace=True)

    def add_level_category(self):
        def fun(x):
            level_from_edge = x["level_from_edge"]
            is_tran_top_near_vbm = x["is_tran_top_near_vbm"]
            vbm_max_el = x["vbm_max_el"]
            cbm_max_el = x["cbm_max_el"]
            level_source_specie = x["level_source_specie"]
            cat = None
            if level_from_edge == "both":
                cat = "1"
            elif level_from_edge == "vbm":
                if is_tran_top_near_vbm is True:
                    cat = "2"
                else:
                    cat = "2'"
            elif level_from_edge == "cbm":
                if is_tran_top_near_vbm is False:
                    cat = "3"
                else:
                    cat = "3'"
            else:
                if vbm_max_el == cbm_max_el:
                    if vbm_max_el != level_source_specie:
                        cat = "4"
                else:
                    cat = "None"

            return cat
        defect_df = self.input_df
        defect_df["level_cat"] = defect_df.apply(lambda x: fun(x), axis=1)
        defect_df.fillna("None", inplace=True)

    def add_transition_wavevlength(self):
        def fun(x): # in nm
            up_tran_en = x["up_tran_en"]
            dn_tran_en = x["dn_tran_en"]
            up_tran_wavelength = round(12400/up_tran_en/10, 0) if up_tran_en !=0 else "None"
            dn_tran_wavelength = round(12400/dn_tran_en/10, 0) if dn_tran_en !=0 else "None"
            return up_tran_wavelength, dn_tran_wavelength
        defect_df = self.input_df
        defect_df["up_tran_wavelength"] = defect_df.apply(lambda x: fun(x)[0], axis=1)
        defect_df["dn_tran_wavelength"] = defect_df.apply(lambda x: fun(x)[1], axis=1)
        defect_df.fillna("None", inplace=True)

    def add_hse_fworker(self):
        db = HSEQubitDefect
        def fun(x):
            taskid = x["task_id"]
            fworker = None
            try:
                entry = db.collection.find_one({"task_id": taskid})
                if "efrc" in entry["dir_name"]:
                    fworker = "efrc"
                elif "scratch" in entry["dir_name"]:
                    fworker = "owls"
                elif "global" in entry["dir_name"]:
                    fworker = "nersc"
            except Exception as err:
                print(err)
            return fworker
        defect_df = self.input_df
        defect_df["fworker"] = defect_df.apply(lambda x: fun(x), axis=1)
        defect_df.fillna("None", inplace=True)
        return defect_df

    def add_cdft_occupations(self):
        def fun(x):
            nbands = x["nbands"]
            up_lumo_band_index = x["up_tran_lumo_homo_band_index"]
            dn_lumo_band_index = x["dn_tran_lumo_homo_band_index"]
            up_lumo_homo_deg = x["up_tran_lumo_homo_deg"]
            dn_lumo_homo_deg = x["dn_tran_lumo_homo_deg"]

            cdft_occs = {}
            for lumo_band_index, lumo_homo_deg, spin in zip([up_lumo_band_index, dn_lumo_band_index],
                                                            [up_lumo_homo_deg, dn_lumo_homo_deg], ["1", "-1"]):
                if lumo_homo_deg == (1, 2) and lumo_band_index:
                    cdft_occ = "{}*1 1*0.5 1*0.5 1*1 {}*0".format(lumo_band_index[0]-3, nbands-lumo_band_index[0])
                elif lumo_homo_deg == (2, 1) and lumo_band_index:
                    cdft_occ = "{}*1 1*0 1*0.5 1*0.5 {}*0".format(lumo_band_index[0]-2,
                                                                  nbands-(lumo_band_index[0]+1))
                elif lumo_homo_deg == (1, 1) and lumo_band_index:
                    cdft_occ = "{}*1 1*0 1*1 {}*0".format(lumo_band_index[0]-2, nbands-lumo_band_index[0])
                elif lumo_homo_deg == (2, 2) and lumo_band_index:
                    cdft_occ = "{}*1 1*0.5 1*0.5 1*0.5 1*0.5 {}*0".format(lumo_band_index[0]-3,
                                                                          nbands-(lumo_band_index[0]+1))
                elif lumo_homo_deg == (3, 2) and lumo_band_index:
                    cdft_occ = "{}*1 1*0.5 1*0.5 1*0.333333 1*0.333333 1*0.333333 {}*0".format(
                        lumo_band_index[0]-3, nbands-(lumo_band_index[0]+2))
                else:
                    cdft_occ = None
                cdft_occs[spin] = cdft_occ
            return cdft_occs["1"], cdft_occs["-1"]

        defect_df = self.input_df
        defect_df["up_tran_cdft_occ"] = defect_df.apply(lambda x: fun(x)[0], axis=1)
        defect_df["dn_tran_cdft_occ"] = defect_df.apply(lambda x: fun(x)[1], axis=1)
        defect_df.fillna("None", inplace=True)


class CDFT:
    def __init__(self):
        self.cdft_db = HSEQubitCDFT
        self.calc_db = HSEQubitDefect
        self.foundation_df = None

    def update_entry_with_occ(self):
        from pymatgen.io.vasp.inputs import Incar
        # run in db1
        es = self.cdft_db.collection.find()
        main_path = '/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/cdft-pbe_pc/'
        for e in es:
            if type(e["cdft_occ"]) == dict:
                continue

            print(e["task_id"])
            src_entry, up_occ, dn_occ = None, None, None
            if e["task_label"] == "CDFT-D-HSE_scf":
                src_entry = self.cdft_db.collection.find_one({"task_id": e["prev_fw_taskid"]})

            else:
                src_entry = e

            dir_name = src_entry["dir_name"].split("/")[-1]
            path = os.path.join(main_path, dir_name)
            incar = Incar.from_file(os.path.join(path, "INCAR"))
            up_occ = incar["FERWE"]
            dn_occ = incar["FERDO"]

            self.cdft_db.collection.update_one(
                {"task_id": e["task_id"]},
                {"$set": {"cdft_occ": {"up": up_occ, "dn": dn_occ}}},
            )


    def get_data_sheet(self):
        es = self.cdft_db.collection.find()
        data = {"prototype": [], "pc_from": [], "gs_taskid": [], "task_id": [], "defect_name": [],
                "charge_state":[], "up_cdft_occ":[],"dn_cdft_occ": [], "task_label":[], "energy": []}
        for e in es:
            print(e["task_id"])
            gs_taskid = int(e["source_entry"].split("/")[-1])
            gs_e = self.calc_db.collection.find_one({"task_id": gs_taskid})
            data["task_id"].append(e["task_id"])
            data["defect_name"].append(gs_e["defect_entry"]["name"])
            data["charge_state"].append(gs_e["charge_state"])
            data["pc_from"].append(gs_e["pc_from"].split("/")[-1])
            print(gs_e["pc_from"].split("/")[-1].split("-")[-2])
            data["prototype"].append(gs_e["pc_from"].split("/")[-1].split("-")[-2])
            data["up_cdft_occ"].append(e["cdft_occ"]["up"])
            data["dn_cdft_occ"].append(e["cdft_occ"]["dn"])
            data["gs_taskid"].append(gs_taskid)
            data["task_label"].append(e["task_label"])
            data["energy"].append(e["output"]["energy"])

        for gs in self.cdft_db.collection.distinct("source_entry"):
            gs_taskid = int(gs.split("/")[-1])
            gs_e = self.calc_db.collection.find_one({"task_id": gs_taskid})
            data["task_id"].append(gs_taskid)
            data["defect_name"].append(gs_e["defect_entry"]["name"])
            data["charge_state"].append(gs_e["charge_state"])
            data["pc_from"].append(gs_e["pc_from"].split("/")[-1])
            data["prototype"].append(gs_e["pc_from"].split("/")[-1].split("-")[-2])
            data["up_cdft_occ"].append("ground-state")
            data["dn_cdft_occ"].append("ground-state")
            data["gs_taskid"].append(gs_taskid)
            data["task_label"].append("CDFT-A-HSE_scf")
            gs_energy = gs_e["output"]["energy"]
            data["energy"].append(gs_energy)

        print(data)
        df = pd.DataFrame(data)
        self.foundation_df = df
        IOTools(cwd=save_xlsx_path, pandas_df=df).to_excel("hse_screened_qubit_cdft_2021-11-18")
        return df

    def get_zpl_data(self, excel_file=None):
        df = self.foundation_df.copy() if not excel_file else IOTools(cwd=input_path, excel_file=excel_file).read_excel()
        data = {"gs_taskid":[], "prototype": [], "host": [], "defect":[], "charge": [], "AB": [],
        "BC": [], "CD": [], "DA": [], "ZPL": []}

        for src_taskid in df["gs_taskid"].unique():
            print(src_taskid)
            entry = df.loc[df["gs_taskid"]==src_taskid]
            up_occ = entry["up_cdft_occ"].unique().tolist()
            up_occ.remove('ground-state')
            dn_occ = entry["dn_cdft_occ"].unique().tolist()
            dn_occ.remove("ground-state")

            for up_cdft_occ, dn_cdft_occ in zip(up_occ, dn_occ):
                try:
                    A = entry.loc[(entry["task_label"] == "CDFT-A-HSE_scf"), "energy"].iloc[0]
                except Exception:
                    A = None
                try:
                    B = entry.loc[(entry["task_label"] == "CDFT-B-HSE_scf") & (entry["up_cdft_occ"] == up_cdft_occ) &
                        (entry["dn_cdft_occ"] == dn_cdft_occ), "energy"].iloc[0]
                except Exception:
                    B = None
                try:
                    C = entry.loc[(entry["task_label"] == "CDFT-C-HSE_relax") & (entry["up_cdft_occ"] == up_cdft_occ) &
                                  (entry["dn_cdft_occ"] == dn_cdft_occ), "energy"].iloc[0]
                except Exception:
                    C = None
                try:
                    D = entry.loc[(entry["task_label"] == "CDFT-D-HSE_scf") & (entry["up_cdft_occ"] == up_cdft_occ) &
                                  (entry["dn_cdft_occ"] == dn_cdft_occ), "energy"].iloc[0]
                except Exception:
                    D = None

                data["AB"].append(B-A if A and B else None)
                data["BC"].append(B-C if B and C else None)
                data["CD"].append(C-D if C and D else None)
                data["DA"].append(D-A if D and A else None)
                data["ZPL"].append(C-A if C and A else None)
                data["gs_taskid"].append(src_taskid)
                data["prototype"].append(entry["prototype"].iloc[0])
                data["host"].append(entry["pc_from"].iloc[0].split("-")[0])
                data["defect"].append(entry["defect_name"].iloc[0])
                data["charge"].append(entry["charge_state"].iloc[0])
                print(data)
        result_df = pd.DataFrame(data)
        result_df["ZPL_wavelength"] = result_df.agg(lambda x: 12400/x["ZPL"]/10 if x["ZPL"] else None, axis=1)
        result_df = result_df.round(3)
        result_df["valid"] = result_df.agg(lambda x: True if x["AB"]>=0 and x["BC"]>=0 and x["CD"]>=0 and x["DA"]>=0
        else False,
                                           axis=1)
        result_df.fillna("None", inplace=True)
        result_df.sort_values(["valid", "ZPL"], inplace=True, ascending=False)
        return result_df
def main():
    # a = Defect(defect_xlsx="defect_2021-11-18_2022-01-10")
    # a = Defect()
    # a.get_defect_df_v2_hse()
    # back_process = BackProcess(a.defect_df)
    # back_process.add_band_edges_and_defects()
    # back_process.add_level_category()
    # back_process.add_transition_wavevlength()
    # back_process.add_cdft_occupations()
    # back_process.add_hse_fworker()
    # back_process.df_to_excel(excel_name="hse_screened_qubit_2021-11-18")

    a = Defect(defect_xlsx="defect_2021-11-18")
    screened = a.get_screened_defect_df()
    IOTools(cwd=save_xlsx_path, pandas_df=screened).to_excel("triplet_2021-11-18.xlsx")

    # a.statistics(screened_df)
    # a.get_host_df()
    # a.get_host_gp()
    # a.get_defect_df()
    # a.get_defect_df_gp(screened, "triplet_max-tran_gp_2021-11-18")

    # a = Defect(defect_xlsx="hse_screened_qubit_2021-11-18_2022-01-10")
    # Defect.bar_plot_defect_levels_v2(a.defect_df)
#
if __name__ == '__main__':
    df = main()
    # a = Defect(defect_xlsx="hse_screened_qubit_2021-11-18")
    # show = a.defect_df.loc[(a.defect_df["up_tran_en"]==0) & (a.defect_df["dn_tran_en"]!=0)]
    # # print(a.defect_df.dn_tran_lumo_homo_deg.unique())
    # show = show.loc[:, ["task_id", "up_tran_lumo_homo_band_index", "dn_tran_lumo_homo_band_index",
    #                            "up_tran_lumo_homo_deg",
    #                            "dn_tran_lumo_homo_deg"]]
    # # show = show.loc[(show["dn_tran_lumo_homo_deg"] == (1, 2)) | (show["dn_tran_lumo_homo_deg"] == (3, 2))]

    # IOTools(cwd=save_xlsx_path, pandas_df=df).to_excel("test")
    # a, b, c = Tools.extract_defect_levels_v2_hse(1911)

    # cdft = CDFT()
    # cdft.update_entry_with_occ()
    # cdft.get_data_sheet()
    # df1 = cdft.get_zpl_data()
    # IOTools(cwd=save_xlsx_path, pandas_df=df1).to_excel("hse_screened_qubit_cdft_zpl_2021-11-18")
