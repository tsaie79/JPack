import os, subprocess, datetime

# from analysis.deep_defect import Potential
from qubitPack.defect_formation_energy_correction.defect_formation_2d_corr import FormationEnergy2D
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

from ase.units import _hplanck, _eps0, _c
from VaspBandUnfolding.vasp_constant import *
DEBYETOSI = 3.335640952e-30

SCAN2dMat = get_db("Scan2dMat", "calc_data",  user="Jeng_ro", password="qimin", port=12347)
SCAN2dDefect = get_db("Scan2dDefect", "calc_data",  user="Jeng_ro", password="qimin", port=12347)
SCAN2dIR = get_db("Scan2dDefect", "ir_data", port=12347)

HSEQubitDefect = get_db("HSE_triplets_from_Scan2dDefect", "calc_data-pbe_pc", port=12347, user="Jeng")
HSEQubitIR = get_db("HSE_triplets_from_Scan2dDefect", "ir_data-pbe_pc", port=12347)
HSEQubitCDFT = get_db("HSE_triplets_from_Scan2dDefect", "cdft-pbe_pc", port=12347, user="Jeng_ro")

C2DB_IR = get_db("C2DB_IR", "calc_data", port=12345, user="Jeng_ro")
C2DB_IR_vacancy_HSE = get_db("C2DB_IR_vacancy_HSE", "calc_data", port=12347, user="Jeng_ro")

DB1_PATH = "/home/qimin/sdb_tsai/site-packages/JPack_independent/projects/defectDB"

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
    def extract_defect_levels_v2(cls, defect_taskid, localisation=0.05, edge_tol=(0.25, 0.25), selected_bands=None):
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
                "all",
                None,
                None,  #(host_db, host_taskid, 0, vbm_dx, cbm_dx),
                localisation,  #0.2
                locpot_c2db=None,  #(c2db, c2db_uid, 0)
                is_vacuum_aligment_on_plot=True,
                edge_tol=edge_tol,
                ir_db=None,
                ir_entry_filter={"pc_from_id": pc_from_id, "defect_name": defect_name, "charge_state": charge_state},
                selected_bands=selected_bands,
            )
            tot, proj, d_df, levels, defect_levels = state
            level_info = d_df.to_dict("records")[0]
        except Exception as er:
            tot = None
            level_info = {}
            levels = {}
            defect_levels = {}
            proj = None
        return tot, level_info, levels, defect_levels, proj

    @classmethod
    def extract_defect_levels_v2_hse(cls, defect_taskid, localisation=0.2, edge_tol=(0.5, 0.5), selected_bands=None,
                                     fig_on="all"):
        from qubitPack.qc_searching.analysis.main import get_defect_state_v3
        from qubitPack.tool_box import get_db

        from pymatgen import Structure
        import os


        defect_db = HSEQubitDefect
        ir_db = HSEQubitIR

        defect = defect_db.collection.find_one({"task_id": defect_taskid})
        level_info, levels, defect_levels = None, None, None
        try:
            state = get_defect_state_v3(
                defect_db,
                {"task_id": defect_taskid},
                -10, 10,
                None,
                fig_on,
                None,
                None,  #(host_db, host_taskid, 0, vbm_dx, cbm_dx),
                localisation,  #0.2
                locpot_c2db=None,  #(c2db, c2db_uid, 0)
                is_vacuum_aligment_on_plot=True,
                edge_tol= edge_tol, # defect state will be picked only if it's above vbm by
                # 0.025 eV
            # and below
                # cbm by 0.025 eV
                ir_db=ir_db, #ir_db,
                ir_entry_filter={"prev_fw_taskid": defect_taskid},
                selected_bands=selected_bands,
            )
            tot, proj, d_df, levels, defect_levels = state
            level_info = d_df.to_dict("records")[0]
        except Exception as er:
            print(er)
            tot = None
            level_info = {}
            levels = {}
            defect_levels = {}
            proj = None
        return tot, level_info, levels, defect_levels, proj

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

    def  get_defect_df_v2(self):
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
            info.update({"localisation_threshold": localisation})
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

    def get_defect_df_v2_hse(self, condition=None):
        data = []
        col = HSEQubitDefect.collection
        condition = {"task_label": "HSE_scf", "nupdown_set": 2} if condition is None else condition
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

            localisation = None
            if e["task_id"] in [575]:
                localisation = 0.35
            elif e["task_id"] in [1088, 2590, 585, 545, 571, 2563, 603, 644, ]:
                localisation = 0.1
            elif e["task_id"] in [605, 2569]:
                localisation = 0.05
            else:
                localisation = 0.2
            d_df, levels, in_gpa_levels, proj_df = Tools.extract_defect_levels_v2_hse(e["task_id"],
                                                                                 localisation=localisation)

            info.update(d_df)
            info.update(levels)
            info.update(in_gpa_levels)
            data.append(info)

        self.defect_df = pd.DataFrame(data)
        self.defect_df.fillna("None", inplace=True)
        self.defect_df.replace({"up_tran_en": "None"}, 0, inplace=True)
        self.defect_df.replace({"dn_tran_en": "None"}, 0, inplace=True)

        print(self.defect_df)
        # IOTools(cwd=save_xlsx_path, pandas_df=self.defect_df).to_excel("defects_36_groups")

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
            print(f"set_prototype: {set_prototype}")
            conditions = (defect_df_copy["prototype"] == set_prototype)
            defect_df_copy = defect_df_copy.loc[conditions, :]
            if defect_df_copy.empty:
                print("defect_df_copy is empty!")
                continue
            host_taskids = defect_df_copy.loc[:, "host_taskid"].unique()
            uids = defect_df_copy.loc[:, "uid"].unique()
            print(f"host_taskids: {host_taskids}, uids: {uids}")
            for host_taskid, uid in zip(host_taskids, uids):
                plot_df = defect_df_copy.loc[(defect_df_copy["host_taskid"] == host_taskid) &
                                             (defect_df_copy["level_cbm"] != "None") &
                                             (defect_df_copy["level_vbm"] != "None"), :]
                if plot_df.empty:
                    continue
                plot_df = plot_df.replace("None", np.nan)

                defects = []
                taskid_labels = []
                mags =[]
                defect_name_labels = []
                charge_labels = []
                print(f"task_id :{plot_df['task_id'].tolist()}")
                for idx, taskid in enumerate(plot_df["task_id"].tolist()):
                    print("=="*20, "\n", f"taskid: {taskid}, host_taskid: {host_taskid}, idx: {idx}")
                    defect = plot_df.loc[(plot_df["task_id"] == taskid), :]
                    if len(defect) == 1:
                        idx = 0
                    defect = defect.iloc[idx, :]
                    mag = defect["mag"]

                    # append all entries as dict in defect to defects list
                    defects.append(defect.to_dict())

                    taskid_labels.append(defect["task_id"])
                    mags.append(defect["mag"])
                    defect_name_labels.append(defect["defect_name"])
                    charge_labels.append(defect["charge"])

                vbm_lim, cbm_lim = plot_df["level_vbm"].median(), plot_df["level_cbm"].median()

                fig, ax = plt.subplots(figsize=(12,11), dpi=300)
                fig_y_height = 5, 5
                ax.set_ylim(vbm_lim-fig_y_height[0], cbm_lim+fig_y_height[1])

                print("len of defects: {}".format(len(defects)))
                for defect, x in zip(defects, np.arange(0, len(defects)*5, 5)):
                    # readout all the data from defect dict
                    vbm = defect["level_vbm"]
                    cbm = defect["level_cbm"]
                    level_edge_ir = defect["level_edge_ir"],
                    up = defect["up_in_gap_level"]
                    dn = defect["dn_in_gap_level"]
                    up_deg = defect["up_in_gap_deg"]
                    dn_deg = defect["dn_in_gap_deg"]
                    up_ir = defect["up_in_gap_ir"]
                    dn_ir = defect["dn_in_gap_ir"]
                    up_occ = defect["up_in_gap_occ"]
                    dn_occ = defect["dn_in_gap_occ"]
                    up_id = defect["up_in_gap_band_id"]
                    dn_id = defect["dn_in_gap_band_id"]
                    up_band = defect["up_in_gap_band_index"]
                    dn_band = defect["dn_in_gap_band_index"]
                    vbm_max_el = defect["vbm_max_el"]
                    vbm_max_proj_orbital = defect["vbm_max_proj_orbital"]
                    cbm_max_el = defect["cbm_max_el"]
                    cbm_max_proj_orbital = defect["cbm_max_proj_orbital"]
                    up_tran_bottom = defect["up_tran_bottom"] + defect["level_vbm"]
                    up_tran_en = defect["up_tran_en"]
                    dn_tran_bottom = defect["dn_tran_bottom"] + defect["level_vbm"]
                    dn_tran_en = defect["dn_tran_en"]

                    transition_from = defect.get("transition_from", None)
                    #print transition_from
                    print(f"transition_from: {transition_from}")

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
                        print(f"level: {level}, occ: {occ}, deg: {deg}, id: {id}")
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
                    ax.text(x, cbm+3, "{}-{}-{}|{}".format(cbm_max_el, cbm_max_proj_orbital, level_edge_ir[0][1], cbm),
                            bbox=dict(
                        facecolor='white', edgecolor='black'), size=12)
                    if up_tran_en and (up_tran_bottom or up_tran_bottom != "None"):
                        if transition_from == "up" or not transition_from:
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

                    ax.text(x, vbm_lim-fig_y_height[0]+4, "{}-{}-{}|{}".format(vbm_max_el, vbm_max_proj_orbital,
                                                                            level_edge_ir[0][0], vbm),
                            bbox=dict(facecolor='white', edgecolor='black'), size=12)
                    if dn_tran_en and dn_tran_bottom:
                        if transition_from == "dn" or not transition_from:
                            ax.arrow(x=x+1.25, y=dn_tran_bottom, dx=0, dy=dn_tran_en, width=0.02,
                                     length_includes_head=True, facecolor="gray")

                    # if dn_tran_en and dn_tran_bottom and transition_from == "dn":
                    #     ax.arrow(x=x+1.25, y=dn_tran_bottom, dx=0, dy=dn_tran_en, width=0.02,
                    #              length_includes_head=True, facecolor="gray")

                ax.yaxis.set_minor_locator(AutoMinorLocator(5))
                for tick in ax.get_yticklabels():
                    tick.set_fontsize(15)

                ax.tick_params(axis="x", bottom=False, labelbottom=True)
                ax.tick_params(axis="y", right=False, which="both")

                ax.set_xticks(np.arange(0, len(defects)*5, 5))
                print(f"taskid_label: {taskid_labels}, mags: {mags}")
                ax.set_xticklabels(["{}\nmag:{}".format(taskid_labels[i], mags[i]) for i in range(len(mags))],
                                   fontdict={"fontsize": 10}, rotation=0)

                for x, name, q in zip(np.arange(0, len(defects)*5, 5), defect_name_labels, charge_labels):
                    ax.text(x-2, cbm_lim+fig_y_height[1]+0.1, "{}-{}".format(q, name), rotation=0, size=10)

                ax.set_ylabel("Energy relative to vacuum (eV)", fontsize=15)
                ax.text(-0.1, 1.05, "{}-{}".format(uid, host_taskid), size=30, transform=ax.transAxes)

                save_plt_path = "/Users/jeng-yuantsai/Research/project/Scan2dDefect/latest_results/plt/2021-11-18" \
                                "/defect/table_1"
                os.makedirs(os.path.join(save_plt_path, set_prototype), exist_ok=True)
                fig.savefig(
                    os.path.join(
                        save_plt_path,
                        set_prototype,
                        "host_taskid-{}.png".format(host_taskid)
                    )
                )
                # plt.show()
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

    def add_LUDL_HODL_info(self):
        def dn_fun(x):
            gap_tran_bottom_to_vbm = x["dn_tran_bottom"]
            if x["level_gap"] != "None" and x["dn_tran_top"] != "None":
                gap_tran_top_to_cbm = x["level_gap"] - x["dn_tran_top"]
            else:
                gap_tran_top_to_cbm = "None"
            return gap_tran_bottom_to_vbm, gap_tran_top_to_cbm

        def up_fun(x):
            gap_tran_bottom_to_vbm = x["up_tran_bottom"]
            if x["level_gap"] != "None" and x["up_tran_top"] != "None":
                gap_tran_top_to_cbm = x["level_gap"] - x["up_tran_top"]
            else:
                gap_tran_top_to_cbm = "None"
            return gap_tran_bottom_to_vbm, gap_tran_top_to_cbm

        df = self.input_df
        df["up_tran_bottom_to_vbm"] = df.apply(lambda x: up_fun(x)[0], axis=1)
        df["up_tran_top_to_cbm"] = df.apply(lambda x: up_fun(x)[1], axis=1)
        df["dn_tran_bottom_to_vbm"] = df.apply(lambda x: dn_fun(x)[0], axis=1)
        df["dn_tran_top_to_cbm"] = df.apply(lambda x: dn_fun(x)[1], axis=1)

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

    def add_singlet_triplet_en_diff(self):
        from JPack_independent.projects.defectDB.analysis.analysis_api import DefectAnalysis
        da = DefectAnalysis(self.input_df)
        da.get_singlet_triplet_en_diff()
        self.input_df = self.input_df.merge(
            da.singlet_df, on=["C2DB_uid", "defect_name", "charge", "task_id"], how="left")

    def add_defect_symmetry(self):
        from JPack_independent.projects.defectDB.analysis.analysis_api import DefectAnalysis
        da = DefectAnalysis(self.input_df)
        da.get_defect_symmetry(self.input_df)
        self.input_df = self.input_df.merge(da.defect_symmetry_df, on=["task_id"], how="left")

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
                    cdft_occ = np.nan
                cdft_occs[spin] = cdft_occ
            return cdft_occs["1"], cdft_occs["-1"]

        defect_df = self.input_df
        defect_df["up_tran_cdft_occ"] = defect_df.apply(lambda x: fun(x)[0], axis=1)
        defect_df["dn_tran_cdft_occ"] = defect_df.apply(lambda x: fun(x)[1], axis=1)
        defect_df.fillna("None", inplace=True)

    def add_tdm_input(self):
        def fun(x):
            mode = ["1*0 1*1", "1*0.5 1*0.5 1*1", "1*0 1*0.5 1*0.5", "1*0.5 1*0.5 1*0.5 1*0.5",
                    "1*0.5 1*0.5 1*0.333333 1*0.333333 1*0.333333"]
            up = " ".join(x["up_cdft_occ"].split(" ")[1:-1])
            dn = " ".join(x["dn_cdft_occ"].split(" ")[1:-1])
            up_first = int(x["up_cdft_occ"].split("*")[0])
            dn_first = int(x["dn_cdft_occ"].split("*")[0])
            T = {"up": [], "dn": []}
            for cdft_occ, first, T_name in zip((up, dn), (up_first, dn_first), ("up", "dn")):
                if cdft_occ == mode[0]:
                    T[T_name].append((first+1, first+2))
                elif cdft_occ == mode[1]:
                    T[T_name].append((first+1, first+3))
                    T[T_name].append((first+2, first+3))
                elif cdft_occ == mode[2]:
                    T[T_name].append((first+1, first+2))
                    T[T_name].append((first+1, first+3))
                elif cdft_occ == mode[3]:
                    T[T_name].append((first+1, first+3))
                    T[T_name].append((first+1, first+4))
                    T[T_name].append((first+2, first+3))
                    T[T_name].append((first+2, first+4))
                elif cdft_occ == mode[4]:
                    T[T_name].append((first+1, first+3))
                    T[T_name].append((first+1, first+4))
                    T[T_name].append((first+1, first+5))
                    T[T_name].append((first+2, first+3))
                    T[T_name].append((first+2, first+4))
                    T[T_name].append((first+2, first+5))
                else:
                    T[T_name] = None
            T_up, T_dn = None, None
            if T["up"]:
                T_up = tuple(T["up"])
            if T["dn"]:
                T_dn = tuple(T["dn"])
            return T_up, T_dn
        self.input_df["up_TDM_input"] = self.input_df.loc[self.input_df["task_label"] == "CDFT-B-HSE_scf"].apply(
            lambda x: fun(x)[0], axis=1)
        self.input_df["dn_TDM_input"] = self.input_df.loc[self.input_df["task_label"] == "CDFT-B-HSE_scf"].apply(
            lambda x: fun(x)[1], axis=1)

    def add_IPR(self):
        from JPack_independent.projects.defectDB.analysis.analysis_api import DefectAnalysis
        da = DefectAnalysis(self.input_df)
        da.get_IPR(self.input_df)
        self.input_df = self.input_df.merge(da.ipr_df, on=["task_id"], how="left")

    def add_IPR_average(self):
        def fun(x):
            dn_in_gap_ipr = x["dn_in_gap_ipr"]
            up_in_gap_ipr = x["up_in_gap_ipr"]
            dn_in_gap_ipr_avg, up_in_gap_ipr_avg = 0, 0
            if dn_in_gap_ipr != ():
                dn_in_gap_ipr_avg = np.mean(dn_in_gap_ipr)
            if up_in_gap_ipr != ():
                up_in_gap_ipr_avg = np.mean(up_in_gap_ipr)
            return up_in_gap_ipr_avg, dn_in_gap_ipr_avg
        self.input_df["up_in_gap_ipr_avg"], self.input_df["dn_in_gap_ipr_avg"] = zip(*self.input_df.apply(fun, axis=1))


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


    def get_data_sheet(self, filter, read_c2db_uid_key=False):
        es = self.cdft_db.collection.find(filter)
        data = {"prototype": [], "pc_from": [], "gs_taskid": [], "task_id": [], "defect_name": [],
                "charge_state":[], "up_cdft_occ":[],"dn_cdft_occ": [], "task_label":[], "energy": [],
                "up_TDM_sq_X": [], "up_TDM_sq_Y": [], "up_TDM_sq_Z": [],
                "dn_TDM_sq_X": [], "dn_TDM_sq_Y": [], "dn_TDM_sq_Z": []
        }
        for e in es:
            gs_taskid = int(e["source_entry"].split("/")[-1])
            gs_e = self.calc_db.collection.find_one({"task_id": gs_taskid})
            data["task_id"].append(e["task_id"])
            data["defect_name"].append(gs_e["defect_entry"]["name"])
            data["charge_state"].append(gs_e["charge_state"])
            if read_c2db_uid_key:
                try:
                    data["prototype"].append(gs_e["pc_from"].split("/")[-1].split("-")[-2])
                    data["pc_from"].append(gs_e["pc_from"].split("/")[-1])
                except Exception as err:
                    print(err, gs_taskid)
                    data["prototype"].append(gs_e["c2db_uid"].split("-")[-2])
                    data["pc_from"].append(gs_e["c2db_uid"])
            else:
                data["pc_from"].append(gs_e["pc_from"].split("/")[-1])
                data["prototype"].append(gs_e["pc_from"].split("/")[-1].split("-")[-2])
            data["up_cdft_occ"].append(e["cdft_occ"]["up"])
            data["dn_cdft_occ"].append(e["cdft_occ"]["dn"])
            data["gs_taskid"].append(gs_taskid)
            data["task_label"].append(e["task_label"])
            data["energy"].append(e["output"]["energy"])
            for idx, component in enumerate(["X", "Y", "Z"]):
                try:
                    data["up_TDM_sq_{}".format(component)].append(e["TDM_transition"]["up_TDM_squared"][idx])
                except Exception:
                    data["up_TDM_sq_{}".format(component)].append(None)
                try:
                    data["dn_TDM_sq_{}".format(component)].append(e["TDM_transition"]["dn_TDM_squared"][idx])
                except Exception:
                    data["dn_TDM_sq_{}".format(component)].append(None)

        # get unique gs_taskid from es
        gs_taskids = list(set(data["gs_taskid"]))
        for gs_taskid in gs_taskids:
            gs_e = self.calc_db.collection.find_one({"task_id": gs_taskid})
            data["task_id"].append(gs_taskid)
            data["defect_name"].append(gs_e["defect_entry"]["name"])
            data["charge_state"].append(gs_e["charge_state"])
            if read_c2db_uid_key:
                try:
                    data["prototype"].append(gs_e["pc_from"].split("/")[-1].split("-")[-2])
                    data["pc_from"].append(gs_e["pc_from"].split("/")[-1])
                except Exception as err:
                    print(err, gs_taskid)
                    data["prototype"].append(gs_e["c2db_uid"].split("-")[-2])
                    data["pc_from"].append(gs_e["c2db_uid"])
            else:
                data["pc_from"].append(gs_e["pc_from"].split("/")[-1])
                data["prototype"].append(gs_e["pc_from"].split("/")[-1].split("-")[-2])
            data["up_cdft_occ"].append("ground-state")
            data["dn_cdft_occ"].append("ground-state")
            data["gs_taskid"].append(gs_taskid)
            data["task_label"].append("CDFT-A-HSE_scf")
            gs_energy = gs_e["output"]["energy"]
            data["energy"].append(gs_energy)
            for idx, component in enumerate(["X", "Y", "Z"]):
                data["up_TDM_sq_{}".format(component)].append(None)
                data["dn_TDM_sq_{}".format(component)].append(None)

        df = pd.DataFrame(data)
        BackProcess(df).add_tdm_input()
        self.foundation_df = df
        print(f"\nCDFT_data:\n{df}")
        # IOTools(cwd=os.path.join(DB1_PATH, save_xlsx_path), pandas_df=df).to_excel("hse_screened_qubit_cdft_2021-11-18")
        return df

    def get_zpl_data(self, excel_file=None):
        df = self.foundation_df.copy() if not excel_file else IOTools(cwd=input_path, excel_file=excel_file).read_excel()
        data = {"gs_taskid":[], "prototype": [], "host": [], "defect":[], "charge": [], "AB": [],
                "BC": [], "CD": [], "DA": [], "ZPL": [],
                "up_TDM_sq_X": [], "up_TDM_sq_Y": [], "up_TDM_sq_Z": [],
                "dn_TDM_sq_X": [], "dn_TDM_sq_Y": [], "dn_TDM_sq_Z": []
        }

        for src_taskid in df["gs_taskid"].unique():
            base_entry = df.loc[df["gs_taskid"]==src_taskid]
            up_occ = base_entry["up_cdft_occ"].unique().tolist()
            up_occ.remove('ground-state')
            dn_occ = base_entry["dn_cdft_occ"].unique().tolist()
            dn_occ.remove("ground-state")

            for up_cdft_occ, dn_cdft_occ in zip(up_occ, dn_occ):
                entry = base_entry.loc[(base_entry["up_cdft_occ"] == up_cdft_occ) &
                                       (base_entry["dn_cdft_occ"] == dn_cdft_occ), :].copy()
                try:
                    A = base_entry.loc[(base_entry["task_label"] == "CDFT-A-HSE_scf"), "energy"].iloc[0]
                except Exception:
                    A = None
                try:
                    B = entry.loc[(entry["task_label"] == "CDFT-B-HSE_scf"), "energy"].iloc[0]
                except Exception:
                    B = None
                try:
                    C = entry.loc[(entry["task_label"] == "CDFT-C-HSE_relax"), "energy"].iloc[0]
                except Exception:
                    C = None
                try:
                    D = entry.loc[(entry["task_label"] == "CDFT-D-HSE_scf"), "energy"].iloc[0]
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

                for idx, component in enumerate(["X", "Y", "Z"]):
                    data["up_TDM_sq_{}".format(component)].append(
                        entry.loc[(entry["task_label"] == "CDFT-B-HSE_scf"), "up_TDM_sq_{}".format(component)].iloc[0])
                    data["dn_TDM_sq_{}".format(component)].append(
                        entry.loc[(entry["task_label"] == "CDFT-B-HSE_scf"), "dn_TDM_sq_{}".format(component)].iloc[0])

        def calculate_TDM_rate(TDM_sq, E, refractive_index=1):
            if TDM_sq is None or E is None:
                return None
            else:
                TDM_rate = refractive_index * E**3 * TDM_sq / 3 / PI / _eps0 / _c**3 / (_hplanck/(2*PI))**4
                return round(TDM_rate*1e-6, 3)

        result_df = pd.DataFrame(data)

        for spin in ["up", "dn"]:
            for component in ["X", "Y", "Z"]:
                result_df["{}_TDM_rate_{}".format(spin, component)] = result_df.apply(
                    lambda x: calculate_TDM_rate(x["{}_TDM_sq_{}".format(spin, component)], x["ZPL"]*EVTOJ)
                    if x["{}_TDM_sq_{}".format(spin, component)] is not None and x["ZPL"] is not None else None, axis=1)

                # result_df["{}_TDM_rate_{}".format(spin, component)] = \
                #     calculate_TDM_rate(result_df["{}_TDM_sq_{}".format(spin, component)], result_df["ZPL"]*EVTOJ)

        result_df["total_TDM_rate_up"] = sum([result_df["{}_TDM_rate_{}".format(spin, component)]
                                              for spin in ["up"] for component in ["X", "Y", "Z"]])
        result_df["total_TDM_rate_dn"] = sum([result_df["{}_TDM_rate_{}".format(spin, component)]
                                              for spin in ["dn"] for component in ["X", "Y", "Z"]])

        # iterate over all the rows
        for idx, row in result_df.iterrows():
            # get the up TDM rate
            up_TDM_rate_X = row["up_TDM_rate_X"]
            up_TDM_rate_Y = row["up_TDM_rate_Y"]
            up_TDM_rate_Z = row["up_TDM_rate_Z"]
            print(f"\n==gs_taskid: {row['gs_taskid']}==")
            print(f"\nup_TDM_rate_X: {up_TDM_rate_X}\nup_TDM_rate_Y: {up_TDM_rate_Y}\nup_TDM_rate_Z: {up_TDM_rate_Z}")
            # print("up_TDM_rate_X: {}".format(up_TDM_rate_X), "up_TDM_rate_Y: {}".format(up_TDM_rate_Y),
            #       "up_TDM_rate_Z: {}".format(up_TDM_rate_Z))
            if up_TDM_rate_X != 0 and up_TDM_rate_Y == 0 and up_TDM_rate_Z == 0:
                result_df.loc[idx, "up_polarization"] = "X"
                result_df.loc[idx, "up_allowed"] = True
                result_df.loc[idx, "up_transition"] =  True
            elif up_TDM_rate_X == 0 and up_TDM_rate_Y != 0 and up_TDM_rate_Z == 0:
                result_df.loc[idx, "up_polarization"] = "Y"
                result_df.loc[idx, "up_allowed"] = True
                result_df.loc[idx, "up_transition"] =  True
            elif up_TDM_rate_X == 0 and up_TDM_rate_Y == 0 and up_TDM_rate_Z != 0:
                result_df.loc[idx, "up_polarization"] = "Z"
                result_df.loc[idx, "up_allowed"] = True
                result_df.loc[idx, "up_transition"] =  True
            elif up_TDM_rate_X != 0 and up_TDM_rate_Y != 0 and up_TDM_rate_Z == 0:
                result_df.loc[idx, "up_polarization"] = "XY"
                result_df.loc[idx, "up_allowed"] = True
                result_df.loc[idx, "up_transition"] =  True
            elif up_TDM_rate_X != 0 and up_TDM_rate_Y == 0 and up_TDM_rate_Z != 0:
                result_df.loc[idx, "up_polarization"] = "XZ"
                result_df.loc[idx, "up_allowed"] = True
                result_df.loc[idx, "up_transition"] =  True
            elif up_TDM_rate_X == 0 and up_TDM_rate_Y != 0 and up_TDM_rate_Z != 0:
                result_df.loc[idx, "up_polarization"] = "YZ"
                result_df.loc[idx, "up_allowed"] = True
                result_df.loc[idx, "up_transition"] =  True
            elif up_TDM_rate_X != 0 and up_TDM_rate_Y != 0 and up_TDM_rate_Z != 0:
                if pd.isna(up_TDM_rate_X) and pd.isna(up_TDM_rate_Y) and pd.isna(up_TDM_rate_Z):
                    result_df.loc[idx, "up_polarization"] = np.nan
                    result_df.loc[idx, "up_allowed"] = False
                    result_df.loc[idx, "up_transition"] =  False
                elif abs(up_TDM_rate_X - up_TDM_rate_Y) < 1:
                    if up_TDM_rate_X - up_TDM_rate_Z > 0:
                        result_df.loc[idx, "up_polarization"] = "XY"
                        result_df.loc[idx, "up_allowed"] = True
                        result_df.loc[idx, "up_transition"] =  True
                    else:
                        result_df.loc[idx, "up_polarization"] = "Z"
                        result_df.loc[idx, "up_allowed"] = True
                        result_df.loc[idx, "up_transition"] =  True
                else:
                    result_df.loc[idx, "up_polarization"] = "XYZ"
                    result_df.loc[idx, "up_allowed"] = True
                    result_df.loc[idx, "up_transition"] =  True
            else:
                result_df.loc[idx, "up_polarization"] = np.nan
                result_df.loc[idx, "up_allowed"] = False
                result_df.loc[idx, "up_transition"] =  True
            print(f"\nup_polarization: {result_df.loc[idx, 'up_polarization']}")
            # print("up_polarization: {}, {}".format(result_df.loc[idx, "up_polarization"], type(result_df.loc[idx, "up_polarization"])))
            # get the dn TDM rate
            dn_TDM_rate_X = row["dn_TDM_rate_X"]
            dn_TDM_rate_Y = row["dn_TDM_rate_Y"]
            dn_TDM_rate_Z = row["dn_TDM_rate_Z"]
            print(f"\ndn_TDM_rate_X: {dn_TDM_rate_X}\ndn_TDM_rate_Y: {dn_TDM_rate_Y}\ndn_TDM_rate_Z: {dn_TDM_rate_Z}")
            # print("dn_TDM_rate_X: {}".format(dn_TDM_rate_X), "dn_TDM_rate_Y: {}".format(dn_TDM_rate_Y),
            #       "dn_TDM_rate_Z: {}".format(dn_TDM_rate_Z))
            if dn_TDM_rate_X != 0 and dn_TDM_rate_Y == 0 and dn_TDM_rate_Z == 0:
                result_df.loc[idx, "dn_polarization"] = "X"
                result_df.loc[idx, "dn_allowed"] = True
                result_df.loc[idx, "dn_transition"] =  True
            elif dn_TDM_rate_X == 0 and dn_TDM_rate_Y != 0 and dn_TDM_rate_Z == 0:
                result_df.loc[idx, "dn_polarization"] = "Y"
                result_df.loc[idx, "dn_allowed"] = True
                result_df.loc[idx, "dn_transition"] =  True
            elif dn_TDM_rate_X == 0 and dn_TDM_rate_Y == 0 and dn_TDM_rate_Z != 0:
                result_df.loc[idx, "dn_polarization"] = "Z"
                result_df.loc[idx, "dn_allowed"] = True
                result_df.loc[idx, "dn_transition"] =  True
            elif dn_TDM_rate_X != 0 and dn_TDM_rate_Y != 0 and dn_TDM_rate_Z == 0:
                result_df.loc[idx, "dn_polarization"] = "XY"
                result_df.loc[idx, "dn_allowed"] = True
                result_df.loc[idx, "dn_transition"] =  True
            elif dn_TDM_rate_X != 0 and dn_TDM_rate_Y == 0 and dn_TDM_rate_Z != 0:
                result_df.loc[idx, "dn_polarization"] = "XZ"
                result_df.loc[idx, "dn_allowed"] = True
                result_df.loc[idx, "dn_transition"] =  True
            elif dn_TDM_rate_X == 0 and dn_TDM_rate_Y != 0 and dn_TDM_rate_Z != 0:
                result_df.loc[idx, "dn_polarization"] = "YZ"
                result_df.loc[idx, "dn_allowed"] = True
                result_df.loc[idx, "dn_transition"] =  True
            elif dn_TDM_rate_X != 0 and dn_TDM_rate_Y != 0 and dn_TDM_rate_Z != 0:
                if pd.isna(dn_TDM_rate_X) and pd.isna(dn_TDM_rate_Y) and pd.isna(dn_TDM_rate_Z):
                    result_df.loc[idx, "dn_polarization"] = np.nan
                    result_df.loc[idx, "dn_allowed"] = False
                    result_df.loc[idx, "dn_transition"] =  False
                elif abs(dn_TDM_rate_X - dn_TDM_rate_Y) < 1:
                    if dn_TDM_rate_X - dn_TDM_rate_Z > 0:
                        result_df.loc[idx, "dn_polarization"] = "XY"
                        result_df.loc[idx, "dn_allowed"] = True
                        result_df.loc[idx, "dn_transition"] =  True
                    else:
                        result_df.loc[idx, "dn_polarization"] = "Z"
                        result_df.loc[idx, "dn_allowed"] = True
                        result_df.loc[idx, "dn_transition"] =  True
                else:
                    result_df.loc[idx, "dn_polarization"] = "XYZ"
                    result_df.loc[idx, "dn_allowed"] = True
                    result_df.loc[idx, "dn_transition"] =  True
            else:
                result_df.loc[idx, "dn_polarization"] = np.nan
                result_df.loc[idx, "dn_allowed"] = False
                result_df.loc[idx, "dn_transition"] =  True
            print("\ndn_polarization: {}".format(result_df.loc[idx, "dn_polarization"]))
            # print("dn_polarization: {}, {}".format(result_df.loc[idx, "dn_polarization"], type(result_df.loc[idx, "dn_polarization"])))
            result_df.loc[idx, "allowed"] = result_df.loc[idx, "up_allowed"] or result_df.loc[idx, "dn_allowed"]

            if result_df.loc[idx, "up_transition"] and not result_df.loc[idx, "dn_transition"]:
                result_df.loc[idx, "transition_from"] = "up"
            elif result_df.loc[idx, "dn_transition"] and not result_df.loc[idx, "up_transition"]:
                result_df.loc[idx, "transition_from"] = "dn"
            elif result_df.loc[idx, "up_transition"] and result_df.loc[idx, "dn_transition"]:
                result_df.loc[idx, "transition_from"] = "both"
            else:
                result_df.loc[idx, "transition_from"] = np.nan
        result_df.drop(columns=["up_allowed", "dn_allowed"], inplace=True)

        result_df.drop(["up_TDM_sq_{}".format(component) for component in ["X", "Y", "Z"]], axis=1, inplace=True)
        result_df.drop(["dn_TDM_sq_{}".format(component) for component in ["X", "Y", "Z"]], axis=1, inplace=True)

        result_df["ZPL_wavelength"] = result_df.agg(lambda x: 12400/x["ZPL"]/10 if x["ZPL"] else None, axis=1)
        result_df = result_df.round(3)
        result_df["valid"] = result_df.agg(
            lambda x: True if (x["AB"]>=0 and x["BC"]>=0 and x["CD"]>=0 and x["DA"]>=0) else False, axis=1)
        result_df.sort_values(["valid", "ZPL"], inplace=True, ascending=False)
        result_df.fillna("None", inplace=True)
        return result_df


class COHP:
    # LOBSTER_DB = get_db("HSE_triplets_from_Scan2dDefect", "lobster",  user="Jeng_ro", password="qimin", port=12347)
    LOBSTER_DB = get_db("Scan2dDefect", "lobster",  user="Jeng_ro", password="qimin", port=12347)

    @classmethod
    def find_label_in_cohp(cls, completecohp, nn, defect_type="vacancy"):
        nn_pair = None
        if defect_type == "antisite":
            nn_pair = [{nn[-1], i} for i in nn] # for antisite NN[-1] only
        elif defect_type == "vacancy":
            nn_pair = [{i, j} for i in nn for j in nn] # for vacancy

        label_list = []
        for label, cohp_pair in completecohp.bonds.items():
            print()
            sites_index = tuple(completecohp.structure.sites.index(site) for site in cohp_pair["sites"])
            if set(sites_index) in nn_pair:
                print(label, cohp_pair["length"])
                label_list.append(label)
        return label_list

    @classmethod
    def analysis(cls):
        from pymatgen.electronic_structure.cohp import CompleteCohp
        from pymatgen.electronic_structure.plotter import CohpPlotter
        from pymatgen.electronic_structure.core import Orbital

        def orbital_resolved(label, orbitals):
            #search for the number of the COHP you would like to plot in ICOHPLIST.lobster (the numbers in COHPCAR.lobster are different!)
            label = str(label)
            cp = CohpPlotter()
            #get a nicer plot label
            plotlabel=str(completecohp.bonds[label]['sites'][0].species_string)+'-'+str(completecohp.bonds[label]['sites'][1].species_string+"{}".format(orbitals))

            cp.add_cohp(plotlabel,completecohp.get_orbital_resolved_cohp(label=label, orbitals=orbitals))
            #check which COHP you are plotting

            print("This is a COHP between the following sites: "+str(completecohp.bonds[label]['sites'][0])+' and '+ str(completecohp.bonds[label]['sites'][1]))

            x = cp.get_plot(integrated=False)
            x.ylim([-5, 5])
            x.show()

        def sum_lobster_orbitals(lobster_list, orbitals, plot_label):
            cp = CohpPlotter()
            #"W75-W50+W55+W56 bonds"
            cp.add_cohp(plot_label, completecohp.get_summed_cohp_by_label_and_orbital_list(lobster_list, orbitals))
            x = cp.get_plot(integrated=False)
            x.ylim([-5, 5])
            x.show()

        def summed(completecohp, lobster_list, plot_label, title, is_coop):
            cp = CohpPlotter(are_coops=is_coop)
            # get a nicer plot label
            cp.add_cohp(plot_label, completecohp.get_summed_cohp_by_label_list(label_list=lobster_list, divisor=1))
            x = cp.get_plot(integrated=False, ylim=[-5, 5])
            x.legend(loc="upper right")
            # x.ylim([-5, 5])
            # x.xlim([-1, 1])
            x.title(title, fontsize=20)
            # x.savefig("{}.png".format(title))
            x.show()

        for lobster in cls.LOBSTER_DB.collection.find({"task_id": {"$in": [6162]}}):
            is_coop = True
            if is_coop:
                f = "COOPCAR.lobster.gz"
            else:
                f = "COHPCAR.lobster.gz"
            dir_name = lobster["dir_name"]
            print(dir_name)
            os.chdir(dir_name)
            subprocess.call(["gunzip", "POSCAR.gz"])
            completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename=f,
                                                  structure_file="POSCAR", are_coops=is_coop)
            label_list = cls.find_label_in_cohp(completecohp, lobster["NN"])
            summed(completecohp, label_list, "{}".format(lobster["NN"]), "{}".format(lobster["prev_fw_taskid"]),
                   is_coop=is_coop)
    # @classmethod
    # def scp_COHP_diagram(cls):
    #     for lobster in cls.LOBSTER_DB.collection.find({}):
    #         dir_name = lobster["dir_name"]
    #         title = lobster["prev_fw_taskid"]
    #         print(dir_name)
    #         subprocess.call(["scp", dir_name+"/{}.png".format(title)])

class TransitionLevels:
    @classmethod
    def regular_antisite(cls, db_files, compounds, prototype, Cs, defect_name):
        q0_mag1 = [25,
                   29,
                   45,
                   61,
                   112,
                   119,
                   125,
                   129,
                   133,
                   178,
                   184,
                   197,
                   208,
                   217,
                   225,
                   226,
                   228,
                   235,
                   238,
                   241,
                   256,
                   294,
                   296,
                   307,
                   309,
                   310,
                   317,
                   320,
                   324,
                   339,
                   346,
                   355,
                   360,
                   366,
                   380,
                   386,
                   391,
                   403,
                   422,
                   425,
                   440,
                   449,
                   461,
                   475,
                   485,
                   486,
                   507,
                   511,
                   528,
                   545,
                   552,
                   556,
                   569,
                   574,
                   580,
                   588,
                   617,
                   634,
                   636,
                   641,
                   644,
                   655,
                   657,
                   662,
                   664,
                   666,
                   673,
                   676,
                   681,
                   733,
                   736,
                   744,
                   747,
                   749,
                   756,
                   767,
                   768,
                   772,
                   776,
                   778]
        q1_mag0 = [36,
                   88,
                   103,
                   106,
                   111,
                   126,
                   131,
                   135,
                   180,
                   187,
                   202,
                   212,
                   216,
                   246,
                   255,
                   308,
                   318,
                   328,
                   335,
                   341,
                   342,
                   344,
                   354,
                   358,
                   378,
                   385,
                   390,
                   400,
                   418,
                   421,
                   424,
                   443,
                   452,
                   454,
                   481,
                   484,
                   505,
                   512,
                   513,
                   525,
                   530,
                   534,
                   546,
                   558,
                   571,
                   575,
                   590,
                   774]

        bulk_k = [[1,1,1]]
        vbm_k = [[0, 0, 0]]
        sites = [[33, 31], [51, 49]] if defect_name == "M_on_X" else [[31, 33], [49, 51]]
        charge_spin_states = [[-1, 0, 1], [2, 1, 0]] if defect_name == "M_on_X" else [[1, 0, -1, -2], [2, 3, 1]]
        results = []
        raw_data_df = []
        for compound, c in zip(compounds, Cs):
            ionize_en = FormationEnergy2D(
                bulk_db_files=db_files,
                bulk_entry_filter={
                    # "nsites": {"$ne":48},
                    "chemsys": compound,
                    "input.structure.lattice.c": {"$in": c},
                    "calcs_reversed.input.kpoints.kpoints.0": {"$in": bulk_k},
                    "defect_type": "host",
                    "task_label": "HSE_scf",
                    "c2db_uid": {"$regex": prototype}
                },
                bk_vbm_bg_db_files=db_files,
                bk_vbm_bg_filter={
                    # "nsites": {"$ne":48},
                    "chemsys": compound,
                    "calcs_reversed.input.kpoints.kpoints.0": {"$in": vbm_k},
                    "input.structure.lattice.c": {"$in": c},
                    "defect_type": "vbm",
                    "task_label": "HSE_scf",
                    "c2db_uid": {"$regex": prototype},
                },
                defect_db_files=db_files,
                defect_entry_filter={
                    "task_id": {"$nin": q0_mag1+q1_mag0},
                    "charge_state": {"$in": charge_spin_states[0]},
                    "nupdown_set": {"$in": charge_spin_states[1]},
                    "chemsys": compound,
                    "$or": [{ "formula_pretty": {"$regex": "[a-zA-Z]+{}[a-zA-Z]+{}".format(*sites[0])}},
                            {"formula_pretty": {"$regex": "[a-zA-Z]+{}[a-zA-Z]+{}".format(*sites[1])}}],
                    "calcs_reversed.input.kpoints.kpoints.0": {"$in": bulk_k},
                    "input.structure.lattice.c": {"$in": c},
                    "defect_type": "defect",
                    "task_label": "HSE_scf",
                    "c2db_uid": {"$regex": prototype},
                }
            )

            # replace problematic entry 1168 with 585 from HSEQubitDefect
            ionize_en.defect_entries = [i for i in ionize_en.defect_entries if i["task_id"]!=1168]
            ionize_en.defect_entries.append(HSEQubitDefect.collection.find_one({"task_id": 585}))

            raw_data, tls = ionize_en.get_transition_levels()
            name = f"{compound.split('-')[0]}_on_{compound.split('-')[1]}" if defect_name == "X_on_M" \
                else f"{compound.split('-')[1]}_on_{compound.split('-')[0]}"
            for a in ionize_en.get_ionized_energy(fig_title=f"{name}-{prototype}", tls=tls):
                a.update({"chemsys": compound, "prototype": prototype, "defect_name": defect_name})
                results.append(a)

            raw_data_df.append(raw_data)

        print("$$"*50, "all")
        tls = pd.DataFrame(results)
        tls["IE0"] = tls.apply(lambda x: x["IE0"][0], axis=1)
        tls["tl"] = tls.apply(lambda x: x["IE0"] if np.any(np.array(x["charge"])<0) else x["Bandgap_avg"] - x["IE0"],
                              axis=1)
        # tls["tl"] = tls.apply(lambda x: None if x["IE0"] > x["Bandgap_avg"] else x["tl"], axis=1)
        # groupby chemsys and subtract tl each row from the next row
        tls["window"] = tls["tl"].diff()

        raw_data_df = pd.concat(raw_data_df)
        display(raw_data_df)
        display(tls)
        return tls, raw_data_df

    def get_tls_plot(self, df, tls_df):
            # plot a bar chart containing the defect levels for a mutiple defects
            vbm_df = df.loc[df["defect_type"] == "vbm"]
            #groupby c2db_uid, chemsys, and formula; and get mean of vbm, cbm
            vbm_df = vbm_df.groupby(["c2db_uid", "chemsys", "formula", "prototype"]).agg({"vbm": "mean",
                                                                             "cbm": "mean"}).reset_index()

            tls_df = tls_df.merge(vbm_df, on=["prototype", "chemsys"], how="left")

            fig, ax = plt.subplots(figsize=(12,11), dpi=300)
            fig_y_height = 1, 1.5
            font_size = 15 #7.25
            vbm_lim, cbm_lim = tls_df["vbm"].min(), tls_df["cbm"].max()
            ax.set_ylim(vbm_lim-fig_y_height[0], cbm_lim+fig_y_height[1])

            defect_labels = []
            # plot mutilple bars

            # for (defect, df), x in zip(tls_df.groupby(["prototype", "chemsys", "defect_name"]),
            #                       np.arange(0, len(tls_df.groupby(["prototype", "chemsys", "defect_name"]))*5, 5)):
            for (defect, df), x in zip(tls_df.groupby(["defect_name", "chemsys"]),
                                       np.arange(0, len(tls_df.groupby(["defect_name", "chemsys"]))*5, 5)):
                print(f"defect: {defect}")
                # name = f"{defect[1].split('-')[1]}_on_{defect[1].split('-')[0]}" if defect[2] == "X_on_M" \
                #     else f"{defect[1].split('-')[0]}_on_{defect[1].split('-')[1]}"
                # defect_labels.append(rf'$\mathrm{{{name.split("_")[0]}_{{{name.split("_")[-1]}}}}}$')
                name = f"{defect[1].split('-')[1]}_on_{defect[1].split('-')[0]}" if defect[0] == "X_on_M" \
                    else f"{defect[1].split('-')[0]}_on_{defect[1].split('-')[1]}"
                defect_labels.append(rf'$\mathrm{{{name.split("_")[0]}_{{{name.split("_")[-1]}}} @ '
                                     rf'{"".join(defect[1].split("-"))}}}$')
                # readout all the data from defect dict
                vbm = df["vbm"].mean()
                cbm = df["cbm"].mean()
                # level_edge_ir = defect["level_edge_ir"],
                tl = df["tl"]+df["vbm"]
                charge = df["charge"]

                ax.bar(x, vbm-(vbm_lim-fig_y_height[0]), 4, (vbm_lim-fig_y_height[0]), color="deepskyblue")
                ax.bar(x, (cbm_lim+fig_y_height[1])-cbm, 4, cbm, color="orange")

                dx = 0.35
                for level, q in zip(tl, charge):
                    print(f"charge: {q}, level: {level}")
                    q_labels = []
                    for q_label in q:
                        if q_label < 0:
                            q_labels.append(f"-" if q_label == -1 else f"{-1*q_label}-")
                        elif q_label > 0:
                            q_labels.append(f"+" if q_label == 1 else f"{q_label}+")
                        else:
                            q_labels.append("0")

                    ax.text(x+0.2, level, rf"$\mathrm{{\epsilon({q_labels[1]}/{q_labels[0]})}}$",
                            fontsize=font_size,
                            color="black")
                    color = "black"
                    ax.hlines(level, x-1.5, x-0.1, colors=color)

                # ax.text(x-1.5, cbm+0.4, f"{cbm}",bbox=dict(facecolor='white', edgecolor='black'), size=font_size)
                # ax.text(x-1.5, vbm-0.55, f"{vbm}",bbox=dict(facecolor='white', edgecolor='black'), size=font_size)

            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            for tick in ax.get_yticklabels():
                tick.set_fontsize(15)

            ax.tick_params(axis="x", bottom=False, labelbottom=True)
            ax.tick_params(axis="y", right=False, which="both", length=5, width=1)


            ax.set_xticks(np.arange(0, len(tls_df.groupby(["prototype", "chemsys", "defect_name"]))*5, 5), )
            ax.set_xticklabels(defect_labels, fontdict={"fontsize": 20}, rotation=270)

            # remove minor and major xticks
            ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=True)
            ax.set_ylabel("Energy relative to vacuum (eV)", fontsize=20)

            return fig, ax



def main():
    # a = Defect(defect_xlsx="hse_screened_qubit_cdft_2021-11-18_2022-01-09")
    # back_process = BackProcess(a.defect_df)
    # back_process.add_tdm_input()
    # back_process.df_to_excel(excel_name="test")

    # a = Defect()
    # a.get_defect_df_v2_hse()
    # back_process = BackProcess(a.defect_df)
    # back_process.add_band_edges_and_defects()
    # back_process.add_level_category()
    # back_process.add_transition_wavevlength()
    # back_process.add_cdft_occupations()
    # back_process.add_hse_fworker()
    # back_process.df_to_excel(excel_name="test")
    #
    # a = Defect(defect_xlsx="defect_2021-11-18")
    # screened = a.get_screened_defect_df()
    # IOTools(cwd=save_xlsx_path, pandas_df=screened).to_excel("triplet_2021-11-18.xlsx")

    # a.statistics(screened_df)
    # a.get_host_df()
    # a.get_host_gp()
    # a.get_defect_df()
    # a.get_defect_df_gp(screened, "triplet_max-tran_gp_2021-11-18")

    # read excel file by IOTools
    df = IOTools(cwd=input_path, excel_file="Table_1_df_2022-04-15").read_excel()
    # from JPack_independent.projects.defectDB.analysis.analysis_api import hse_candidate_df
    Defect.bar_plot_defect_levels_v2(df)

    # cdft = CDFT()
    # cdft.update_entry_with_occ()
    # cdft.get_data_sheet()
    # df1 = cdft.get_zpl_data()
    # IOTools(cwd=os.path.join(DB1_PATH, save_xlsx_path), pandas_df=df1).to_excel("hse_screened_qubit_cdft_zpl_2021-11-18")


if __name__ == '__main__':
    # main()
    # a = Defect(defect_xlsx="hse_screened_qubit_2021-11-18")
    # show = a.defect_df.loc[(a.defect_df["up_tran_en"]==0) & (a.defect_df["dn_tran_en"]!=0)]
    # # print(a.defect_df.dn_tran_lumo_homo_deg.unique())
    # show = show.loc[:, ["task_id", "up_tran_lumo_homo_band_index", "dn_tran_lumo_homo_band_index",
    #                            "up_tran_lumo_homo_deg",
    #                            "dn_tran_lumo_homo_deg"]]
    # # show = show.loc[(show["dn_tran_lumo_homo_deg"] == (1, 2)) | (show["dn_tran_lumo_homo_deg"] == (3, 2))]

    # IOTools(cwd=save_xlsx_path, pandas_df=df).to_excel("test")


    # cdft = CDFT()
    # cdft.update_entry_with_occ()
    # cdft.get_data_sheet()
    # df1 = cdft.get_zpl_data()
    # IOTools(cwd=os.path.join(DB1_PATH, save_xlsx_path), pandas_df=df1).to_excel("hse_screened_qubit_cdft_zpl_2021-11-18")

    # a = [839, 898, 903, 957, 1023, 1056]
    # for i in [2569]:
    #     Tools.extract_defect_levels_v2_hse(i, localisation=0.05, tol=(0.5, 0.5))
    # COHP.analysis()


    # a, b, c, d, e = Tools.extract_defect_levels_v2_hse(2590, localisation=0.2)
    a, b, c, d, e = Tools.extract_defect_levels_v2(89, localisation=0.1)


