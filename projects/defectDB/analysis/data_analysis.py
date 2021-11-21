import os, datetime

from analysis.deep_defect import Potential

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
                False,
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
        # condition = {"task_label": "SCAN_scf", "host_info.c2db_info.uid": "BaBr2-CdI2-NM"}
        condition = {"task_label": "SCAN_scf", "task_id": {"$lte": 6084}}
        for e in list(col.find(condition)):
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
                "defect_pot_a": defect_pot_a[0],
                "site_oxi_state": tuple([tuple(i) for i in e["host_info"]["scan_bs"]["site_oxi_state"]]),
                "number_NN": len(e["NN"])
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


    def get_screened_defect_df(self):
        defect_df = self.defect_df.copy()
        triplet_plt_df = defect_df.loc[(defect_df["mag"] == 2) & (defect_df["level_cbm"] != "None")
                                    & (defect_df["level_vbm"] != "None"), :]
        def triplet_plt():
            return triplet_plt_df

        def qubit_candidate():
            # triplet_df.replace({"up_tran_en": "None"}, 0, inplace=True)
            # triplet_df.replace({"dn_tran_en": "None"}, 0, inplace=True)
            screened_df = triplet_plt_df.loc[(triplet_plt_df["up_tran_en"] >= 0.5) | (triplet_plt_df["dn_tran_en"]
                                                                                       >= 0.5), :]
            return screened_df

        def vacancy():
            screened_df = triplet_plt_df.loc[(triplet_plt_df["defect_type"] == "vacancy"), :]
            return screened_df

        def band_edges():
            screened_df = triplet_plt_df.loc[(triplet_plt_df["is_vbm_cbm_from_same_element"] == False), :]

            return screened_df
        return band_edges()

    def get_defect_df_gp(self, df, output_gp_name):
        # all_triplet_df = defect_df.loc[~defect_df["reduced_site_sym"].isin([("-6m2", "-6m2"), ("3m", "3m"), ("4mm", "4mm")]), :]

        grouping = ["prototype", "pmg_spg", "reduced_site_sym", "reduced_site_specie", "uid", "host_taskid",
                    "site_oxi_state"]
        grouping.extend([
            "vbm_up_max_el", "vbm_up_max_proj_orbital", "vbm_down_max_el",
            "vbm_down_max_proj_orbital", "cbm_up_max_el", "cbm_up_max_proj_orbital",
            "cbm_down_max_el", "cbm_down_max_proj_orbital", "is_vbm_cbm_from_same_element"
        ])
        grouping.extend(["defect_type", "charge", "defect_name", "mag", "number_NN"])
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
        all_protos = ["AgBr3", "BN", "BiI3", "CH", "CdI2", "FeSe", "GaS", "GaSe", "GeS2", "GeSe", "MoS2", "TiCl3"]
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

                defects = []
                taskid_labels = []
                mags =[]
                defect_name_labels = []
                charge_labels = []
                for taskid in plot_df["task_id"].tolist():
                    defect = plot_df.loc[(plot_df["task_id"] == taskid) , :]
                    mag = plot_df["mag"].iloc[0]
                    print("defect_taskid:{}, host_taskid:{}".format(defect["task_id"].iloc[0],
                                                                    defect["host_taskid"].iloc[0]))

                    defects.append(
                        {
                            "cbm": defect["level_cbm"].iloc[0],
                            "vbm": defect["level_vbm"].iloc[0],
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
                    vbm, cbm, up, dn, up_deg, dn_deg, up_ir, dn_ir, up_occ, dn_occ, up_id, dn_id = \
                        defect["vbm"], defect["cbm"], defect["up"], defect["dn"], defect["up_deg"], \
                        defect["dn_deg"], \
                        defect["up_ir"], \
                        defect["dn_ir"], \
                        defect["up_occ"],\
                        defect["dn_occ"],\
                        defect["up_id"], \
                        defect["dn_id"]

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
                            "{}/ {}/ {}/ {}/{}".format(level, occ, deg, ir, band_id) for
                            level, occ, deg, ir, band_id in zip(up, up_occ, up_deg, up_ir, up_id)
                        ]
                    )
                    ax.text(x-1.5, cbm+0.45, up_info, bbox=dict(facecolor='white', edgecolor='black'), size=8)

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
                            "{}/ {}/ {}/ {}/{}".format(level, occupied, deg, ir, band_id) for
                            level, occupied, deg, ir, band_id in zip(dn, dn_occ, dn_deg, dn_ir, dn_id)
                        ]
                    )
                    ax.text(x-1.5, vbm_lim-fig_y_height[0]+0.1, dn_info, bbox=dict(facecolor='white',
                                                                                   edgecolor='black'),
                            size=8)

                ax.yaxis.set_minor_locator(AutoMinorLocator(5))
                for tick in ax.get_yticklabels():
                    tick.set_fontsize(15)

                ax.tick_params(axis="x", bottom=False, labelbottom=True)
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

    def statistics(self, df):
        df = df.loc[df["defect_type"] == "antisite", :]
        df = df.groupby("charge")
        print(df.agg({"task_id": ["count"]}))
        # print(df["task_id"].count())

    def back_process(self):
        def add_number_occ_electrons():
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
            defect_df = self.defect_df.copy()
            defect_df["occ_e"] = defect_df.apply(lambda x: fun(x), axis=1)
            defect_df["triplet_e"] = defect_df.apply(lambda x: x["occ_e"][-1] - x["charge"], axis=1)
            return defect_df

        def band_edges_and_defects():
            pass
        return add_number_occ_electrons()

def main():
    a = Defect(defect_xlsx="defect_2021-11-18")
    # print(a.defect_df["up_in_gap_occ"].iloc[1])
    screened_df = a.get_screened_defect_df()
    # a.statistics(screened_df)
    # a.get_defect_df_v2()
    # a.get_host_df()
    # a.get_host_gp()
    # a.get_defect_df()
    a.get_defect_df_gp(screened_df, "triplet_plt")

    # print(a.defect_df["up_in_gap_occ"])
    # Defect.bar_plot_defect_levels_v2(screened_df)

if __name__ == '__main__':
    print(os.getcwd())
    main()
