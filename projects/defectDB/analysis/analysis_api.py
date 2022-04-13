import pandas as pd
import numpy as np
from qubitPack.tool_box import IOTools, get_db

from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os
from IPython.display import display
from subprocess import check_output
from pymatgen.io.vasp.outputs import Wavecar
from pymatgen.io.vasp.inputs import Poscar, Element
from pathlib import Path

from JPack_independent.projects.defectDB.analysis.data_analysis import Defect
from JPack_independent.projects.defectDB.new_entry_prep.data_preparation import *


MODULE_DIR = Path(__file__).resolve().parent
# This is neccessary for importing this module, while when running this module in shell, comment it first.

input_path = os.path.join(MODULE_DIR, "input")
save_xlsx_path = os.path.join(MODULE_DIR, "output/xlsx")
save_plt_path = os.path.join(MODULE_DIR, "output/plt")

defect_df = IOTools(cwd=input_path, excel_file="defect_2021-11-18_2022-01-10").read_excel()
triplet_df = IOTools(cwd=input_path, excel_file="triplet_2021-11-18").read_excel()
qubit_candidate_df = IOTools(cwd=input_path, excel_file="qubit_2021-11-18_2022-01-10").read_excel()
hse_qubit_df = IOTools(cwd=input_path, excel_file="hse_qubit_2022-04-13").read_excel()
hse_screened_qubit_df = IOTools(cwd=input_path, excel_file="hse_screened_qubit_2021-11-18").read_excel()
hse_screened_qubit_cdft_df = IOTools(cwd=input_path, excel_file="hse_screened_qubit_cdft_2021-11-18_2022-01-31").read_excel()
hse_screened_qubit_cdft_zpl_df = IOTools(cwd=input_path, excel_file="hse_screened_qubit_cdft_zpl_2021-11-18_2022-01-31").read_excel()
defects_36_groups_df = IOTools(cwd=input_path, excel_file="defects_36_groups_2022-03-07").read_excel()
defects_36_groups_zpl_zfs_df = IOTools(cwd=input_path, excel_file="defects_36_groups_zpl_zfs").read_excel()
tmd_antisites_df = IOTools(cwd=input_path, excel_file="tmd_antisites_2022-04-13").read_excel()
hse_candidate_df = IOTools(cwd=input_path, excel_file="Table_6_df_2022-04-07").read_excel()

# os.environ["PATH"] += os.pathsep + "/Library/TeX/texbin/latex"
# plt.style.use(["science", "nature", "light"])


class Toolbox:
    @classmethod
    def get_poscar_and_wavecar(cls, taskid):
        from qubitPack.tool_box import get_db
        col = get_db("HSE_triplets_from_Scan2dDefect", "calc_data-pbe_pc", port=12347).collection
        def scp(taskid):
            for i in col.find({"task_id": taskid}):
                task_id = i["task_id"]
                dir_name = i["dir_name"].split("/")[-1]
                os.makedirs(os.path.join("structures", str(task_id)), exist_ok=True)
                for f in ["POSCAR", "WAVECAR"]:
                    check_output("scp -P 12348 qimin@localhost:{}/{}/{}.gz structures/{}/".format(
                        "/mnt/sdc/tsai/Research/projects/backup/HSE_triplets_from_Scan2dDefect/calc_data-pbe_pc",
                        dir_name, f, task_id).split(" "))
        scp(taskid)
    @classmethod
    def get_parcharg(cls, taskid, band, spin):
        try:
            check_output("gunzip structures/{}/WAVECAR.gz".format(taskid).split(" "), input=b"Y")
            check_output("gunzip structures/{}/POSCAR.gz".format(taskid).split(" "), input=b"Y")
        except Exception:
            pass
        wv = Wavecar("structures/{}/WAVECAR".format(taskid))
        wv.get_parchg(Poscar.from_file("structures/{}/POSCAR".format(taskid)), 0, band-1, spin, phase=True).\
            write_file("structures/{}/{}-{}_spin-{}.vasp".format(taskid, taskid, band, spin))
        check_output("open structures/{}/{}-{}_spin-{}.vasp".format(taskid, taskid, band, spin).split(" "))


    @classmethod
    def plot_defect_levels(cls, task_ids, is_hse):
        if is_hse == "36_group":
            df = defects_36_groups_df.copy()
            print(df.count())
        elif is_hse is True:
            df = hse_qubit_df.copy()
        else:
            df = defect_df.copy()
        Defect.bar_plot_defect_levels_v2(df.loc[df["task_id"].isin(task_ids)])

    @classmethod
    def coop_query(cls):
        from qubitPack.tool_box import get_db
        scan_cohp_db = get_db("Scan2dDefect", "lobster", user="Jeng_ro", password="qimin", port=12347)
        scan_mat_db = get_db("Scan2dMat", "calc_data", user="Jeng_ro", password="qimin", port=12347)


        data = {"pair": [], "task_id": [], "level_cat": [], "defect_name": [], "uid": [], "prev_fw_taskid": []}

        for e in scan_cohp_db.collection.find():
            data["pair"].append(tuple(e["pair"]))
            data["task_id"].append(e["task_id"])
            data["level_cat"].append(e["level_cat"])
            data["defect_name"].append(e["defect_name"])
            pc = scan_mat_db.collection.find_one({"task_id": int(e["pc_from"].split("/")[-1])})
            data["uid"].append(pc["c2db_info"]["uid"])
            data["prev_fw_taskid"].append(e["prev_fw_taskid"])

        return pd.DataFrame(data)

    # df = coop_query()
    # df_gp = df.groupby(["uid", "pair"]).agg({"task_id": "unique","level_cat":"unique", "defect_name": "unique"})


class Task1:
    """
    Mechanism of forming defect levels for cat2 and cat3
    """
    def __init__(self, df=defect_df):
        self.df = df
        self.grouping = ["prototype", "pmg_spg", "pmg_pg", "uid", "host_taskid", "reduced_site_sym",
                         "reduced_site_specie",
                         "charge", "defect_type", "defect_name", "mag", "level_cat"]
        # grouping.extend(["vbm_max_el", "cbm_max_el", "vbm_max_proj_orbital", "cbm_max_proj_orbital", "level_edge_ir"])

    def get_cat_gp(self, prototype=None, chemsys=None, c2db_uid=None, charge=None, level_gap=0, level_cat=None,
                   mag=None, defect_type=None, halogenide=None):
        """
        Stats for neutral defects that are composed by antisites and vacancies from cat 2 and 2'
        :return:
        """
        df = self.df
        if chemsys:
            df = df.loc[df["chemsys"].isin(chemsys)]

        if prototype:
            df = df.loc[df["prototype"].isin(prototype)]

        if c2db_uid:
            df = df.loc[df["uid"].isin(host)]

        if charge:
            df = df.loc[df["charge"].isin(charge)]

        if level_cat:
            df = df.loc[df["level_cat"].isin(level_cat)]

        if mag:
            df = df.loc[df["mag"].isin(mag)]

        if defect_type:
            df = df.loc[df["defect_type"].isin(defect_type)]

        if halogenide is True:
            df = df.loc[(df["chemsys"].str.contains("Cl")) | (df["chemsys"].str.contains("Br")) |
                        (df["chemsys"].str.contains("I-")) | (df["chemsys"].str.contains("-I"))]
        elif halogenide is False:
            df = df.loc[~((df["chemsys"].str.contains("Cl")) | (df["chemsys"].str.contains("Br")) |
                        (df["chemsys"].str.contains("I-")) | (df["chemsys"].str.contains("-I")))]
        elif halogenide is None:
            pass

        df = df.loc[df["level_gap"] >= level_gap]
        df.fillna("None", inplace=True)

        df_gp = df.groupby(self.grouping).agg({"task_id": "unique", "level_gap": "unique", "gap_hse_nosoc": "unique",
                                               "gap_hse": "unique", "mag": "unique",}) #"localisation_threshold":
        # "unique"})
        # IOTools(cwd="output_excel", pandas_df=df_gp).to_excel("task-1_neutral_cat2-2'_defect-gp_2021-11-18", index=True)
        return df_gp


    def get_specie_contribution(self, level_cat=("2", ), halogenide=None):
        '''
        Task
            Make a table of quantity ratio of cation and anion for a defect from different category.
        Conclusion
            75% of hosts in Cat 2 and 3 consist of halogen-based VBM.
        :return:
        '''
        df = self.df.copy()
        def get_sitesym(x):
            defect_specie = x["defect_name"].split("_")[-1]
            reduced_site_specie = x["reduced_site_specie"]
            reduced_site_sym = x["reduced_site_sym"]
            specie_sitesym_dict = dict(zip(reduced_site_specie, reduced_site_sym))
            return (defect_specie, specie_sitesym_dict[defect_specie])

        df["defect_site_sym"] = df.apply(lambda x: get_sitesym(x)[1], axis=1)
        df["defect_specie"] = df.apply(lambda x: get_sitesym(x)[0], axis=1)
        df.fillna("None", inplace=True)
        df = df.loc[df["level_cat"].isin(level_cat)]

        if halogenide is True:
            df = df.loc[(df["chemsys"].str.contains("Cl")) | (df["chemsys"].str.contains("Br")) |
                        (df["chemsys"].str.contains("I-")) | (df["chemsys"].str.contains("-I"))]
        elif halogenide is False:
            df = df.loc[~((df["chemsys"].str.contains("Cl")) | (df["chemsys"].str.contains("Br")) |
                          (df["chemsys"].str.contains("I-")) | (df["chemsys"].str.contains("-I")))]
        elif halogenide is None:
            pass

        defect_specie_df = df.groupby(["defect_specie"]).agg({"task_id": ["unique", "count"], "defect_name": "unique"})
        vbm_specie_df = df.groupby(["vbm_max_el"]).agg({"task_id": ["unique", "count"]})
        return defect_specie_df, vbm_specie_df

class Task4:
    """
    ZPL section
    Task
        Data presentation of relaxation energy
    Conclusion

    """
    @classmethod
    def get_hist_relax_en(cls):
        df = hse_screened_qubit_cdft_zpl_df.copy()
        df = df.loc[df["valid"]==True]
        df = df[(df["ZPL_wavelength"].map(float)<=2500)]

        df["AB/ZPL"] = df["AB"]/df["ZPL"]
        df["CD/ZPL"] = df["CD"]/df["ZPL"]


        ax = df.loc[:, ["AB/ZPL", "CD/ZPL"]].plot.hist(bins=25, edgecolor="black")
        # ax.set_xlim((df1["min-wavelength"].min()-0.1, df1["min-wavelength"].max()+0.1))
        ax.set_xlim(0.75, 1.25)
        ax.set_ylim(0, 50)

        ax.tick_params(which="both", bottom=True, top=False, right=False)
        # ax.tick_params(which="minor", left=False)
        ax.set_xlabel("Ratio of relaxation energy and ZPL")
        ax.set_ylabel("Count")
        ax.legend(labels=[r"$\mathrm{\frac{Absorption}{ZPL}}$", r"$\mathrm{\frac{Emission}{ZPL}}$"])
        ax.set_title("Distribution of ratio of relaxation energy and ZPL")
        plt.show()


class Task3Qubit:
    """
    Pick nice candidates for promising qubit and emitter
    """
    def __init__(self, df=hse_screened_qubit_cdft_zpl_df.loc[hse_screened_qubit_cdft_zpl_df["valid"]==True]):
        self.df = df

    @property
    def data_preparation(self):
        df = self.df.copy()

        def add_chemsys(x):
            try:
                gs_taskid = x["task_id"]
            except:
                gs_taskid = x["gs_taskid"]

            chemsys = hse_screened_qubit_df.loc[hse_screened_qubit_df["task_id"] == gs_taskid, "chemsys"].iloc[0]
            return chemsys
        # df["chemsys"] = df.apply(lambda x: add_chemsys(x), axis=1)

        def add_host_d_f_nature(x):
            chemsys = x["chemsys"]
            is_host_d_f_element = None
            for el in chemsys.split("-"):
                host_el = Element(el)
                if not host_el.is_transition_metal and not host_el.is_rare_earth_metal:
                    is_host_d_f_element = False
                else:
                    is_host_d_f_element = True
                    break
            return is_host_d_f_element
        # df = df.loc[~df["defect"].str.contains("vac")]
        # df["host_has_d_f_element"] = df.apply(lambda x: add_host_d_f_nature(x), axis=1)

        def add_zfs_data(x):
            zfs_collection = get_db("HSE_triplets_from_Scan2dDefect", "zfs_data-pbe_pc", port=12347).collection
            try:
                gs_taskid = x["task_id"]
            except:
                gs_taskid = x["gs_taskid"]
            zfs_e = zfs_collection.find_one({"prev_fw_taskid": gs_taskid})
            D, E = "None", "None"
            try:
                D, E = zfs_e["pyzfs_out"]["D"], zfs_e["pyzfs_out"]["E"]
            except Exception:
                print(gs_taskid)
            return D, E
        df["zfs_D"] = df.apply(lambda x: add_zfs_data(x)[0], axis=1)
        df["zfs_E"] = df.apply(lambda x: add_zfs_data(x)[1], axis=1)
        return df

    def get_non_halogenides(self, df):
        non_halogen_df = df.loc[~(
                (df["chemsys"].str.contains("I-")) |
                (df["chemsys"].str.contains("-I")) |
                (df["chemsys"].str.contains("Br")) |
                (df["chemsys"].str.contains("Cl"))
        )]
        return non_halogen_df

    def get_halogenides(self, df):
        halogen_df = df.loc[(df["chemsys"].str.contains("I-")) |
                                      (df["chemsys"].str.contains("-I")) |
                                      (df["chemsys"].str.contains("Br")) |
                                      (df["chemsys"].str.contains("Cl"))]
        return halogen_df

    @property
    def get_filtered_df(self):
        """
        Task
            Qubit defect levels
        Conditions
            1. ZPL: >= 0.8 eV (1,550 nm/C-band) (NV-: 1.946 eV (638 nm) @ 10 K)
            2. ZFS: >= 1 GHz (NV-: 2.87 GHz)
            3. Light atom as antisite atom (not d or f elements)

        :return:
        """
        df = self.data_preparation
        df = df.loc[(df["zfs_D"] != "None") & (df["zfs_E"] != "None")]

        df["zfs_D"] = df.apply(lambda x: round(x["zfs_D"]/1000, 2), axis=1)
        df["ZPL_wavelength"] = df.apply(lambda x: round(x["ZPL_wavelength"]), axis=1)

        df = df.loc[abs(df["zfs_E"].map(float)) < 1]
        df = df.loc[abs(df["zfs_D"].map(float)) >= 1]
        df = df.loc[df["ZPL"] >= 0.8]
        return df

    @property
    def show_non_halogen_gp(self):
        df = self.get_non_halogenides(self.get_filtered_df)
        df_gp = df.groupby(["host_has_d_f_element",
                            "prototype", "charge", "host",
                            "chemsys", "defect"]).agg(
            {"ZPL":["unique"], "zfs_D": ["unique"],
             "ZPL_wavelength": "unique", "gs_taskid": "unique"})
        # IOTools(cwd=".", pandas_df=non_halogen_gp).to_excel("draft_hse_selective_nonhalogen", index=True)
        return df_gp

    @property
    def picked_non_halogenides(self):
        df = self.get_non_halogenides(self.get_filtered_df)
        picked_non_halogen_df = df.loc[df["gs_taskid"].isin(
            [929, 162, 575, 151, 182, 171, 186, 211, 229, 196, 347])]
        picked_non_halogen_gp = picked_non_halogen_df.groupby(
            ["prototype", "charge", "host", "chemsys", "defect"]).agg(
            {"ZPL":["unique"], "zfs_D": ["unique"],
             "ZPL_wavelength": "unique", "total_TDM_rate": "unique",  "gs_taskid": "unique"}
        )
        return picked_non_halogen_df, picked_non_halogen_gp

    @property
    def picked_halogenides(self):
        df = self.get_halogenides(self.get_filtered_df)
        df.groupby(["host_has_d_f_element", "prototype", "charge", "host", "defect"]).agg(
            {"ZPL":["unique"],
             "zfs_D": ["unique"], "ZPL_wavelength": "unique",
             "gs_taskid": "unique"
            })
        picked_halogen_df = df.loc[(df["host_has_d_f_element"] == False) & (~df[
            "prototype"].isin(["TiCl3", "AgBr3"]))]
        picked_halogen_gp = picked_halogen_df.groupby(["prototype", "charge", "host", "defect"]).agg(
            {"ZPL":["unique"],
             "zfs_D": ["unique"], "ZPL_wavelength": "unique",
             "gs_taskid": "unique"
            })
        return picked_halogen_df, picked_halogen_gp

    def get_plot_correlation_zfs_zpl(self):
        h = self.picked_halogenides[0].copy()
        h["ZFS-D"] = h["zfs_D"]
        nh = self.picked_non_halogenides[0].copy()
        nh["ZFS-D"] = nh["zfs_D"]

        plt.scatter(x=nh["ZPL"], y=nh["ZFS-D"], label="Non-halogenides")
        plt.scatter(x=h["ZPL"], y=h["ZFS-D"], label="Halogenides")
        plt.xlabel("ZPL (eV)")
        plt.ylabel("ZFS-D (GHz)")
        plt.legend()
        plt.show()


class Task3Emitter(Task3Qubit):
    def __init__(self):
        super(Task3Emitter, self).__init__()
        self.grouping = ["prototype", "charge", "host", "defect", "total_TDM_rate", "ZPL", "ZPL_wavelength",
                         "up_TDM_rate_X", "up_TDM_rate_Y", "up_TDM_rate_Z",
                         "dn_TDM_rate_X", "dn_TDM_rate_Y", "dn_TDM_rate_Z"]
    @property
    def get_filtered_df(self):
        """
        Conditions
            1. Total TDM transition >= 7 MHz (check paper review: "Low dipole transition rate”)
            2. ZPL >= 0.8
            4. host_has_d_f_element is FALSE
        :return:
        """
        df = self.data_preparation
        df = df.loc[(df["zfs_D"] != "None") & (df["zfs_E"] != "None")]
        df["zfs_D"] = df.apply(lambda x: round(x["zfs_D"]/1000, 2), axis=1)
        df["ZPL_wavelength"] = df.apply(lambda x: round(x["ZPL_wavelength"]), axis=1)
        df["AB/ZPL"] = df["AB"]/df["ZPL"]
        df["CD/ZPL"] = df["CD"]/df["ZPL"]

        df = df[df["host_has_d_f_element"]==False] # reduce intersystem-crossing rate
        df = df[(df["ZPL"] >= 0.8)]
        # df1 = df1[df1["ZPL_wavelength"] <= 2500]
        df = df.loc[df["total_TDM_rate"] >= 7]
        return df

    @property
    def get_large_TDM_rate(self):
        lg_rate_df = self.get_filtered_df()
        lg_rate_df.sort_values(["total_TDM_rate", "ZPL"], ascending=False)
        lg_rate_gp = lg_rate_df.groupby(self.grouping).agg({"gs_taskid": ["unique", "count"], "defect": "unique"})
        return lg_rate_df, lg_rate_gp

    @property
    def picked_non_halogenides(self):
        gs_taskids = [109, 349, 151, 242, 196, 314, 917, 35, 975, 924, 306, 1234, 1423]
        df = self.get_non_halogenides(self.get_filtered_df)
        df = df.loc[df["gs_taskid"].isin(gs_taskids)]
        df_gp = df.groupby(self.grouping).agg({"gs_taskid": ["unique", "count"], "zfs_E": "unique"})
        display(df)
        display(df_gp)
        return df, df_gp

    @property
    def picked_halogenides(self):
        gs_taskids = [109, 349, 151, 242, 196, 314, 917, 35, 975, 924, 306, 1234, 1423]
        df = self.get_halogenides(self.get_filtered_df)
        df = df.loc[df["gs_taskid"].isin(gs_taskids)]
        df_gp = df.groupby(self.grouping).agg({"gs_taskid": ["unique", "count"], "zfs_E": "unique"})
        display(df)
        display(df_gp)
        return df, df_gp


class Task5(Task3Qubit):
    def __init__(self):
        super(Task5, self).__init__()

    def extra_data_preparation(self, df):
        df["source-specie"] = df.apply(
            lambda x: hse_screened_qubit_df.loc[
                hse_screened_qubit_df["task_id"] == x["gs_taskid"], "level_source_specie"].iloc[0], axis=1)
        df["non-source-specie"] = df.apply(
            lambda x: x["chemsys"].split("-")[1]
            if x["chemsys"].split("-").index(x["source-specie"])==0 else x["chemsys"].split("-")[0], axis=1)

        df["source-group"] = df.apply(lambda x: int(Element(x["source-specie"]).group), axis=1)
        df["non-source-group"] = df.apply(lambda x: int(Element(x["non-source-specie"]).group), axis=1)
        df["Z"] = df.apply(lambda x: int(Element(x["source-specie"]).Z), axis=1)
        return df

    def get_zpl_plots(self, non_source_group):
        df = self.data_preparation
        df = self.extra_data_preparation(df)
        zpl_df = df.loc[(df["ZPL"] <= 15)]

        zpl_df = zpl_df.loc[zpl_df["non-source-group"] == non_source_group]

        display(zpl_df.groupby(["source-group", "Z", "prototype", "host","defect", "source-specie",
                          "non-source-specie"]).agg({"gs_taskid": ["unique", "count"]}))

        groups = zpl_df["source-group"].unique().tolist()

        fig, axes = plt.subplots(len(groups[:]), 1, figsize=(10,12))
        for idx, group in enumerate(groups[:]):
            print(group)
            plot_df = zpl_df.loc[zpl_df["source-group"] == group]
            x_name = "non-source-specie"
            try:
                ax = axes[idx]
            except Exception:
                ax = axes
            for proto in plot_df["source-specie"].unique():
                tgt_df = plot_df.loc[(plot_df["source-specie"] == proto) & ~(plot_df["gs_taskid"].isin([934]))]
                x = [int(Element(el).row) for el in tgt_df[x_name].tolist()]
                ax.scatter(x=x, y=tgt_df["ZPL"], label="Source:{}, Source-specie:{}".format(group, proto))
                for x, y, host, defect, charge, gs_taskid, prototype in zip(x, tgt_df["ZPL"], tgt_df["host"], tgt_df["defect"], tgt_df["charge"], tgt_df["gs_taskid"], tgt_df["prototype"]):
                    ax.text(x, y, "{},{},{},{},{},{}".format(x, host, "-".join(defect.split("_")), charge, gs_taskid, prototype))


            ax.set_xlabel(x_name+" row")
            ax.set_ylabel("ZPL (eV)")
            ax.set_title("Non-source specie: {}".format(non_source_group))
            ax.legend()

        plt.show()

    def get_focus_zpl_plots(self, non_source_group, source_group):
        df = self.data_preparation
        zpl_df = df.loc[(self.df["ZPL"]<=15)]

        zpl_df = zpl_df.loc[(zpl_df["non-source-group"] == non_source_group) & (zpl_df["source-group"] == source_group)]

        display(zpl_df.groupby(["source-group", "Z", "prototype", "host","defect", "source-specie", "non-source-specie"]).agg({"gs_taskid": ["unique", "count"]}))

        groups = zpl_df["source-specie"].unique().tolist()

        fig, axes = plt.subplots(len(groups[:]), 1, figsize=(10,10))
        for idx, group in enumerate(groups[:]):
            print(group)
            plot_df = zpl_df.loc[zpl_df["source-specie"] == group]
            x_name = "non-source-specie"
            try:
                ax = axes[idx]
            except Exception:
                ax = axes
            for proto in plot_df["prototype"].unique():
                tgt_df = plot_df.loc[(plot_df["prototype"] == proto)]
                x = [int(Element(el).row) for el in tgt_df[x_name].tolist()]
                ax.scatter(x=x, y=tgt_df["ZPL"], label="Source:{}, Prototype:{}".format(group, proto))
                for x, y, host, defect, charge, gs_taskid, prototype in zip(x, tgt_df["ZPL"], tgt_df["host"], tgt_df["defect"], tgt_df["charge"], tgt_df["gs_taskid"], tgt_df["prototype"]):
                    ax.text(x, y, "{},{},{},{},{},{}".format(x, host, "-".join(defect.split("_")), charge, gs_taskid, prototype))


            ax.set_xlabel(x_name+" row")
            ax.set_ylabel("ZPL (eV)")
            ax.set_title("source specie: {}".format(group))
            ax.legend()

        plt.show()

    def get_zfs_plots(self, non_source_group):
        df = self.data_preparation
        zfs_df = df.loc[(self.df["zfs_D"] != "None") & (df["zfs_E"] != "None")]
        zfs_df = zfs_df.loc[zfs_df["non-source-group"] == non_source_group]
        zfs_df["zfs_D"] = zfs_df.apply(lambda x: round(abs(x["zfs_D"]/1000), 2), axis=1)

        display(zfs_df.groupby(["source-group", "Z", "prototype", "host","defect",
                                "source-specie", "non-source-specie"]).agg({"gs_taskid": ["unique", "count"], "zfs_D": "unique"}))

        groups = zfs_df["source-group"].unique().tolist()

        fig, axes = plt.subplots(len(groups[:]), 1, figsize=(10,12))
        for idx, group in enumerate(groups[:]):
            print(group)
            plot_df = zfs_df.loc[zfs_df["source-group"] == group]
            x_name = "non-source-specie"
            try:
                ax = axes[idx]
            except Exception:
                ax = axes
            for proto in plot_df["source-specie"].unique():
                tgt_df = plot_df.loc[(plot_df["source-specie"] == proto) & ~(plot_df["gs_taskid"].isin([934]))]
                x = [int(Element(el).row) for el in tgt_df[x_name].tolist()]
                ax.scatter(x=x, y=tgt_df["zfs_D"], label="Source:{}, Source-specie:{}".format(group, proto))
                for x, y, host, defect, charge, gs_taskid, prototype in zip(x, tgt_df["zfs_D"], tgt_df["host"], tgt_df["defect"], tgt_df["charge"], tgt_df["gs_taskid"], tgt_df["prototype"]):
                    ax.text(x, y, "{},{},{},{},{},{}".format(x, host, "-".join(defect.split("_")), charge, gs_taskid, prototype))


            ax.set_xlabel(x_name+" row")
            ax.set_ylabel("ZFS (GHz)")
            ax.set_title("Non-source specie: {}".format(non_source_group))
            ax.legend()

        plt.show()

    def get_focus_zfs_plots(self, non_source_group, source_group):
        df = self.data_preparation
        zfs_df = df.loc[(df["zfs_D"] != "None") & (df["zfs_E"] != "None")]
        zfs_df["zfs_D"] = zfs_df.apply(lambda x: round(abs(x["zfs_D"]/1000), 2), axis=1)

        zfs_df = zfs_df.loc[(zfs_df["non-source-group"] == non_source_group) & (zfs_df["source-group"] == source_group)]

        display(zfs_df.groupby(["source-group", "Z", "prototype", "host","defect", "source-specie", "non-source-specie"]).agg({"gs_taskid": ["unique", "count"]}))

        groups = zfs_df["source-specie"].unique().tolist()

        fig, axes = plt.subplots(len(groups[:]), 1, figsize=(10,10))
        for idx, group in enumerate(groups[:]):
            print(group)
            plot_df = zfs_df.loc[zfs_df["source-specie"] == group]
            x_name = "non-source-specie"
            try:
                ax = axes[idx]
            except Exception:
                ax = axes
            for proto in plot_df["prototype"].unique():
                tgt_df = plot_df.loc[(plot_df["prototype"] == proto)]
                x = [int(Element(el).row) for el in tgt_df[x_name].tolist()]
                ax.scatter(x=x, y=tgt_df["zfs_D"], label="Source:{}, Prototype:{}".format(group, proto))
                for x, y, host, defect, charge, gs_taskid, prototype in zip(x, tgt_df["zfs_D"], tgt_df["host"], tgt_df["defect"], tgt_df["charge"], tgt_df["gs_taskid"], tgt_df["prototype"]):
                    ax.text(x, y, "{},{},{},{},{},{}".format(x, host, "-".join(defect.split("_")), charge, gs_taskid, prototype))


            ax.set_xlabel(x_name+" row")
            ax.set_ylabel("ZFS (GHz)")
            ax.set_title("source specie: {}".format(group))
            ax.legend()

        plt.show()

    def get_zpl_nonsource_halogenides(self):
        self.get_zpl_plots(17)

    def get_zpl_nonsource_nonhalogenides(self):
        els = [16, 15, 14, 13, 2, 1, 3, 12]
        for el in els:
            self.get_zpl_plots(el)

    def get_zfs_nonsource_halogenides(self):
        self.get_zfs_plots(17)

    def get_zfs_nonsource_nonhalogenides(self):
        els = [16, 15, 14, 13, 2, 1, 3, 12]
        for el in els:
            self.get_zfs_plots(el)


class Task6(Task3Qubit):
    def __init__(self):
        super(Task6, self).__init__(df=defect_df)

    @property
    def get_well_def_df(self):
        df = self.df.loc[
            (defect_df["level_cbm"] != "None") & (defect_df["level_vbm"] != "None") & (defect_df["level_gap"] != 0) &
            (defect_df["up_in_gap_band_id"] != ()) & (defect_df["dn_in_gap_band_id"] != ()) &
            (defect_df["up_in_gap_band_id"] != None) & (defect_df["dn_in_gap_band_id"] != None)]
        return df

    @property
    def data_preparation(self):
        df = self.get_well_def_df

        def add_top_level_to_edge(x):
            level_cbm = x["level_cbm"]
            level_vbm = x["level_vbm"]
            level_gap = x["level_gap"]

            up_in_gap_level = x["up_in_gap_level"]
            up_in_gap_occ = x["up_in_gap_occ"]
            dn_in_gap_level = x["dn_in_gap_level"]
            dn_in_gap_occ = x["dn_in_gap_occ"]

            high_unocc_level, high_unocc_level_to_vbm = None, None
            up_unocc_level = [en for en, occ in zip(up_in_gap_level, up_in_gap_occ) if occ == 0]
            dn_unocc_level = [en for en, occ in zip(dn_in_gap_level, dn_in_gap_occ) if occ == 0]
            if up_unocc_level+dn_unocc_level != []:
                high_unocc_level = max(up_unocc_level+dn_unocc_level)
                high_unocc_level_to_vbm = high_unocc_level-level_vbm

            low_occ_level, low_occ_level_to_cbm = None, None
            up_occ_level = [en for en, occ in zip(up_in_gap_level, up_in_gap_occ) if occ != 0]
            dn_occ_level = [en for en, occ in zip(dn_in_gap_level, dn_in_gap_occ) if occ != 0]
            if up_occ_level + dn_occ_level != []:
                low_occ_level = min(up_occ_level+dn_occ_level)
                low_occ_level_to_cbm = level_cbm - low_occ_level


            max_up_in_gap_level = max(up_in_gap_level)
            max_dn_in_gap_level = max(dn_in_gap_level)
            max_in_gap_level = max(max_up_in_gap_level, max_dn_in_gap_level)

            d_level_to_vbm = max_in_gap_level - level_vbm
            d_level_to_cbm = level_cbm - max_in_gap_level
            top_level_to_edge = None
            dist_from_vbm_ratio = 0.20
            if 0 < d_level_to_vbm <= dist_from_vbm_ratio*level_gap:
                top_level_to_edge = "vbm"
            elif  dist_from_vbm_ratio*level_gap < d_level_to_vbm <= (1-dist_from_vbm_ratio)*level_gap:
                top_level_to_edge = "mid"
            elif (1-dist_from_vbm_ratio)*level_gap < d_level_to_vbm <= level_gap:
                top_level_to_edge = "cbm"
            return top_level_to_edge, high_unocc_level, low_occ_level, high_unocc_level_to_vbm, low_occ_level_to_cbm

        df["top_level_close_to"] = df.apply(lambda x: add_top_level_to_edge(x)[0], axis=1)
        df["high_unocc_level"] = df.apply(lambda x: add_top_level_to_edge(x)[1], axis=1)
        df["low_occ_level"] = df.apply(lambda x: add_top_level_to_edge(x)[2], axis=1)
        df["high_unocc_level_to_vbm"] = df.apply(lambda x: add_top_level_to_edge(x)[3], axis=1)
        df["low_occ_level_to_cbm"] = df.apply(lambda x: add_top_level_to_edge(x)[4], axis=1)

        def add_elec_neg(x):
            from pymatgen.core.molecular_orbitals import MolecularOrbitals
            uid = x["uid"]
            formula = uid.split("-")[0]
            mo = MolecularOrbitals(formula)
            elec_neg = mo.elec_neg
            return elec_neg

        df["elec_neg"] = df.apply(lambda x: add_elec_neg(x), axis=1)

        df["non_source_specie"] = df.apply(
            lambda x:
            x["chemsys"].split("-")[1] if x["chemsys"].split("-").index(x["level_source_specie"])==0
            else x["chemsys"].split("-")[0], axis=1
        )

        df["level_from_anion_cation"] = df.apply(
            lambda x: "cation" if Element(x["level_source_specie"]).X < Element(x["non_source_specie"]).X
            else "anion", axis=1)
        return df

    @property
    def get_cation_vac_df(self):
        df = self.data_preparation
        cation_vac_df = df.loc[(self.df["level_from_anion_cation"] == "anion") & (df["defect_type"] == "vacancy")]
        cation_vac_df = cation_vac_df.fillna("None")
        return cation_vac_df

    @property
    def get_anion_vac_df(self):
        df = self.data_preparation
        anion_vac_df = df.loc[(self.df["level_from_anion_cation"] == "cation") & (df["defect_type"] == "vacancy")]
        anion_vac_df = anion_vac_df.fillna("None")
        return anion_vac_df

    @property
    def get_anion_antisite_df(self):
        df = self.data_preparation
        anion_antisite_df = df.loc[(df["level_from_anion_cation"] == "cation") & (df["defect_type"] == "antisite")]
        anion_antisite_df = anion_antisite_df.fillna("None")
        return anion_antisite_df

    @property
    def get_cation_antisite_df(self):
        df = self.data_preparation
        cation_antisite_df = df.loc[(df["level_from_anion_cation"] == "anion") & (df["defect_type"] == "antisite")]
        cation_antisite_df = cation_antisite_df.fillna("None")
        return cation_antisite_df

    def get_plot_lowest_donor_electroneg_in_anion_vac(self):
        """
        Get plot of the lowest donor state to electronegativity for anion vacancy
        :return:
        """
        df = self.get_anion_vac_df
        df = df.loc[(df["low_occ_level_to_cbm"] != "None") & (df["elec_neg"] != "None")]

        df["low_occ_level_to_cbm_r"] = df["low_occ_level_to_cbm"] / df["level_gap"]

        fig, axes = plt.subplots(3,2, figsize=(7,6))

        for idx, charge in enumerate([-1, 0, 1]):
            charge_df = df.loc[df["charge"] == charge]
            x_data = charge_df["elec_neg"]
            y_data = charge_df["low_occ_level_to_cbm_r"]
            z_data = charge_df["charge"]

            ax = axes[idx][0]
            ax.scatter(x_data, y_data, label="q={}".format(charge))
            ax.set_ylabel("low-occ-level-to-cbm/gap")
            ax.set_xlabel("elec-neg")
            ax.legend()

            from scipy.optimize import curve_fit
            y_f = lambda x, a, b: a*x+b
            popt, pcov = curve_fit(y_f, x_data, y_data, p0=(0, 0))
            print(popt, pcov)
            ax.plot(np.linspace(min(x_data), max(x_data), 100),
                    [y_f(x, popt[0], popt[1]) for x in np.linspace(min(x_data), max(x_data), 100)], color="r")

            ax = axes[idx][1]
            im = ax.imshow(np.abs(pcov))
            fig.colorbar(im, ax=ax, orientation='vertical')
        plt.show()

    def get_plot_highest_acceptor_electroneg_in_cation_vac(self):
        """
         Get plot of the highest acceptor state to electronegativity for cation vacancy
         :return:
         """
        df = self.get_cation_vac_df
        df = df.loc[(df["high_unocc_level_to_vbm"] != "None") & (df["elec_neg"] != "None")]
        df["high_unocc_level_to_vbm_r"] = df["high_unocc_level_to_vbm"] / df["level_gap"]

        fig, axes = plt.subplots(3,2, figsize=(7,6))
        for idx, charge in enumerate([-1, 0, 1]):
            charge_df = df.loc[df["charge"] == charge]
            x_data = charge_df["elec_neg"]
            y_data = charge_df["high_unocc_level_to_vbm_r"]
            z_data = charge_df["charge"]

            ax = axes[idx][0]
            ax.scatter(x_data, y_data, label="q={}".format(charge))
            ax.set_ylabel("high-unocc-level-to-vbm/gap")
            ax.set_xlabel("elec-neg")
            ax.legend()

            from scipy.optimize import curve_fit
            y_f = lambda x, a, b: a*x+b
            popt, pcov = curve_fit(y_f, x_data, y_data, p0=(0, 0))
            print(popt, pcov)
            ax.plot(np.linspace(min(x_data), max(x_data), 100), [y_f(x, popt[0], popt[1]) for x in np.linspace(min(x_data), max(x_data), 100)], color="r")

            ax = axes[idx][1]
            im = ax.imshow(np.abs(pcov))
            fig.colorbar(im, ax=ax, orientation='vertical')
        plt.show()

    def get_plot_lowest_donor_electroneg_in_anion_antisite(self):
        df = self.get_anion_antisite_df
        df = df.loc[(df["low_occ_level_to_cbm"] != "None") & (df["elec_neg"] != "None")]
        df["low_occ_level_to_cbm_r"] = df["low_occ_level_to_cbm"] / df["level_gap"]

        fig, axes = plt.subplots(3,2, figsize=(7,6))

        for idx, charge in enumerate([-1, 0, 1]):
            charge_df = df.loc[df["charge"] == charge]
            x_data = charge_df["elec_neg"]
            y_data = charge_df["low_occ_level_to_cbm_r"]
            z_data = charge_df["charge"]

            ax = axes[idx][0]
            im = ax.scatter(x_data, y_data, label="q={}".format(charge))
            ax.set_ylabel("low-occ-level-to-cbm/gap")
            ax.set_xlabel("elec-neg")
            ax.legend()

            from scipy.optimize import curve_fit
            y_f = lambda x, a, b: a*x+b
            popt, pcov = curve_fit(y_f, x_data, y_data, p0=(0, 0))
            print(popt, pcov)
            ax.plot(np.linspace(min(x_data), max(x_data), 100), [y_f(x, popt[0], popt[1]) for x in np.linspace(min(x_data), max(x_data), 100)], color="r")

            ax = axes[idx][1]
            im = ax.imshow(np.abs(pcov))
            fig.colorbar(im, ax=ax, orientation='vertical')
        plt.show()

    def get_plot_highest_acceptor_electroneg_in_cation_antisite(self):
        df = self.get_cation_antisite_df
        df = df.loc[(df["high_unocc_level_to_vbm"] != "None") & (df["elec_neg"] != "None")]
        df["high_unocc_level_to_vbm_r"] = df["high_unocc_level_to_vbm"] / df["level_gap"]

        fig, axes = plt.subplots(3,2, figsize=(7,6))

        for idx, charge in enumerate([-1, 0, 1]):
            charge_df = df.loc[df["charge"] == charge]
            x_data = charge_df["elec_neg"]
            y_data = charge_df["high_unocc_level_to_vbm_r"]
            z_data = charge_df["charge"]

            ax = axes[idx][0]
            im = ax.scatter(x_data, y_data, label="q={}".format(charge))
            ax.set_ylabel("high-unocc-level-to-vbm/gap")
            ax.set_xlabel("elec-neg")
            ax.legend()


            from scipy.optimize import curve_fit
            y_f = lambda x, a, b: a*x+b
            popt, pcov = curve_fit(y_f, x_data, y_data, p0=(0, 0))
            print(popt, pcov)
            ax.plot(np.linspace(min(x_data), max(x_data), 100), [y_f(x, popt[0], popt[1]) for x in np.linspace(min(x_data), max(x_data), 100)], color="r")

            ax = axes[idx][1]
            im = ax.imshow(np.abs(pcov))
            fig.colorbar(im, ax=ax, orientation='vertical')
        plt.show()

class Defect36Group:
    def __init__(self, df=defects_36_groups_df):
        self.df = df

    @property
    def data_preparation(self):
        df = self.df.copy()
        def add_zfs_data(x):
            zfs_collection = get_db("HSE_triplets_from_Scan2dDefect", "zfs_data-pbe_pc", port=12347).collection
            try:
                gs_taskid = x["task_id"]
            except:
                gs_taskid = x["gs_taskid"]
            zfs_e = zfs_collection.find_one({"prev_fw_taskid": gs_taskid})
            D, E = "None", "None"
            try:
                D, E = zfs_e["pyzfs_out"]["D"], zfs_e["pyzfs_out"]["E"]
            except Exception:
                print(gs_taskid)
            return D, E


        def add_gap_diff():
            db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12345, user="readUser", password="qiminyan")
            for uid in df["uid"].unique():
                print(uid)
                gap_hse = db.collection.find_one({"uid": uid})["gap_hse"]
                gap_hse_nosoc = db.collection.find_one({"uid": uid})["gap_hse_nosoc"]
                df.loc[df["uid"] == uid, "gap_hse"] = gap_hse
                df.loc[df["uid"] == uid, "gap_hse_nosoc"] = gap_hse_nosoc
                df.loc[df["uid"] == uid, "gap_diff"] = df.loc[df["uid"] == uid, "level_gap"] - gap_hse_nosoc
            df["gap_diff"] = df["gap_diff"].round(3)
            return df

        df = add_gap_diff()
        df["zfs_D"] = df.apply(lambda x: add_zfs_data(x)[0], axis=1)
        df["zfs_E"] = df.apply(lambda x: add_zfs_data(x)[1], axis=1)
        return df

class SheetCollection:
    def __init__(self):
        self.scan_defect_df = None
        self.scan_defect_df_gp = None

        self.scan_triplet_df = None
        self.scan_triplet_df_gp = None

        self.scan_pre_candidate_df = None
        self.scan_pre_candidate_df_gp = None

        self.hse_triplet_df = None
        self.hse_triplet_df_gp = None

        self.hse_candidate_df = None
        self.hse_candidate_df_gp = None

        self.grouping = [
            "prototype", "spacegroup", "C2DB_uid", "good_site_sym", "good_site_specie",
            "defect_type", "defect_name", "charge", "perturbed_vbm", "perturbed_cbm",
            "perturbed_bandgap", "level_edge_category",
            "vertical_transition_up",  "vertical_transition_down"
        ]

        self.grouping_hse_candidate = [
            "prototype", "spacegroup", "C2DB_uid", "defect_type", "defect_name", "charge",
            "perturbed_bandgap","zfs_D", "zfs_E",  "vertical_transition_up",
            "vertical_transition_down", "transition_from", "gap_tran_bottom_to_vbm", "gap_tran_top_to_cbm",
            "ZPL", "AB", "BC", "CD", 'DA',
            "up_TDM_rate_X", "up_TDM_rate_Y", "up_TDM_rate_Z", "dn_TDM_rate_X", "dn_TDM_rate_Y", "dn_TDM_rate_Z",
            "up_polarization", "dn_polarization", "allowed"
        ]


    def get_diff_btw_two_dfs(self, df1, df2):
        a = df1.copy()
        b = df2.copy()
        df_diff = pd.concat([a.loc[:, ["uid", "defect_name", "charge"]],
                             b.loc[:, ["uid", "defect_name", "charge"]]],
                            axis=0).drop_duplicates(keep=False)
        return df_diff

    def cat_revision(self, x):
        level_cat = x["level_cat"]
        if level_cat == "2'":
            level_cat = 2
        elif level_cat == "3'":
            level_cat = 3
        return level_cat

    def check_ingap_level(self, df, vbm_distance=-0.25, cbm_distance=0.25):
        df["in_gap_up"] = df.apply(lambda x: np.all(np.array(x["up_in_gap_level"]) - x["level_vbm"] >= vbm_distance)
                                             and np.all(np.array(x["up_in_gap_level"]) - x["level_cbm"] <= cbm_distance),
                                   axis=1)
        df["in_gap_dn"] = df.apply(lambda x: np.all(np.array(x["dn_in_gap_level"]) - x["level_vbm"] >= vbm_distance)
                                         and np.all(np.array(x["dn_in_gap_level"]) - x["level_cbm"] <= cbm_distance),
                                   axis=1)

        df["in_gap"] = df.apply(lambda x: x["in_gap_up"] or x["in_gap_dn"], axis=1)
        return df

    def check_vertical_transition(self, df, vertical_transition=0.5):
        df["exist_vertical_transition"] = df.apply(
            lambda x: x["vertical_transition_up"] >= vertical_transition or x["vertical_transition_down"] >=
                      vertical_transition, axis=1)
        return df

    def get_scan_defect_df(self):
        df = defect_df.copy()
        df["level_cat"] = df.apply(lambda x: self.cat_revision(x), axis=1)

        df["spacegroup"] = df["pmg_spg"]
        df["good_site_sym"] = df["reduced_site_sym"]
        df["good_site_specie"] = df["reduced_site_specie"]
        df["C2DB_uid"] = df["uid"]

        df["perturbed_vbm"] = df["level_vbm"]
        df["perturbed_cbm"] = df["level_cbm"]
        df["perturbed_bandgap"] = df["level_gap"]
        df["level_edge_category"] = df["level_cat"]
        df["vertical_transition_up"] = df["up_tran_en"]
        df["vertical_transition_down"] = df["dn_tran_en"]

        df = df.fillna("None")

        # Remove bad bandedges
        df = df.loc[(df["perturbed_vbm"] != "None") & (df["perturbed_cbm"] != "None") ]

        # remove duplicated defects (same defect name, same prototype, same spacegroup, same charge, same host_taskid)
        # and keep the one with the highest task_id (the one with the most recent scan) and remove the other one (the
        # one with the older scan)  (this is to avoid duplicated defects in the df)
        df = df.groupby(["defect_name", "prototype", "spacegroup", "C2DB_uid", "host_taskid",  "charge", ]).apply(
            lambda x: x.loc[x["task_id"].idxmax()])
        df = df.reset_index(drop=True)

        self.scan_defect_df = df
        self.scan_defect_df_gp = self.scan_defect_df.groupby(self.grouping).agg({"task_id": "unique"})

        return self.scan_defect_df
    def get_scan_triplet_df(self, df=None):
        df = self.get_scan_defect_df().copy() if df is None else df
        self.scan_triplet_df = df.loc[df["mag"] == 2]
        self.scan_triplet_df_gp = self.scan_triplet_df.groupby(self.grouping).agg({"task_id": "unique"})

        return self.scan_triplet_df

    def get_scan_pre_candidate_df(self, df=None):
        df = self.get_scan_triplet_df().copy() if df is None else df
        # get df that the vertical transition up is larger than or equal to 0.5 eV or the vertical transition down is
        # larger than or equal to 0.5 eV
        df = self.check_vertical_transition(df, 0.5)
        self.scan_pre_candidate_df = df.loc[df["exist_vertical_transition"]]
        self.scan_pre_candidate_df_gp = self.scan_pre_candidate_df.groupby(self.grouping).agg({"task_id": "unique"})

        return self.scan_pre_candidate_df

    def get_hse_triplet_df(self):
        df = hse_qubit_df.copy()
        # print number of defects of df
        print("Number of defects in hse_qubit_df: {}".format(len(df)))
        no_36group_df = df.loc[
            ~(
                    (df["uid"].isin(defects_36_groups_df["uid"])) &
                    (df["defect_name"].isin(defects_36_groups_df["defect_name"])) &
                    (df["charge"].isin(defects_36_groups_df["charge"]))
            )]
        print(f"{len(no_36group_df)} defects in 36 groups")
        no_36group_and_no_6group_tmd_df = no_36group_df.loc[
            ~(
                    (no_36group_df["uid"].isin(tmd_antisites_df["uid"])) &
                    (no_36group_df["defect_name"].isin(tmd_antisites_df["defect_name"])) &
                    (no_36group_df["charge"].isin(tmd_antisites_df["charge"]))
            )]

        print(f"{len(no_36group_and_no_6group_tmd_df)} defects in 36 groups and no 6 group tmd")
        clean_36_group_df = defects_36_groups_df.loc[defects_36_groups_df["task_id"] != 246]
        df = pd.concat([no_36group_and_no_6group_tmd_df, clean_36_group_df, tmd_antisites_df])
        df["level_cat"] = df.apply(lambda x: self.cat_revision(x), axis=1)

        df["spacegroup"] = df["pmg_spg"]
        df["good_site_sym"] = df["reduced_site_sym"]
        df["good_site_specie"] = df["reduced_site_specie"]
        df["C2DB_uid"] = df["uid"]

        df["perturbed_vbm"] = df["level_vbm"]
        df["perturbed_cbm"] = df["level_cbm"]
        df["perturbed_bandgap"] = df["level_gap"]
        df["level_edge_category"] = df["level_cat"]
        df["vertical_transition_up"] = df["up_tran_en"]
        df["vertical_transition_down"] = df["dn_tran_en"]
        df = df.fillna("None")

        # Remove bad bandedges
        df = df.loc[(df["perturbed_vbm"] != "None") & (df["perturbed_cbm"] != "None")]

        # Make sure that the definition of in-gap levels are followed (it is no a filter, but a definition)
        df = self.check_ingap_level(df, -0.25, 0.25)
        df = df.loc[df["in_gap"]]

        self.hse_triplet_df = df
        self.hse_triplet_df_gp = df.groupby(self.grouping).agg({"task_id": ["unique", "count"]})

        return self.hse_triplet_df

    def get_hse_pre_candidate_df(self, df=None):
        #note that this sheet can have duplicates
        df = self.get_hse_triplet_df().copy() if df is None else df
        df = self.check_vertical_transition(df, 0.5)
        self.hse_pre_candidate_df = df.loc[df["exist_vertical_transition"]]

        self.hse_pre_candidate_df_gp = self.hse_pre_candidate_df.groupby(self.grouping).agg({"task_id": ["unique",
                                                                                                         "count"]})

        return self.hse_pre_candidate_df


    def get_hse_candidate_df(self, df=None):
        df = self.get_hse_pre_candidate_df().copy() if df is None else df
        def get_zpl_zfs_df(df):
            task3 = Task3Qubit(df)
            zfs_df = task3.data_preparation
            zfs_df = zfs_df.set_index("task_id")

            zpl = DataPrepCDFT(defect_entry_df=df)
            zpl.one_shot_zpl_df()
            zpl_df = zpl.zpl_df
            zpl_df = zpl_df.set_index("task_id")

            zfs_zpl_df = zfs_df.join(zpl_df)
            zfs_zpl_df = zfs_zpl_df.reset_index()
            def fun(x):
                zfs_D = x["zfs_D"]
                zfs_E = x["zfs_E"]

                if type(zfs_D) == float:
                    zfs_D = round(zfs_D/1000, 2)
                if type(zfs_E) == float:
                    zfs_E = round(abs(zfs_E), 0)

                return zfs_D, zfs_E


            zfs_zpl_df["zfs_D"] = zfs_zpl_df.apply(lambda x: fun(x)[0], axis=1)
            zfs_zpl_df["zfs_E"] = zfs_zpl_df.apply(lambda x: fun(x)[1], axis=1)
            return zfs_zpl_df

        df = get_zpl_zfs_df(df)
        def get_filtered_df(df, D_zfs=1, wavelength_zpl=2500, homo_to_vmb=0.1, lumo_to_cbm=0.1):
            def func(x):
                tran_from = x["transition_from"]
                print(x["task_id"])
                if tran_from == "dn":
                    gap_tran_top_to_cbm = x["level_gap"] - x["dn_tran_top"]
                    print(x["task_id"], x["transition_from"], x["level_gap"], x["dn_tran_top"])
                else:
                    gap_tran_top_to_cbm = x["level_gap"] - x["up_tran_top"]
                    print(x["task_id"], x["transition_from"], x["level_gap"], x["up_tran_top"])
                return gap_tran_top_to_cbm

            df = df.loc[df["zfs_D"] != "None"]
            df = df.loc[df["zfs_D"] >= D_zfs]

            df = df.loc[df["allowed"] == True]
            df = df.loc[df["ZPL"] != "None"]
            df = df.loc[df["ZPL_wavelength"] <= wavelength_zpl]

            df["gap_tran_bottom_to_vbm"] = df.apply(
                lambda x: x["dn_tran_bottom"] if x["transition_from"] == "dn" else x["up_tran_bottom"], axis=1)
            df["gap_tran_top_to_cbm"] = df.apply(
                lambda x: x["level_gap"] - x["dn_tran_top"] if x["transition_from"] == "dn" else x["level_gap"] - x["up_tran_top"], axis=1)

            df = df.loc[df["gap_tran_bottom_to_vbm"] >= homo_to_vmb]
            df = df.loc[df["gap_tran_top_to_cbm"] >= lumo_to_cbm]
            return df

        df = get_filtered_df(df)
        # Exclude some double countings
        def exclude_double_counting(df):
            double_counting_taskids_list = [186, 283, 211, 229, 242]
            df = df.loc[~(df["task_id"].isin(double_counting_taskids_list))]

            df = df.groupby(["defect_name", "prototype", "spacegroup", "C2DB_uid", "host_taskid", "charge",
                             "transition_from"]).apply(
                lambda x: x.loc[x["zfs_D"].astype(float).idxmax()])
            df = df.reset_index(drop=True)
            return df

        self.hse_candidate_df = exclude_double_counting(df)
        self.hse_candidate_df_gp = self.hse_candidate_df.groupby(self.grouping_hse_candidate).agg(
            {"task_id": ["unique", "count"]})
        return df


class FigureAndStats(SheetCollection):
    def __init__(self):
        super().__init__()
        self.figures = []

    def get_zpl_fig(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        # get a histogram plot of ZPL_wave values for df
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #rename the columns ZPL_wavelength to ZPL-wavelength
        df.rename(columns={"ZPL_wavelength": "ZPL-wavelength"}, inplace=True)
        ax.hist(df["ZPL-wavelength"], bins=25, edgecolor="black")
        ax.set_xlabel("ZPL-wavelength")
        ax.set_ylabel("Counts")
        ax.set_title("ZPL-wavelength histogram")
        self.figures.append(fig)
        return fig

    def get_zpl_stat(self, df=None):
        df1 = self.get_hse_candidate_df().copy() if df is None else df.copy()
        df1 = df1.loc[:, "ZPL_wavelength"]
        # ref Appl. Phys. Rev. 7, 031308 (2020); doi: 10.1063/5.0006075
        visible = df1.loc[(df1 >= 380) & (df1 < 700)]
        nir = df1.loc[(df1 >= 700) & (df1 <= 2500)]

        bio1 = df1.loc[(df1 >= 650) & (df1 <= 950)] #650
        bio2 = df1.loc[(df1 >= 1000) & (df1 <= 1350)]

        O_band = df1.loc[(df1 >= 1260) & (df1 < 1360)]
        E_band = df1.loc[(df1 >= 1360) & (df1 < 1460)]
        S_band = df1.loc[(df1 >= 1460) & (df1 < 1530)]
        C_band = df1.loc[(df1 >= 1530) & (df1 < 1565)]
        other_band = df1.loc[(df1 >= 1565) & (df1 <= 2500)]

        set_visible_df = pd.DataFrame(
            {
                "name": ["visible", "nir"],
                "frequency": [i.count() for i in [visible, nir]],
            }
        )
        set_visible_df.set_index("name", inplace=True)
        set_visible_df["fraction"] = round(set_visible_df.frequency / set_visible_df.frequency.sum()*100, 2)
        display(set_visible_df)


        set_bio_df = pd.DataFrame(
            {
                "name": ["bio1", "bio2"],
                "frequency": [i.count() for i in [bio1, bio2]],
            }
        )
        set_bio_df.set_index("name", inplace=True)
        set_bio_df["fraction"] = round(set_bio_df.frequency / set_visible_df.loc["nir", "frequency"]*100, 2)
        # set_bio_df["fraction"] = round(set_bio_df.frequency / set_visible_df.frequency.sum()*100, 2)
        set_bio_df.loc["sum", :] = set_bio_df.sum()


        display(set_bio_df)

        set_tele_df = pd.DataFrame(
            {
                "name": ["O_band", "E_band", "S_band", "C_band", "other_band"],
                "frequency": [i.count() for i in [O_band, E_band, S_band, C_band, other_band]],
            }
        )
        set_tele_df.set_index("name", inplace=True)
        set_tele_df["fraction"] = round(set_tele_df.frequency / set_visible_df.loc["nir", "frequency"]*100, 2)
        # set_tele_df["fraction"] = round(set_tele_df.frequency /  set_visible_df.frequency.sum()*100, 2)
        set_tele_df.loc["sum", :] = set_tele_df.iloc[0:3].sum()
        display(set_tele_df)




    def get_zfs_fig(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        # get a histogram plot of ZFS_D values for df
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        ax.hist(df["zfs_D"], bins=20)
        ax.set_xlabel("ZFS_D")
        ax.set_ylabel("Count")
        ax.set_title("ZFS_D")
        self.figures.append(fig)
        return fig

    def test_valid_zpl(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        df["valid_zpl"] = df.apply(lambda x: x["ZPL"] <= x["vertical_transition_up"] if x["transition_from"] == "up"
        else x["ZPL"] <= x["vertical_transition_down"], axis=1)
        df["diff_zpl_vertical_transition"] = df.loc[df["valid_zpl"] == False].apply(lambda x: x["ZPL"] - x[
            "vertical_transition_up"] if x["transition_from"] == "up" else x["ZPL"] - x["vertical_transition_down"], axis=1)

        return df.loc[df["valid_zpl"] == False,
                      ["task_id",  "C2DB_uid", "ZPL", "vertical_transition_up", "vertical_transition_down",
                       "diff_zpl_vertical_transition", "AB", "BC", "CD", "DA", "transition_from", "up_tran_cdft_occ",
                       "dn_tran_cdft_occ", "dn_in_gap_deg", "up_in_gap_deg"]]



if __name__ == '__main__':
    pass



