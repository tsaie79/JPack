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

# hse_qubit_df = IOTools(cwd=input_path, excel_file="hse_qubit_2022-04-13").read_excel()
hse_qubit_df = IOTools(cwd=input_path, excel_file="hse_qubit_2022-04-21").read_excel()

tmd_antisites_df = IOTools(cwd=input_path, excel_file="tmd_antisites_2022-04-13").read_excel()

hse_screened_qubit_df = IOTools(cwd=input_path, excel_file="hse_screened_qubit_2021-11-18").read_excel()
hse_screened_qubit_cdft_df = IOTools(cwd=input_path, excel_file="hse_screened_qubit_cdft_2021-11-18_2022-01-31").read_excel()
hse_screened_qubit_cdft_zpl_df = IOTools(cwd=input_path, excel_file="hse_screened_qubit_cdft_zpl_2021-11-18_2022-01-31").read_excel()

# defects_36_groups_df = IOTools(cwd=input_path, excel_file="defects_36_groups_2022-03-07").read_excel()
defects_36_groups_df = IOTools(cwd=input_path, excel_file="defects_36_groups_2022-04-27").read_excel()


defects_36_groups_zpl_zfs_df = IOTools(cwd=input_path, excel_file="defects_36_groups_zpl_zfs").read_excel()

# hse_candidate_df = IOTools(cwd=input_path, excel_file="Table_5_df_2022-04-15").read_excel()
hse_candidate_df = IOTools(cwd=input_path, excel_file="Table_5_df_2022-04-27").read_excel()

table1_df = IOTools(cwd=input_path, excel_file="Table_1_df_2022-04-15").read_excel()



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
                # print(f"\n=={gs_taskid} has ZFS data==")
                # print(f"\nD: {D}, E: {E}")
            except Exception:
                print(f"\n=={gs_taskid} has no ZFS data==")
            return D, E
        df["zfs_D"] = df.apply(lambda x: add_zfs_data(x)[0], axis=1)
        df["zfs_E"] = df.apply(lambda x: add_zfs_data(x)[1], axis=1)
        print(f"\nZFS data added to dataframe:\n{df.loc[:, ['task_id', 'zfs_D', 'zfs_E']]}")
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
            "LUMO-HOMO_up",  "LUMO-HOMO_down"
        ]

        self.grouping_hse_candidate = [
            "prototype", "spacegroup", "C2DB_uid", "defect_type", "defect_name", "charge",
            "perturbed_bandgap","zfs_D", "zfs_E",  "LUMO-HOMO_up", "LUMO-HOMO_down", "LUMO_HOMO_deg", "transition_from",
            "gap_tran_bottom_to_vbm", "gap_tran_top_to_cbm",
            "ZPL", "AB", "BC", "CD", 'DA',
            "up_TDM_rate_X", "up_TDM_rate_Y", "up_TDM_rate_Z", "dn_TDM_rate_X", "dn_TDM_rate_Y", "dn_TDM_rate_Z",
            "TDM_rate", "up_polarization", "dn_polarization", "allowed", "class"
        ]


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

    def get_zpl_zfs_df(self, df):
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
        self.zfs_zpl_df = zfs_zpl_df
        return zfs_zpl_df

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
        df["LUMO-HOMO_up"] = df["up_tran_en"]
        df["LUMO-HOMO_down"] = df["dn_tran_en"]

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

    def get_hse_triplet_df(self, df=None, remove_duplicates_based_on=["uid", "defect_name", "charge"],
                           remove_taskid_in_36group=[186, 1020, 1000, 283, 211, 644, 229, 548, 242, 603]):

        df = df if df is not None else hse_qubit_df.copy()
        # print number of defects of df
        print("Number of defects in hse_qubit_df: {}".format(len(df)))
        no_36group_df = df.loc[
            ~(
                    (df["uid"].isin(defects_36_groups_df["uid"])) &
                    (df["defect_name"].isin(defects_36_groups_df["defect_name"])) &
                    (df["charge"].isin(defects_36_groups_df["charge"]))
            )]
        print(f"{len(no_36group_df)} defects not in 36-groups")
        no_36group_and_no_6group_tmd_df = no_36group_df.loc[
            ~(
                    (no_36group_df["uid"].isin(tmd_antisites_df["uid"])) &
                    (no_36group_df["defect_name"].isin(tmd_antisites_df["defect_name"])) &
                    (no_36group_df["charge"].isin(tmd_antisites_df["charge"]))
            )]

        print(f"{len(no_36group_and_no_6group_tmd_df)} defects not in 36-groups and not in 6-group tmd")
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
        df["LUMO-HOMO_up"] = df["up_tran_en"]
        df["LUMO-HOMO_down"] = df["dn_tran_en"]

        df = df.fillna("None")

        # Remove bad bandedges
        df = df.loc[(df["perturbed_vbm"] != "None") & (df["perturbed_cbm"] != "None")]

        # Make sure that the definition of in-gap levels are followed (it is no a filter, but a definition)
        df = self.check_ingap_level(df, -0.50, 0.50)
        # df = df.loc[df["in_gap"]]

        if remove_taskid_in_36group:
            df = df.loc[~df["task_id"].isin(remove_taskid_in_36group)]

        if remove_duplicates_based_on:
            df = df.drop_duplicates(
                subset=remove_duplicates_based_on, keep="first").sort_values("task_id", ascending=True)
        self.hse_triplet_df = df
        self.hse_triplet_df_gp = df.groupby(self.grouping).agg({"task_id": ["unique", "count"]})

        return self.hse_triplet_df

    def get_hse_pre_candidate_df(self, df=None):
        #note that this sheet can have duplicates
        df = self.get_hse_triplet_df(remove_duplicates_based_on=None).copy() if df is None else df
        df = self.check_vertical_transition(df, 0.5)
        self.hse_pre_candidate_df = df.loc[df["exist_vertical_transition"]]

        self.hse_pre_candidate_df_gp = self.hse_pre_candidate_df.groupby(self.grouping).agg({"task_id": ["unique",
                                                                                                         "count"]})

        return self.hse_pre_candidate_df


    def get_hse_candidate_df(self, df=None):
        df = self.get_hse_pre_candidate_df().copy() if df is None else df

        df = self.get_zpl_zfs_df(df)
        df["LUMO_HOMO_deg"] = df.apply(lambda x: x["dn_tran_lumo_homo_deg"] if x["transition_from"] == "dn" else x[
            "up_tran_lumo_homo_deg"], axis=1)
        df["TDM_rate"] = df.apply(lambda x: x["total_TDM_rate_up"] if x["transition_from"] == "up" else x[
            "total_TDM_rate_dn"], axis=1)

        def fun(x):
            if x["transition_from"] == "dn":
                gap_tran_bottom_to_vbm = x["dn_tran_bottom"]
                gap_tran_top_to_cbm = x["level_gap"] - x["dn_tran_top"]
            elif x["transition_from"] == "up":
                gap_tran_bottom_to_vbm = x["up_tran_bottom"]
                gap_tran_top_to_cbm = x["level_gap"] - x["up_tran_top"]
            else:
                if x["dn_tran_bottom"] == "None":
                    gap_tran_bottom_to_vbm = x["up_tran_bottom"]
                    gap_tran_top_to_cbm = x["level_gap"] - x["up_tran_top"]
                else:
                    gap_tran_bottom_to_vbm = x["dn_tran_bottom"]
                    gap_tran_top_to_cbm = x["level_gap"] - x["dn_tran_top"]
            return gap_tran_bottom_to_vbm, gap_tran_top_to_cbm

        df["gap_tran_bottom_to_vbm"] = df.apply(lambda x: fun(x)[0], axis=1)
        df["gap_tran_top_to_cbm"] = df.apply(lambda x: fun(x)[1], axis=1)

        # df["gap_tran_bottom_to_vbm"] = df.apply(
        #     lambda x: x["dn_tran_bottom"] if x["transition_from"] == "dn" else x["up_tran_bottom"], axis=1)
        # df["gap_tran_top_to_cbm"] = df.apply(
        #     lambda x: x["level_gap"] - x["dn_tran_top"] if x["transition_from"] == "dn" else x["level_gap"] - x["up_tran_top"], axis=1)

        def get_filtered_df(df, D_zfs=0.95, wavelength_zpl=2500, homo_to_vmb=0.095, lumo_to_cbm=0.095):
            df = df.loc[df["zfs_D"] != "None"]
            df = df.loc[df["zfs_D"] >= D_zfs]

            df = df.loc[df["allowed"] == True]
            df = df.loc[df["ZPL"] != "None"]
            df = df.loc[df["ZPL_wavelength"] <= wavelength_zpl]

            df = df.loc[df["gap_tran_bottom_to_vbm"] >= homo_to_vmb]
            df = df.loc[df["gap_tran_top_to_cbm"] >= lumo_to_cbm]

            problematic_taskids = FigureAndStats().test_valid_opt_tran(df).loc[:, "task_id"]
            df = df.loc[~(df["task_id"].isin(problematic_taskids))]
            return df

        passed_df = get_filtered_df(df)
        not_passed_df = df.loc[~(df["task_id"].isin(passed_df["task_id"]))]

        # Exclude some double countings
        def exclude_double_counting(df):
            double_counting_taskids_list = [186, 283, 211, 229, 242]
            df = df.loc[~(df["task_id"].isin(double_counting_taskids_list))]
            df = df.groupby(["defect_name", "prototype", "spacegroup", "C2DB_uid", "host_taskid", "charge",
                             "transition_from"]).apply(lambda x: x.loc[x["zfs_D"].astype(float).idxmax()])
            df = df.reset_index(drop=True)
            return df


        passed_df = exclude_double_counting(passed_df)
        def exclude_weird_zfs(df):
            #exclude task_id = 990 from df
            df = df.loc[df["task_id"] != 990]
            return df

        passed_df = exclude_weird_zfs(passed_df)
        self.hse_candidate_df = passed_df.replace(np.nan, "None").copy()
        self.hse_candidate_df_gp = self.hse_candidate_df.copy().replace(np.nan, "None").groupby(
            self.grouping_hse_candidate).agg({"task_id": "unique"})
        self.hse_candidate_not_passed_df = not_passed_df.replace(np.nan, "None").copy()
        self.hse_candidate_not_passed_df_gp = self.hse_candidate_not_passed_df.copy().replace(np.nan, "None").groupby(
            self.grouping_hse_candidate).agg(
            {"task_id": "unique"})


class FigureAndStats(SheetCollection):
    def __init__(self):
        super().__init__()
        self.figures = []

    def get_grand_defect_level_plot(self, df, sym_name, ylim):
        # plot a bar chart containing the defect levels for a mutiple defects
        plot_df = df if df.shape[0] > 0 else self.hse_candidate_df

        defects = []
        taskid_labels = []
        mags =[]
        defect_name_labels = []
        charge_labels = []
        host_labels = []
        print(f"task_id :{plot_df['task_id'].tolist()}")
        for idx, taskid in enumerate(plot_df["task_id"].tolist()):
            print("=="*20, "\n", f"taskid: {taskid}, idx: {idx}")
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
            if int(defect["charge"]) == -1:
                charge_labels.append("1-")
            elif int(defect["charge"]) == 1:
                charge_labels.append("1+")
            else:
                charge_labels.append(defect["charge"])

            host_labels.append(defect["C2DB_uid"])

        vbm_lim, cbm_lim = plot_df["level_vbm"].median(), plot_df["level_cbm"].median()

        fig, ax = plt.subplots(figsize=(12,11), dpi=300)
        fig_y_height = 5, 5
        font_size = 7.25
        ax.set_ylim(vbm_lim-fig_y_height[0], cbm_lim+fig_y_height[1])

        print("len of defects: {}".format(len(defects)))
        # plot mutilple bars
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
            up_tran_bottom = defect["up_tran_bottom"] + defect["level_vbm"] if defect["up_tran_bottom"] != "None" else \
                None
            up_tran_en = defect["up_tran_en"]
            dn_tran_bottom = defect["dn_tran_bottom"] + defect["level_vbm"] if defect["dn_tran_bottom"] != 'None' else None
            dn_tran_en = defect["dn_tran_en"]
            transition_from = defect["transition_from"]
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
            for level, occ, deg, id in zip(up, up_occ, up_deg, up_band):
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
                    "{}/{}/{}/{}/{}".format(level, occ, deg, ir, band_id) for
                    level, occ, deg, ir, band_id in zip(up, up_occ, up_deg, up_ir, up_band)
                ]
            )
            ax.text(x-1.5, cbm+0.75, up_info, bbox=dict(facecolor='white', edgecolor='black'), size=font_size)
            ax.text(x-1.5, cbm+0.4, "{}/{}/{}/{}".format(cbm_max_el, cbm_max_proj_orbital, level_edge_ir[0][1], cbm),
                    bbox=dict(
                        facecolor='white', edgecolor='black'), size=font_size)
            if up_tran_en and (up_tran_bottom or up_tran_bottom != "None") and transition_from == "up":
                ax.arrow(x=x-1.25, y=up_tran_bottom, dx=0, dy=up_tran_en, width=0.02,
                         length_includes_head=True, facecolor="gray")

            dx = 1
            for level, occ, deg, id in zip(dn, dn_occ, dn_deg, dn_band):
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
                    "{}/{}/{}/{}/{}".format(level, occupied, deg, ir, band_id) for
                    level, occupied, deg, ir, band_id in zip(dn, dn_occ, dn_deg, dn_ir, dn_band)
                ]
            )
            ax.text(x-1.5, vbm-1.05, dn_info, bbox=dict(facecolor='white', edgecolor='black'),
                    size=font_size)
            ax.text(x-1.5, vbm-0.55, "{}/{}/{}/{}".format(vbm_max_el, vbm_max_proj_orbital,
                                                          level_edge_ir[0][0], vbm),
                    bbox=dict(facecolor='white', edgecolor='black'), size=font_size)
            # ax.text(x-1.5, vbm_lim-fig_y_height[0]+3, dn_info, bbox=dict(facecolor='white', edgecolor='black'),
            #         size=font_size)
            # ax.text(x-1.5, vbm_lim-fig_y_height[0]+2.5, "{}/{}/{}/{}".format(vbm_max_el, vbm_max_proj_orbital,
            #                                                            level_edge_ir[0][0], vbm),
            #         bbox=dict(facecolor='white', edgecolor='black'), size=font_size)
            if dn_tran_en and dn_tran_bottom and transition_from == "dn":
                ax.arrow(x=x+1.25, y=dn_tran_bottom, dx=0, dy=dn_tran_en, width=0.02,
                         length_includes_head=True, facecolor="gray")

        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        for tick in ax.get_yticklabels():
            tick.set_fontsize(15)

        ax.tick_params(axis="x", bottom=False, labelbottom=True)
        ax.tick_params(axis="y", right=False, which="both")


        ax.set_xticks(np.arange(0, len(defects)*5, 5))
        # ax.set_xticklabels(["{}".format(taskid_labels[i]) for i in range(len(mags))],
        #                    fontdict={"fontsize": 10}, rotation=0)

        ax.set_xticklabels([rf"$\mathrm{{{defect_name_labels[i].split('_')[2]}_{{"
                            rf"{defect_name_labels[i].split('_')[4]}}}^{{{charge_labels[i]}}}}}$"+
                            rf" @ {''.join(host_labels[i].split('-')[0].split('2'))}"+ "\n"+rf"{taskid_labels[i]}" for
                            i in range(len(mags))], fontdict={"fontsize": 10}, rotation=0)

        # remove minor and major xticks
        ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=True)
        print(f"taskid_label: {taskid_labels}, mags: {mags}")

        ax.set_title(f"Antisite defects in PTMC with {sym_name}", fontsize=15)
        ax.set_ylabel("Energy relative to vacuum (eV)", fontsize=15)
        ax.set_ylim(ylim[0], ylim[1])

        fig.savefig("/Users/jeng-yuantsai/Research/project/Scan2dDefect/manuscript/fig/test.pdf", format="pdf")
        return fig

    def get_zpl_fig(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        # get a histogram plot of ZPL_wave values for df
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #rename the columns ZPL_wavelength to ZPL-wavelength
        df.rename(columns={"ZPL_wavelength": "ZPL-wavelength"}, inplace=True)

        import matplotlib.cm as cm
        c = "darkorange"
        ax.axvspan(650, 950, facecolor=c, alpha=0.5)
        # ax.axvline(650, color="wheat", linestyle="-")
        ax.text(650, 9.3, "650", rotation=90, verticalalignment="center", horizontalalignment="center")
        # ax.axvline(950, color="wheat", linestyle="-")
        ax.text(950, 9.3, "950", rotation=90, verticalalignment="center", horizontalalignment="center")

        c = "darkgreen"
        ax.axvspan(1000, 1350, facecolor=c, alpha=0.5)
        # ax.axvline(1000, color="limegreen", linestyle="-")
        ax.text(1000, 9.3, "1000", rotation=90, verticalalignment="center", horizontalalignment="center")
        # ax.axvline(1350, color="limegreen", linestyle="-")
        ax.text(1350, 9.3, "1350", rotation=90, verticalalignment="center", horizontalalignment="center")

        for wv in range(380, 780, 20):
            cmap = cm.get_cmap('jet')
            c = cmap(float(wv - 380) / (780 - 380))
            ax.axvspan(wv, wv + 20, facecolor=c, alpha=0.5)
        # ax.axvline(380, color="yellow", linestyle="-")
        ax.text(380, 9.3, "380", rotation=90, verticalalignment="center", horizontalalignment="center")
        # ax.axvline(780, color="yellow", linestyle="-")
        ax.text(780, 9.3, "780", rotation=90, verticalalignment="center", horizontalalignment="center")

        c = "grey"
        ax.axvspan(1260, 1565, facecolor=c, alpha=0.5)
        # ax.axvline(1260, color=c, linestyle="-")
        ax.text(1260, 9.3, "1260", rotation=90, verticalalignment="center", horizontalalignment="center")
        # ax.axvline(1630, color=c, linestyle="-")
        ax.text(1565, 9.3, "1565", rotation=90, verticalalignment="center", horizontalalignment="center")

        # ax.hist(df["ZPL-wavelength"], bins=25, edgecolor="black")
        df1 = df.loc[:, "ZPL-wavelength"]
        # plot bar chart of ZPL-wavelength with range of 250 to 2000 with 25
        ax.bar(
            [i+25 for i in range(250, 2000, 50)],
            [df1.loc[(df1 >= i) & (df1 < i + 50)].count() for i in range(250, 2000, 50)],
            width=50,
            edgecolor="black",
        )
        # set count of ZPL-wavelength at the top of each bar except the zero-height bar
        for i in range(250, 2000, 50):
            if df1.loc[(df1 >= i) & (df1 < i + 50)].count() != 0:
                ax.text(i + 25, df1.loc[(df1 >= i) & (df1 < i + 50)].count(),
                        str(df1.loc[(df1 >= i) & (df1 < i + 50)].count()),
                        horizontalalignment="center",
                        verticalalignment="bottom")


        ax.set_xlim(250, 2000)
        ax.set_ylim(0, 10)
        ax.set_xlabel("ZPL wavelength (nm)")
        ax.set_ylabel("Frequency")
        ax.set_title("ZPL wavelength distribution")
        # remove subticks
        ax.tick_params(axis="x", which="minor", bottom=True, top=False, labelbottom=False, direction="out")
        ax.tick_params(axis="x", which="major", bottom=True, top=False, direction="out")
        ax.tick_params(axis="y", which="minor", left=False, right=False, labelbottom=True)
        ax.tick_params(axis="y", which="major", left=True, right=False, labelbottom=True, direction="out")

        self.figures.append(fig)
        return fig

    def get_zpl_stat(self, df=None):
        df1 = self.get_hse_candidate_df().copy() if df is None else df.copy()
        df1 = df1.loc[:, "ZPL_wavelength"]
        # ref Appl. Phys. Rev. 7, 031308 (2020); doi: 10.1063/5.0006075
        visible = df1.loc[(df1 >= 380) & (df1 < 780)]
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
                "wavelength": ["380-780", "700-2500"],
            }
        )
        set_visible_df.set_index("name", inplace=True)
        set_visible_df["fraction"] = round(set_visible_df.frequency / df1.count()*100)
        display(set_visible_df)


        set_bio_df = pd.DataFrame(
            {
                "name": ["bio1", "bio2"],
                "frequency": [i.count() for i in [bio1, bio2]],
                "wavelength": ["650-950", "1000-1350"],
            }
        )
        set_bio_df.set_index("name", inplace=True)
        set_bio_df["fraction"] = round(set_bio_df.frequency / df1.count()*100)
        set_bio_df.loc["sum", :] = set_bio_df.sum()


        display(set_bio_df)

        set_tele_df = pd.DataFrame(
            {
                "name": ["O_band", "E_band", "S_band", "C_band"],
                "frequency": [i.count() for i in [O_band, E_band, S_band, C_band]],
                "wavelength": ["1260-1360", "1360-1460", "1460-1530", "1530-1565"],
            }
        )
        set_tele_df.set_index("name", inplace=True)
        set_tele_df["fraction"] = round(set_tele_df.frequency / df1.count()*100)
        set_tele_df.loc["sum", :] = set_tele_df.iloc[:].sum()
        display(set_tele_df)

    def get_cdft_relax_energy_fig(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        # plot histogram of CDFT relaxation energy BC and DA
        fig = plt.figure()
        ax = fig.add_subplot(111)
        df1 = df.loc[:, "BC"]
        # plot bar chart range from 0 to 0.5 with interval 0.01 of BC and DA

        ax.bar(
            np.arange(0, 0.4-0.01, 0.01)+0.01/2,
            df1.value_counts(bins=np.arange(0, 0.4, 0.01)).sort_index(),
            width=0.01,
            label="BC",
            edgecolor="black",
        )

        # set text of count for each bin at the top of the bar if the count != 0, else remove the text
        for i in range(len(df1.value_counts(bins=np.arange(0, 0.4, 0.01)).sort_index())):
            if df1.value_counts(bins=np.arange(0, 0.4, 0.01)).sort_index().iloc[i] != 0:
                ax.text(
                    np.arange(0, 0.4-0.01, 0.01)[i]+0.01/2,
                    df1.value_counts(bins=np.arange(0, 0.4, 0.01)).sort_index().iloc[i],
                    df1.value_counts(bins=np.arange(0, 0.4, 0.01)).sort_index().iloc[i],
                    horizontalalignment="center",
                    verticalalignment="bottom",
                    )

        # remove subticks
        ax.tick_params(axis="x", which="minor", bottom=True, top=False, labelbottom=False, direction="out")
        ax.tick_params(axis="x", which="major", bottom=True, top=False, direction="out")
        ax.tick_params(axis="y", which="major", left=True, right=False, labelbottom=True, direction="out")
        ax.tick_params(axis="y", which="minor", left=False, right=False, labelbottom=True, direction="out")
        ax.set_xlim(0, 0.4)
        ax.set_ylim(0, 18)


        #plot inset bar chart of DA
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        axins = inset_axes(ax, width="50%", height="50%", loc="upper center")
        df1 = df.loc[:, "DA"]
        axins.bar(
            np.arange(0, 0.2-0.01, 0.01)+0.01/2,
            df1.value_counts(bins=np.arange(0, 0.2, 0.01)).sort_index(),
            width=0.01,
            label="DA",
            edgecolor="black",
            color="orange"
        )

        # set text of count for each bin at the top of the bar if the count != 0, else remove the text
        for i in range(len(df1.value_counts(bins=np.arange(0, 0.2, 0.01)).sort_index())):
            if df1.value_counts(bins=np.arange(0, 0.2, 0.01)).sort_index().iloc[i] != 0:
                axins.text(
                    np.arange(0, 0.2-0.01, 0.01)[i]+0.01 if i < 2 else np.arange(0, 0.2-0.01, 0.01)[i]+0.01/2,
                    df1.value_counts(bins=np.arange(0, 0.2, 0.01)).sort_index().iloc[i],
                    df1.value_counts(bins=np.arange(0, 0.2, 0.01)).sort_index().iloc[i],
                    horizontalalignment="center",
                    verticalalignment="bottom",
                    )

        axins.legend(loc="upper right")
        axins.tick_params(axis="x", which="minor", bottom=True, top=False, labelbottom=False, direction="out")
        axins.tick_params(axis="x", which="major", bottom=True, top=False, direction="out")
        axins.tick_params(axis="y", which="major", left=True, right=False, labelbottom=True, direction="out")
        axins.tick_params(axis="y", which="minor", left=False, right=False, labelbottom=True, direction="out")
        axins.set_xlim(0, 0.2)
        axins.set_ylim(0, 30)

        ax.set_ylim(0, 18)
        ax.set_xlabel("Relaxation energy (eV)")
        ax.set_ylabel("Frequency")
        ax.set_title("CDFT relaxation energy distribution")
        ax.legend(loc="upper right")

        self.figures.append(fig)
        return fig

    def get_cdft_relax_energy_stat(self, df=None):
        """
        Get statistics of CDFT relaxation energy distribution
        """
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        # plot histogram of CDFT relaxation energy BC and DA
        bc, da = 0.05, 0.05
        print("BC<={}, N={}-{}%, DA<={}, N={}-{}%".format(bc, df.loc[df["BC"] <= bc].count().iloc[0],
                                                          df.loc[df["BC"]<= bc].count().iloc[0]/df.count().iloc[0]*100,
                                                          da, df.loc[df["DA"] <= da].count().iloc[0],
                                                          df.loc[df["DA"] <= da].count().iloc[0]/df.count().iloc[
                                                              0]*100))



    def get_zfs_fig(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        df.rename(columns={"zfs_D": "D"}, inplace=True)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        c = "red"
        ax.axvspan(1, 2, facecolor=c, alpha=0.4)
        ax.text(1.5, 13, "L", verticalalignment="center", horizontalalignment="center")

        c="blue"
        ax.axvspan(2, 4, facecolor=c, alpha=0.4)
        ax.text(3, 13, "S", verticalalignment="center", horizontalalignment="center")

        c="green"
        ax.axvspan(4, 8, facecolor=c, alpha=0.4)
        ax.text(6, 13, "C", verticalalignment="center", horizontalalignment="center")

        c="orange"
        ax.axvspan(8, 12, facecolor=c, alpha=0.4)
        ax.text(10, 13, "X", verticalalignment="center", horizontalalignment="center")

        # ax.hist(df["D"], bins=25, edgecolor="black")
        df1 = df.loc[:, "D"]
        ax.bar(
            [i+0.5 for i in range(1, 12, 1)],
            [df1.loc[(df1 >= i) & (df1 < i + 1)].count() for i in range(1, 12, 1)],
            width=1,
            edgecolor="black",
        )

        # set text of count for each bin at the top of the bar and if the count is 0, then remove the text
        for i in range(1, 12, 1):
            count = df1.loc[(df1 >= i) & (df1 < i + 1)].count()
            if count == 0:
                ax.text(i+0.5, 0, "", verticalalignment="center", horizontalalignment="center")
            else:
                ax.text(i+0.5, count+1, str(count), verticalalignment="center", horizontalalignment="center")

        ax.set_xlim(1, 12)
        ax.set_ylim(0, 25)
        ax.set_xlabel("ZFS D (GHz)")
        ax.set_ylabel("Frequency")
        ax.set_title("ZFS parameter D distribution")
        # remove subticks
        ax.tick_params(axis="x", which="minor", bottom=False, top=False, labelbottom=False, direction="out")
        ax.tick_params(axis="x", which="major", bottom=True, top=False, direction="out")
        ax.tick_params(axis="y", which="major", left=True, right=False, labelbottom=True, direction="out")
        ax.tick_params(axis="y", which="minor", left=False, right=False, labelbottom=True, direction="out")

        self.figures.append(fig)
        return fig



    def get_zfs_stat(self, df=None):
        import numpy as np
        ## ref: http://www.sb.fsu.edu/~fajer/Fajerlab/LinkedDocuments/Ralph%20Weber.pdf
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        df1 = df.loc[:, "zfs_D"]
        minor_f = df1.loc[(df1 >= 0) & (df1 < 1)]
        L = df1.loc[(df1 >= 1) & (df1 < 2)]
        S = df1.loc[(df1 >= 2) & (df1 < 4)]
        C = df1.loc[(df1 >= 4) & (df1 < 8)]
        X = df1.loc[(df1 >= 8) & (df1 < 12)]
        Ku = df1.loc[(df1 >= 12) & (df1 < 18)]
        K = df1.loc[(df1 >= 18) & (df1 < 26)]
        Q = df1.loc[(df1 >= 26) & (df1 < 40)]
        V = df1.loc[(df1 >= 40) & (df1 < 75)]
        W = df1.loc[(df1 >= 75) & (df1 < 110)]
        D = df1.loc[(df1 >= 110) & (df1 < 170)]

        zfs_band_df = pd.DataFrame(
            {
                "name": ["minor_f", "L", "S", "C", "X", "Ku", "K", "Q", "V", "W", "D"],
                "frequency": [i.count() for i in [minor_f, L, S, C, X, Ku, K, Q, V, W, D]],
                "wavelength": ["0-1", "1-2", "2-4", "4-8", "8-12", "12-18", "18-26", "26-40", "40-75", "75-110", "110-170"],
            }
        )
        zfs_band_df.set_index("name", inplace=True)
        zfs_band_df["fraction"] = round(zfs_band_df.frequency / zfs_band_df.frequency.sum()*100, 0)
        display(zfs_band_df)
        zfs_band_df.loc[["L", "S", "C", "X"], "fraction"].sum()

        df1 = df.loc[:, "zfs_E"]
        #plot histogram of zfs_E
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(df1[df1<=1], bins=25, edgecolor="black")
        # ax.set_xlim(0, 1)
        # ax.set_ylim(0, 25)
        ax.set_xlabel("ZFS E")
        ax.set_ylabel("Frequency")
        ax.set_title("ZFS parameter E distribution")
        display(fig)




    def get_tdm_fig(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        # add inset plot for the distribution of the data in the main plot
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        axins = inset_axes(ax, width="30%", height="50%", loc="upper center")
        bin1 = df["TDM_rate"].loc[(df["TDM_rate"] >= 0) & (df["TDM_rate"] < 6)]
        bin2 = df["TDM_rate"].loc[(df["TDM_rate"] >= 6) & (df["TDM_rate"] <= 175)]
        axins.bar(
            [0,1],
            [bin1.count(), bin2.count()],
            width=0.5,
            edgecolor="black",
        )
        # set ticks as 0-6, 6-175
        axins.set_xticks([0, 1])
        axins.set_xticklabels(["[0-6)", "[6-175]"])
        axins.set_ylim(0, 30)
        # remove subticks
        axins.tick_params(axis="x", which="minor", bottom=False, top=False, labelbottom=False, direction="out")
        axins.tick_params(axis="x", which="major", bottom=False, top=False, direction="out")
        axins.tick_params(axis="y", which="both", left=True, right=False, labelbottom=True, direction="out")
        # set text of count for each bin at the top of the bar and if the count is 0, then remove the text
        for i in [0, 1]:
            count = bin1.count() if i == 0 else bin2.count()
            axins.text(i, count-1.8, str(count), verticalalignment="center", horizontalalignment="center")

        df1 = df.loc[:, "TDM_rate"]
        # plot bar chart for the distribution of the data in the main plot, with the range from 0 to 175 with increment of 5
        ax.bar(
            [i+2.5 for i in range(0, 175, 5)],
            [df1.loc[(df1 >= i) & (df1 < i + 5)].count() for i in range(0, 175, 5)],
            width=5,
            edgecolor="black",
        )
        ax.set_xlim(0, 175)

        # set text of count for each bin at the top of the bar and if the count is 0, then remove the text
        for i in range(len(ax.patches)):
            if ax.patches[i].get_height() == 0:
                ax.patches[i].set_visible(False)
            else:
                ax.text(
                    ax.patches[i].get_x() + ax.patches[i].get_width() / 2 + 2 if i < 2 else ax.patches[i].get_x() +
                                                                                             ax.patches[i].get_width() / 2,
                    ax.patches[i].get_height(),
                    str(int(ax.patches[i].get_height())),
                    ha="center",
                    va="bottom",
                )

        # ax.set_xlim(0, 175)
        # ax.set_ylim(0, 30)
        ax.set_xlabel("TDM rate (MHz)")
        ax.set_ylabel("Frequency")
        ax.set_title("TDM rate distribution")
        ax.tick_params(axis="x", which="minor", bottom=True, top=False, labelbottom=False, direction="out")
        ax.tick_params(axis="x", which="major", bottom=True, top=False, direction="out")
        ax.tick_params(axis="y", which="minor", left=False, right=False, labelbottom=True, direction="out")
        ax.tick_params(axis="y", which="major", left=True, right=False, labelbottom=True, direction="out")

        self.figures.append(fig)
        return fig

    def get_tdm_stat(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        df1 = df.loc[:, "TDM_rate"]
        bin1 = df1.loc[(df1 >= 0) & (df1 < 25)]
        bin2 = df1.loc[(df1 >= 25) & (df1 < 50)]
        bin3 = df1.loc[(df1 >= 50) & (df1 < 75)]
        bin4 = df1.loc[(df1 >= 75) & (df1 < 100)]
        bin5 = df1.loc[(df1 >= 100) & (df1 < 125)]
        bin6 = df1.loc[(df1 >= 125) & (df1 < 150)]
        bin7 = df1.loc[(df1 >= 150) & (df1 <= 175)]

        tdm_band_df = pd.DataFrame(
            {
                "bin": ["0-25", "25-50", "50-75", "75-100", "100-125", "125-150", "150-175"],
                "frequency": [bin1.count(), bin2.count(), bin3.count(),
                              bin4.count(), bin5.count(), bin6.count(), bin7.count()]
            }
        )
        tdm_band_df.set_index("bin", inplace=True)
        tdm_band_df["fraction"] = round(tdm_band_df.frequency / tdm_band_df.frequency.sum()*100, 0)
        display(tdm_band_df)

        sub_bin1 = df1.loc[(df1 >= 0) & (df1 < 6)]
        sub_bin2 = df1.loc[(df1 >= 6) & (df1 <= 175)]
        tdm_band_df = pd.DataFrame(
            {
                "bin": ["0-6", "6-175"],
                "frequency": [sub_bin1.count(), sub_bin2.count()],
            }
        )
        tdm_band_df.set_index("bin", inplace=True)
        tdm_band_df["fraction"] = round(tdm_band_df.frequency / tdm_band_df.frequency.sum()*100, 0)
        display(tdm_band_df)

        print(f"TDM rate >= 6: {df1.loc[df1 >= 6].count()}, fraction:"
              f" {round(df1.loc[df1 >= 6].count() / df1.count()*100, 2)}%")












    def test_valid_opt_tran(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        df["valid_zpl"] = df.apply(lambda x: x["ZPL"] <= x["vertical_transition_up"] if x["transition_from"] == "up"
        else x["ZPL"] <= x["vertical_transition_down"], axis=1)

        df["diff_zpl_vertical_transition"] = df.apply(lambda x: x["ZPL"] - x[
            "vertical_transition_up"] if x["transition_from"] == "up" else x["ZPL"] - x["vertical_transition_down"], axis=1)

        # validiation of AB>=0, BC>=0, CD>=0, DA>=0
        df["valid_cdft"] = df.apply(lambda x: x["valid_zpl"] and x["AB"] >= 0 and x["BC"] >= 0 and x["CD"] >= 0
                                              and x["DA"] >= 0, axis=1)

        return df.loc[df["valid_cdft"] == False, ["task_id",  "C2DB_uid", "ZPL", "LUMO-HOMO_up", "LUMO-HOMO_down",
                       "diff_zpl_vertical_transition", "AB", "BC", "CD", "DA", "transition_from", "up_tran_cdft_occ",
                       "dn_tran_cdft_occ", "dn_in_gap_deg", "up_in_gap_deg"]]


    @staticmethod
    def get_diff_btw_two_dfs(df1, df2):
        a = df1.copy().drop_duplicates(subset=["uid", "defect_name", "charge"], keep="last")
        b = df2.copy().drop_duplicates(subset=["uid", "defect_name", "charge"], keep="last")

        # get intersection of two dataframes and drop duplicates
        c = pd.merge(a, b, on=["uid", "defect_name", "charge"], how="inner")\
            .drop_duplicates(subset=["uid", "defect_name", "charge"], keep="last")

        # get a-c
        d = a.loc[~a.uid.isin(c.uid) | ~a.defect_name.isin(c.defect_name) | ~a.charge.isin(c.charge)]
        print(f"{len(d)} defects in df1 but not in df2")
        # get b-c
        e = b.loc[~b.uid.isin(c.uid) | ~b.defect_name.isin(c.defect_name) | ~b.charge.isin(c.charge)]
        print(f"{len(e)} defects in df2 but not in df1")

        # get count of df1 following c
        df1_intersection = df1.loc[df1.uid.isin(c.uid) & df1.defect_name.isin(c.defect_name) & df1.charge.isin(
            c.charge)].drop_duplicates(subset=["uid", "defect_name", "charge"])
        print(f"df1_intersection: {df1_intersection.shape}")

        # get count of df2 following c
        df2_intersection = df2.loc[df2.uid.isin(c.uid) & df2.defect_name.isin(c.defect_name) & df2.charge.isin(
            c.charge)].drop_duplicates(subset= ["uid", "defect_name", "charge"])
        print(f"df2_intersection: {df2_intersection.shape}")

        print(f"df1 intersection + (df1 - df1 intersection): {df1_intersection.shape[0] + d.shape[0]}")
        print(f"df2 intersection + (df2 - df2 intersection): {df2_intersection.shape[0] + e.shape[0]}")

        return d, e

#concate defect_qubit_df and defect_qubit_df_new based on task_id

class DefectAnalysis(SheetCollection):
    def __init__(self, df, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hse_candidate_df = df if df is not None else self.get_hse_candidate_df()
        self.ptmcs = self.hse_candidate_df.loc[
            (self.hse_candidate_df.prototype.isin(["GaS", "GaSe"])) &
            (self.hse_candidate_df.chemsys.isin(["Ga-S", "Ga-Se", "Ga-Te", "In-S", "In-Se", "In-Te"]))]

        self.xm = self.ptmcs.loc[self.ptmcs.charge == -1]
        self.mx = self.ptmcs.loc[self.ptmcs.charge == 1]

    def scp_files_from_server(self, taskid, files=None, db=None, base_path=None, open_file=False):
        base_path = base_path if base_path is not None else ""
        files = files if files else ["POSCAR"]
        db = db if db else HSEQubitDefect
        db_name = db.db.name
        col = db.collection
        col_name = col.name
        from subprocess import check_output
        for i in col.find({"task_id": taskid}):
            task_id = i["task_id"]
            dir_name = i["dir_name"].split("/")[-1]
            server_base_path = "/mnt/sdc/tsai/Research/projects"
            os.makedirs(os.path.join("structures",base_path, str(task_id)), exist_ok=True)
            for f in files:
                check_output(f"scp -P 12348 qimin@localhost:{server_base_path}/{db_name}/{col_name}/{dir_name}/"
                             f"{f}.gz structures/{base_path}/{task_id}/".split(" "))
                if f == "POSCAR":
                    check_output(f"gunzip structures/{base_path}/{taskid}/POSCAR.gz".split(" "), input=b"Y")
                    check_output(f"mv structures/{base_path}/{taskid}/POSCAR structures/{base_path}/{taskid}"
                                 f"/{taskid}.vasp".split(" "))
                    if open_file:
                        check_output(f"open structures/{base_path}/{taskid}/{taskid}.vasp".split(" "))


    def get_chg_from_wv(self, taskid, band, spin=0, base_path=None):
        from pymatgen.io.vasp.outputs import Wavecar
        from pymatgen.io.vasp.inputs import Poscar
        base_path = base_path if base_path else ""
        try:
            check_output(f"gunzip structures/{base_path}/{taskid}/WAVECAR.gz".split(" "), input=b"Y")
            check_output(f"gunzip structures/{base_path}/{taskid}/POSCAR.gz".split(" "), input=b"Y")
        except Exception:
            pass
        wv = Wavecar(f"structures/{base_path}/{taskid}/WAVECAR")
        wv.get_parchg(Poscar.from_file(f"structures/{base_path}/{taskid}/{taskid}.vasp"), 0, band-1, spin,
                      phase=True).write_file(
            f"structures/{base_path}/{taskid}/{taskid}-{band}_spin-{spin}.vasp")
        check_output(f"open structures/{base_path}/{taskid}/{taskid}-{band}_spin-{spin}.vasp".split(" "))


    def get_local_env(self):

        # if self.xm doesn't have column NN, self.get_NN() will be called to get NN for self.xm
        self.xm = self.get_NN(self.xm) if not "NN" in list(self.xm.columns) else self.xm

        self.xm_GaS = self.xm.loc[self.xm.prototype == "GaS"]
        self.xm_GaSe = self.xm.loc[self.xm.prototype == "GaSe"]
        def plot_dz_vs_source_specie():
            # plot dz vs source specie for xm
            fig, ax = plt.subplots()
            # group by chemsys
            for i, (name, group) in enumerate(self.xm_GaSe.groupby("chemsys")):
                # plot line for each group
                ax.plot(group.level_source_specie, group.dn_tran_en, label=name, marker="o", linestyle="-")
            ax.set_xlabel("source specie")
            ax.set_ylabel("dz")
            ax.set_title(r"\mathrm{dz} vs \mathrm{source specie} in $\mathrm{X_{M}^{1-}}$")
            # set legend outside of plot
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            return fig


        plot_dz_vs_source_specie()

    def get_proj_analysis(self, df=None):
        df = df if df is not None else self.xm_GaS
        self.proj_analysis = []
        for taskid in df["task_id"]:
            e_df = df.loc[df["task_id"] == taskid]
            nn = e_df.NN.values[0]
            print(nn)
            a, b, c, d = Tools.extract_defect_levels_v2_hse(taskid,
                localisation=float(e_df["localisation_threshold"].iloc[0]), tol=(0.5, 0.5)
        )

            d["vertical_species_proj"] = d.apply(lambda x: x[nn[0]] + x[nn[-1]], axis=1)
            d["horizontal_species_proj"] = d.apply(lambda x: sum([x[i] for i in nn[1:4]]), axis=1)
            d.reset_index(drop=False, inplace=True)
            d = d.loc[:, ["band_index", "spin", "orbital", "vertical_species_proj", "horizontal_species_proj"]]
            # fill a column with task_id
            d["task_id"] = pd.Series([taskid] * len(d))
            # create a new dataframe with all d
            self.proj_analysis.append(d)
        self.proj_analysis = pd.concat(self.proj_analysis, axis=0)
        self.proj_analysis.reset_index(drop=True, inplace=True)

        # plot horizontal and vertical projection


    def get_geo_NN(self, df=None):
        df = df if df is not None else self.ptmcs
        from pymatgen import PeriodicSite
        dist_df = {"task_id": [], "a":[], "b":[], "gamma":[], "NN_dist": [], "NN_index": [], "distinct_NN_dist": [],
                   "host_NN_dist":[], "host_NN_index":[], "distinct_host_NN_dist":[], "gap_hse_line":[],
                   "d1 in host": [], "d2 in host": [], "d1 in defect": [], "d2 in defect": [], "C2DB_uid": [],
                   "defect_name": [], "charge": []}
        # electronic_df = {"bandgap": []}

        for taskid in df.loc[:, "task_id"].tolist():
            dist_df["task_id"].append(taskid)
            dist_df["C2DB_uid"].append(df.loc[df["task_id"] == taskid, "C2DB_uid"].values[0])
            dist_df["defect_name"].append(df.loc[df["task_id"] == taskid, "defect_name"].values[0])
            dist_df["charge"].append(df.loc[df["task_id"] == taskid, "charge"].values[0])

            e = HSEQubitDefect.collection.find_one({"task_id": taskid})

            bulk_site = PeriodicSite.from_dict(e["defect_entry"]["unique_site"])
            host_st = Structure.from_dict(e["defect_entry"]["supercell"]["bulk"])

            NN_in_bulk = host_st.get_neighbors(bulk_site, 3, include_index=True, include_image=False)
            NN_in_bulk_sites = [site.index for site in NN_in_bulk]
            host_NN_dist = np.array([site.nn_distance for site in NN_in_bulk]).round(3)
            dist_df["host_NN_dist"].append(list(host_NN_dist))

            distinct_host_NN_dist, counts = np.unique(host_NN_dist, return_counts=True)
            sort_index = np.argsort(counts)
            distinct_host_NN_dist = distinct_host_NN_dist[sort_index]
            distinct_host_NN_dist = distinct_host_NN_dist[-2:] # take the last two
            dist_df["distinct_host_NN_dist"].append(list(distinct_host_NN_dist))
            dist_df["d1 in host"].append(distinct_host_NN_dist[0])
            dist_df["d2 in host"].append(distinct_host_NN_dist[1] if len(distinct_host_NN_dist) > 1 else "None")
            dist_df["host_NN_index"].append(NN_in_bulk_sites)
            dist_df["a"].append(host_st.lattice.a/e["defect_entry"]["supercell"]["size"][0],)
            dist_df["b"].append(host_st.lattice.b/e["defect_entry"]["supercell"]["size"][1],)
            dist_df["gamma"].append(host_st.lattice.gamma,)

            dist_df["gap_hse_line"].append(e["host_info"]["c2db_ir_hse_line"]["bandgap"])

            # defect
            nn = e["NN"]
            st = Structure.from_dict(e["output"]["structure"])
            NN_in_defect = st.get_neighbors(st.sites[nn[-1]], 3.5, include_index=True, include_image=False)
            NN_in_defect_sites = [site.index for site in NN_in_defect]

            NN_dist = np.array([i.nn_distance for i in NN_in_defect]).round(3)
            dist_df["NN_dist"].append(list(NN_dist))
            # find distinct entries in distinct_NN_dist, and sort them by counts in distinct_NN_dist
            distinct_NN_dist, counts = np.unique(NN_dist, return_counts=True)
            # sort the NN_dist by counts
            sort_index = np.argsort(counts)
            # sort the NN_dist by sort_index
            distinct_NN_dist = distinct_NN_dist[sort_index]
            dist_df["distinct_NN_dist"].append(list(distinct_NN_dist))
            dist_df["d1 in defect"].append(distinct_NN_dist[0])
            dist_df["d2 in defect"].append(distinct_NN_dist[1] if len(distinct_NN_dist) > 1 else "None")
            dist_df["NN_index"].append(list(NN_in_defect_sites))

        dist_df = pd.DataFrame(dist_df)
        # replace the rows of columns "d1 in host", "d2 in host", "distinct_host_NN_dist", "host_NN_dist,
        # and "NN_index" with the last rows that have the same C2DB_uid
        good_host = dist_df.groupby(["C2DB_uid"]).last().reset_index()
        good_host = good_host.loc[:, ["C2DB_uid", "d1 in host", "d2 in host", "distinct_host_NN_dist",
                                       "host_NN_dist", "host_NN_index"]]
        display(good_host)
        # if C2DB_uid in dist_df is same as C2DB_uid in good_host, replace the columns with those in good_host
        dist_df.loc[:, ["d1 in host", "d2 in host", "distinct_host_NN_dist", "host_NN_dist", "host_NN_index"]] = \
            dist_df.loc[:, ["C2DB_uid"]].merge(good_host, on="C2DB_uid", how="left")


        # self.ptmcs = pd.merge(self.ptmcs, dist_df, on="task_id", how="left").round(3)
        return dist_df.round(3)

    def get_host_materials_info(self, df):
        def gap_hse_line_plt():
            fig, ax = plt.subplots()
            GaS_df = df.loc[df["prototype"] == "GaS"]
            GaSe_df = df.loc[df["prototype"] == "GaSe"]
            ax.scatter(GaSe_df["chemsys"].tolist()[::-1], GaS_df["gap_hse_line"].tolist()[::-1], label="GaS-hse-line",)
            ax.scatter(GaSe_df["chemsys"].tolist()[::-1], GaSe_df["gap_hse_line"].tolist()[::-1], label="GaSe-hse-line")
            ax.set_xlabel("Chemical system")
            ax.set_ylabel("Bandgap (eV)")
            ax.set_title("Bandgap from HSE-line/ Lattice constant a")
            # set legend position at the upper left of the plot
            ax.legend(loc="upper left", bbox_to_anchor=(-0.4, 1))

            ax2 = ax.twinx()
            ax2.scatter(GaSe_df["chemsys"].tolist()[::-1], GaS_df["a"].tolist()[::-1], label="a", color="red")
            ax2.set_ylabel("a (Å)")
            ax2.legend(loc="upper right", bbox_to_anchor=(1.2, 1))
            return fig

        return gap_hse_line_plt()


    def get_defects_info(self):
        def singlet_triplet_en_diff():
            triplets = self.ptmcs["task_id"].unique()
            singlets = {"task_id": [], "singlet_taskid": [], "C2DB_uid": [], "defect_name": [], "charge":[],
                        "singlet_energy": [], "triplet_energy": [], "singlet_triplet_energy_diff": [],}
            for t in triplets:
                print(t)
                triplet_e = HSEQubitDefect.collection.find_one({"task_id": int(t)})
                singlet_e = HSEQubitDefect.collection.find_one({"pc_from": triplet_e["pc_from"], "nupdown_set": 0,
                                                                "defect_entry.name": triplet_e["defect_entry"]["name"],
                                                                "charge_state": triplet_e["charge_state"]})

                singlets["singlet_taskid"].append(singlet_e["task_id"])
                singlets["task_id"].append(int(t))
                singlets["singlet_energy"].append(singlet_e["output"]["energy"])

                singlets["C2DB_uid"].append(self.ptmcs.loc[self.ptmcs["task_id"] == t]["C2DB_uid"].tolist()[0])
                singlets["defect_name"].append(self.ptmcs.loc[self.ptmcs["task_id"] == t]["defect_name"].tolist()[0])
                singlets["charge"].append(self.ptmcs.loc[self.ptmcs["task_id"] == t]["charge"].tolist()[0])
                singlets["triplet_energy"].append(triplet_e["output"]["energy"])
                singlets["singlet_triplet_energy_diff"].append(singlet_e["output"]["energy"] - triplet_e["output"]["energy"])
            return pd.DataFrame(singlets).round(3)
        return singlet_triplet_en_diff()




if __name__ == '__main__':
    # df = pd.read_csv("/home/yluo/Documents/defect_qubit_df.csv")
    #
    pass

