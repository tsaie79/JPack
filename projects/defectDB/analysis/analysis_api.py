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
from pymatgen.core.units import ang_to_bohr

from pathlib import Path

from JPack_independent.projects.defectDB.analysis.data_analysis import *
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
# hse_qubit_df = IOTools(cwd=input_path, excel_file="hse_qubit_2022-04-21").read_excel()
hse_qubit_df = IOTools(cwd=input_path, excel_file="hse_qubit_2022-07-10").read_excel()

tmd_antisites_df = IOTools(cwd=input_path, excel_file="tmd_antisites_2022-04-13").read_excel()
pure_c2db_pc_tmd_antisite_df = IOTools(cwd=input_path, excel_file="pure_c2db_TMD_2022-07-08").read_excel()

hse_screened_qubit_df = IOTools(cwd=input_path, excel_file="hse_screened_qubit_2021-11-18").read_excel()
hse_screened_qubit_cdft_df = IOTools(cwd=input_path, excel_file="hse_screened_qubit_cdft_2021-11-18_2022-01-31").read_excel()
hse_screened_qubit_cdft_zpl_df = IOTools(cwd=input_path, excel_file="hse_screened_qubit_cdft_zpl_2021-11-18_2022-01-31").read_excel()

# defects_36_groups_df = IOTools(cwd=input_path, excel_file="defects_36_groups_2022-03-07").read_excel()
defects_36_groups_df = IOTools(cwd=input_path, excel_file="defects_36_groups_2022-04-27").read_excel()


defects_36_groups_zpl_zfs_df = IOTools(cwd=input_path, excel_file="defects_36_groups_zpl_zfs").read_excel()

# hse_candidate_df = IOTools(cwd=input_path, excel_file="Table_5_df_2022-04-15").read_excel()
hse_candidate_df = IOTools(cwd=input_path, excel_file="Table_5_df_2022-12-05").read_excel()

prev_hse_candidate_df = IOTools(cwd=input_path, excel_file="Table_5_df_2022-06-01").read_excel()

table1_df = IOTools(cwd=input_path, excel_file="Table_1_df_2022-04-15").read_excel()
table4_df = IOTools(cwd=input_path, excel_file="Table_4_df_2022-12-05").read_excel()

new_defect_df = IOTools(cwd=input_path, excel_file="new_defect_2022-10-14").read_excel()

# os.environ["PATH"] += os.pathsep + "/Library/TeX/texbin/latex"
# plt.style.use(["science", "nature", "light"])

class Toolbox:
    @classmethod
    def get_poscar_and_wavecar(cls, taskid, hse, files=["POSCAR"]):
        from qubitPack.tool_box import get_db
        if hse:
            col = get_db("HSE_triplets_from_Scan2dDefect", "calc_data-pbe_pc", port=12349).collection
        else:
            col = get_db("Scan2dDefect", "calc_data", port=12349).collection
        def scp(taskid):
            for i in col.find({"task_id": taskid}):
                task_id = i["task_id"]
                dir_name = i["dir_name"].split("/")[-1]
                if hse:
                    src_dir = "/home/qimin/sdc_tsai/Research/projects/HSE_triplets_from_Scan2dDefect/calc_data-pbe_pc/"
                    tgt_dir = str(task_id)
                else:
                    src_dir = "/home/qimin/sdc_tsai/Research/projects/Scan2dDefect/calc_data/scf"
                    tgt_dir = str(task_id)+"_SCAN"
                os.makedirs(os.path.join("structures", tgt_dir), exist_ok=True)
                for f in files:
                    f_out = f"{taskid}_{f}.gz"
                    check_output("scp -P 12348 qimin@localhost:{}/{}/{}.gz structures/{}/{}".format(
                        src_dir,
                        dir_name, f, tgt_dir, f_out).split(" "))
        scp(taskid)
    @classmethod
    def get_parcharg(cls, taskid, band, spin, hse):
        try:
            if hse:
                tgt_dir = str(taskid)
            else:
                tgt_dir = str(taskid)+"_SCAN"
            check_output("gunzip structures/{}/{}_WAVECAR.gz".format(tgt_dir, taskid).split(" "),
                         input=b"Y")
            check_output("gunzip structures/{}/{}_POSCAR.gz".format(tgt_dir, taskid).split(" "), input=b"Y")
        except Exception:
            pass
        wv = Wavecar("structures/{}/{}_WAVECAR".format(tgt_dir, taskid))
        wv.get_parchg(Poscar.from_file("structures/{}/{}_POSCAR".format(tgt_dir, taskid)), 0, band-1, spin, phase=True). \
            write_file("structures/{}/{}-{}_spin-{}.vasp".format(tgt_dir, taskid, band, spin))
        check_output("open structures/{}/{}-{}_spin-{}.vasp".format(tgt_dir, taskid, band, spin).split(" "))
        check_output("gzip structures/{}/{}_WAVECAR".format(tgt_dir, taskid).split(" "))
        check_output("gzip structures/{}/{}_POSCAR".format(tgt_dir, taskid).split(" "))

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
        scan_cohp_db = get_db("Scan2dDefect", "lobster", user="Jeng_ro", password="qimin", port=12349)
        scan_mat_db = get_db("Scan2dMat", "calc_data", user="Jeng_ro", password="qimin", port=12349)


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
            zfs_collection = get_db("HSE_triplets_from_Scan2dDefect", "zfs_data-pbe_pc", port=12349).collection
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
            "prototype",
            "spacegroup", "C2DB_uid", "good_site_sym", "good_site_specie",
            "defect_type", "charge", "perturbed_vbm", "perturbed_cbm",
            "perturbed_bandgap", "level_edge_category",
            "LUMO-HOMO_up",  "LUMO-HOMO_down"
        ]

        self.grouping_hse_candidate = [
            "prototype", "spacegroup", "C2DB_uid", "defect_type", "defect_name", "charge", "defect_symmetry",
            "perturbed_bandgap","zfs_D", "zfs_E",  "LUMO-HOMO_up", "LUMO-HOMO_down", "LUMO_HOMO_deg",
            "LUMO_HOMO_splitting_type",
            "transition_from",
            "gap_tran_bottom_to_vbm", "gap_tran_top_to_cbm",
            "ZPL", "AB", "BC", "CD", 'DA',
            "up_TDM_rate_X", "up_TDM_rate_Y", "up_TDM_rate_Z", "dn_TDM_rate_X", "dn_TDM_rate_Y", "dn_TDM_rate_Z",
            "TDM_rate", "up_polarization", "dn_polarization", "allowed", "class", "singlet_triplet_energy_diff",
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

    def get_scan_defect_df(self, df=None, include_ipr=True):
        # defect_df is generated by data_preparation.py
        df = defect_df.copy() if df is None else df
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
        if include_ipr:
            da = DefectAnalysis(df)
            da.get_IPR(df)
            self.scan_defect_df = df.merge(da.ipr_df, on=["task_id"], how="left")
        else:
            self.scan_defect_df = df
        # use {"task_id": "unique"} to make sure there is no duplicated defects, then use {"task_id": "max"}
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
        df = self.check_vertical_transition(df, 0.05)
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

        #add singlet informaion
        da = DefectAnalysis(df)
        da.get_singlet_triplet_en_diff()
        df = df.merge(da.singlet_df, on=["C2DB_uid", "defect_name", "charge", "task_id"], how="left")
        # compliment TMDs' results
        for tmd_taskid in tmd_antisites_df["task_id"].unique():
            chemsys = tmd_antisites_df.loc[tmd_antisites_df["task_id"] == tmd_taskid, "chemsys"].iloc[0]
            if chemsys == "S-W":
                df.loc[(df["task_id"] == tmd_taskid) & (df["chemsys"] == "S-W"), "singlet_triplet_energy_diff"] = \
                    0.273
            elif chemsys == "Se-W":
                df.loc[(df["task_id"] == tmd_taskid) & (df["chemsys"] == "Se-W"), "singlet_triplet_energy_diff"] = \
                    0.295
            elif chemsys == "Te-W":
                df.loc[(df["task_id"] == tmd_taskid) & (df["chemsys"] == "Te-W"), "singlet_triplet_energy_diff"] = \
                    0.301
            elif chemsys == "Mo-S":
                df.loc[(df["task_id"] == tmd_taskid) & (df["chemsys"] == "Mo-S"), "singlet_triplet_energy_diff"] = \
                    0.349
            elif chemsys == "Mo-Se":
                df.loc[(df["task_id"] == tmd_taskid) & (df["chemsys"] == "Mo-Se"), "singlet_triplet_energy_diff"] = \
                    0.343
            elif chemsys == "Mo-Te":
                df.loc[(df["task_id"] == tmd_taskid) & (df["chemsys"] == "Mo-Te"), "singlet_triplet_energy_diff"] = \
                    0.264

        df = df.loc[df["singlet_triplet_energy_diff"].notnull()]
        df = df.loc[df["singlet_triplet_energy_diff"] > 0]

        #add defect symmetry
        da.get_defect_symmetry(df)
        df = df.merge(da.defect_symmetry_df, on=["task_id"], how="left")

        self.hse_triplet_df = df
        grouping = self.grouping.copy()
        grouping.append("singlet_triplet_energy_diff")
        grouping.append("defect_symmetry")
        self.hse_triplet_df_gp = df.groupby(grouping).agg({"task_id": ["unique", "count"]})

        return self.hse_triplet_df


    def SP_get_hse_triplet_df(self, df=None, remove_duplicates_based_on=["uid", "defect_name", "charge"],
                           remove_taskid_in_36group=[186, 1020, 1000, 283, 211, 644, 229, 548, 242, 603]):
        """
        Get HSE triplet dataframe WITHOUT MATCHING THE TMD ANTISITES FROM OUR NATURE COMM. PAPER.
        Therefore, ALL pc are from C2DB.
        :param df:
        :param remove_duplicates_based_on:
        :param remove_taskid_in_36group:
        :return:
        """
        df = df if df is not None else hse_qubit_df.copy()
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

        #add singlet informaion
        # da = DefectAnalysis(df)
        # da.get_singlet_triplet_en_diff()
        # df = df.merge(da.singlet_df, on=["C2DB_uid", "defect_name", "charge", "task_id"], how="left")

        df = df.loc[df["singlet_triplet_energy_diff"] != "None"]
        df = df.loc[df["singlet_triplet_energy_diff"] > 0]

        #add defect symmetry
        # da = DefectAnalysis(df)
        # da.get_defect_symmetry(df)
        # df = df.merge(da.defect_symmetry_df, on=["task_id"], how="left")

        self.hse_triplet_df = df
        grouping = self.grouping.copy()
        grouping.append("singlet_triplet_energy_diff")
        grouping.append("defect_symmetry")
        self.hse_triplet_df_gp = df.groupby(grouping).agg({"task_id": ["unique", "count"]})

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
            gap_tran_bottom_to_vbm, gap_tran_top_to_cbm = None, None
            if x["transition_from"] == "dn":
                gap_tran_bottom_to_vbm = x["dn_tran_bottom_to_vbm"]
                gap_tran_top_to_cbm = x["dn_tran_top_to_cbm"]
            elif x["transition_from"] == "up":
                gap_tran_bottom_to_vbm = x["up_tran_bottom_to_vbm"]
                gap_tran_top_to_cbm = x["up_tran_top_to_cbm"]
            return gap_tran_bottom_to_vbm, gap_tran_top_to_cbm

        df["gap_tran_bottom_to_vbm"] = df.apply(lambda x: fun(x)[0], axis=1)
        df["gap_tran_top_to_cbm"] = df.apply(lambda x: fun(x)[1], axis=1)

        # df["gap_tran_bottom_to_vbm"] = df.apply(
        #     lambda x: x["dn_tran_bottom"] if x["transition_from"] == "dn" else x["up_tran_bottom"], axis=1)
        # df["gap_tran_top_to_cbm"] = df.apply(
        #     lambda x: x["level_gap"] - x["dn_tran_top"] if x["transition_from"] == "dn" else x["level_gap"] - x["up_tran_top"], axis=1)

        def get_filtered_df(df, D_zfs=0.95, wavelength_zpl=2500, homo_to_vbm=0.050, lumo_to_cbm=0.095):
            """
            :param df:
            :param D_zfs: 0.95
            :param wavelength_zpl: 2500
            :param homo_to_vbm: 0.095
            :param lumo_to_cbm: 0.095
            :return:
            """
            df = df.loc[df["zfs_D"] != "None"]
            df = df.loc[df["zfs_D"] >= D_zfs]

            df = df.loc[df["allowed"] == True]
            df = df.loc[df["ZPL"] != "None"]
            df = df.loc[df["ZPL_wavelength"] <= wavelength_zpl]

            df = df.loc[df["gap_tran_bottom_to_vbm"] >= homo_to_vbm]
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
        self.hse_candidate_df["LUMO_HOMO_splitting_type"] = self.hse_candidate_df.apply(
            lambda x: f"{x['LUMO_HOMO_deg'][1]}-{x['LUMO_HOMO_deg'][0]}_type",axis=1)

        self.hse_candidate_df_gp = self.hse_candidate_df.copy().replace(np.nan, "None").groupby(
            self.grouping_hse_candidate).agg({"task_id": "max"})
        self.hse_candidate_not_passed_df = not_passed_df.replace(np.nan, "None").copy()
        not_passed_grouping = self.grouping_hse_candidate.copy()
        not_passed_grouping.remove("LUMO_HOMO_splitting_type")
        self.hse_candidate_not_passed_df_gp = self.hse_candidate_not_passed_df.copy().replace(np.nan, "None").groupby(
            not_passed_grouping).agg(
            {"task_id": "unique"})

    def get_final_hse_candidate_df(self, hse_candidate_df=None, remove_e_config_summary_cols=True):
        df = hse_candidate_df.copy()
        # remove e_config, h_config, and e2_from. This must be done if get_e_config_summary_df was utilized
        if remove_e_config_summary_cols:
            df = df.drop(["e_config", "h_config", "e2_from"], axis=1)
        # get_e_config
        da = DefectAnalysis(df)
        updated_df = da.get_e_config_summary_df() # these e_config are explicitly identified by eyes

        # screen out some candidates
        updated_df = updated_df.loc[
            (updated_df["e_config"]!="unknown") &
            (updated_df["h_config"]!="unknown") &
            (updated_df["e2_from"]!="unknown")
            ]
        updated_df = updated_df.loc[updated_df["singlet_triplet_energy_diff"]!="None"]

        # remove D3h 109, 129 since they don't have proper intersystem-crossings for qubit operations.
        updated_df = updated_df.loc[~(updated_df["task_id"].isin([109, 129]))]

        #remove D3, task_id == 107 for complex defect levels
        updated_df = updated_df.loc[~(updated_df["task_id"].isin([107]))]

        #remove strange structure task_id == 1234
        updated_df = updated_df.loc[~(updated_df["task_id"].isin([1234]))]

        self.updated_hse_candidate_df = updated_df.copy()
        grouping = self.grouping_hse_candidate.copy()
        grouping.extend(["e_config", "h_config", "e2_from"])
        self.updated_hse_candidate_df_gp = self.updated_hse_candidate_df.groupby(grouping).agg({"task_id": "unique"})




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

    def get_plain_grand_defect_level_plot(self, df, sym_name, ylim):
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
        ax.set_ylim(vbm_lim-fig_y_height[0], cbm_lim+fig_y_height[1])

        print("len of defects: {}".format(len(defects)))
        # plot mutilple bars
        for defect, x in zip(defects, np.arange(0, len(defects)*5, 5)):
            # readout all the data from defect dict
            vbm = defect["level_vbm"]
            cbm = defect["level_cbm"]
            up = defect["up_in_gap_level"]
            dn = defect["dn_in_gap_level"]
            up_deg = defect["up_in_gap_deg"]
            dn_deg = defect["dn_in_gap_deg"]
            up_occ = defect["up_in_gap_occ"]
            dn_occ = defect["dn_in_gap_occ"]
            up_band = defect["up_in_gap_band_index"]
            dn_band = defect["dn_in_gap_band_index"]

            if not up_deg:
                up_deg = ["" for i in up]
            if not up_occ:
                up_occ = ["" for i in up]

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
                # ax.text(x-2, level, id)
                color = None
                if deg == 1:
                    color = "black"
                    ax.hlines(level, x-2, x-1, colors=color)
                    if occ == 1:
                        # plot upward arrow
                        ax.arrow(x-1.5, level-0.05, 0, 0.1, head_width=0.5, head_length=0.1, color=color)

                elif deg == 2:
                    color = "black"
                    ax.hlines(level, x-2, x-1, colors=color)
                    ax.hlines(level, x-0.5, x+0.5, colors=color)
                    if occ == 1:
                        # plot upward arrow
                        ax.arrow(x-1.5, level-0.05, 0, 0.1, head_width=0.5, head_length=0.1, color=color)
                        ax.arrow(x, level-0.05, 0, 0.1, head_width=0.5, head_length=0.1, color=color)
                elif deg == 3:
                    color = "green"
                elif deg == 4:
                    color = "yellow"
                else:
                    color = "black"


            dx = 1
            for level, occ, deg, id in zip(dn, dn_occ, dn_deg, dn_band):
                # ax.text(x+2, level, id)
                color = None
                if deg == 1:
                    color = "black"
                    ax.hlines(level, x+1, x+2, colors=color)
                    if occ == 1:
                        # plot downward arrow
                        ax.arrow(x+1.5, level+0.05, 0, -0.1, head_width=0.5, head_length=0.1, color=color)
                elif deg == 2:
                    color = "black"
                    ax.hlines(level, x-0.5, x+0.5, colors=color)
                    ax.hlines(level, x+1, x+2, colors=color)
                    if occ == 1:
                        ax.arrow(x, level+0.05, 0, -0.1, head_width=0.5, head_length=0.1, color=color)
                        ax.arrow(x+1.5, level+0.05, 0, -0.1, head_width=0.5, head_length=0.1, color=color)
                elif deg == 3:
                    color = "green"
                elif deg == 4:
                    color = "yellow"
                else:
                    color = "black"


        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        for tick in ax.get_yticklabels():
            tick.set_fontsize(15)

        ax.tick_params(axis="x", bottom=False, labelbottom=True)
        ax.tick_params(axis="y", right=False, which="both")


        ax.set_xticks(np.arange(0, len(defects)*5, 5))

        ax.set_xticklabels([rf"$\mathrm{{{defect_name_labels[i].split('_')[2]}_{{"
                            rf"{defect_name_labels[i].split('_')[4]}}}^{{{charge_labels[i]}}}}}$"+
                            rf"@{''.join(host_labels[i].split('-')[0].split('2'))}" for
                            i in range(len(mags))], fontdict={"fontsize": 20}, rotation=270)

        # remove minor and major xticks
        ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=True)
        print(f"taskid_label: {taskid_labels}, mags: {mags}")

        ax.set_ylabel("Energy relative to vacuum (eV)", fontsize=20)
        ax.set_ylim(ylim[0], ylim[1])
        return fig, ax


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
        # for i in range(250, 2000, 50):
        #     if df1.loc[(df1 >= i) & (df1 < i + 50)].count() != 0:
        #         ax.text(i + 25, df1.loc[(df1 >= i) & (df1 < i + 50)].count(),
        #                 str(df1.loc[(df1 >= i) & (df1 < i + 50)].count()),
        #                 horizontalalignment="center",
        #                 verticalalignment="bottom")


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
        # for i in range(len(df1.value_counts(bins=np.arange(0, 0.4, 0.01)).sort_index())):
        #     if df1.value_counts(bins=np.arange(0, 0.4, 0.01)).sort_index().iloc[i] != 0:
        #         ax.text(
        #             np.arange(0, 0.4-0.01, 0.01)[i]+0.01/2,
        #             df1.value_counts(bins=np.arange(0, 0.4, 0.01)).sort_index().iloc[i],
        #             df1.value_counts(bins=np.arange(0, 0.4, 0.01)).sort_index().iloc[i],
        #             horizontalalignment="center",
        #             verticalalignment="bottom",
        #             )

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
        # for i in range(len(df1.value_counts(bins=np.arange(0, 0.2, 0.01)).sort_index())):
        #     if df1.value_counts(bins=np.arange(0, 0.2, 0.01)).sort_index().iloc[i] != 0:
        #         axins.text(
        #             np.arange(0, 0.2-0.01, 0.01)[i]+0.01 if i < 2 else np.arange(0, 0.2-0.01, 0.01)[i]+0.01/2,
        #             df1.value_counts(bins=np.arange(0, 0.2, 0.01)).sort_index().iloc[i],
        #             df1.value_counts(bins=np.arange(0, 0.2, 0.01)).sort_index().iloc[i],
        #             horizontalalignment="center",
        #             verticalalignment="bottom",
        #         )

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
        print("BC<={}, N={},{}%, DA<={}, N={},{}%".format(bc, df.loc[df["BC"] <= bc].count().iloc[0],
                                                          df.loc[df["BC"]<= bc].count().iloc[0]/df.count().iloc[0]*100,
                                                          da, df.loc[df["DA"] <= da].count().iloc[0],
                                                          df.loc[df["DA"] <= da].count().iloc[0]/df.count().iloc[
                                                              0]*100))



    def get_zfs_fig(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()
        df.rename(columns={"zfs_D": "D"}, inplace=True)
        df = df.groupby("task_id").agg({"D": "max"}).reset_index()

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
        # for i in range(1, 12, 1):
        #     count = df1.loc[(df1 >= i) & (df1 < i + 1)].count()
        #     if count == 0:
        #         ax.text(i+0.5, 0, "", verticalalignment="center", horizontalalignment="center")
        #     else:
        #         ax.text(i+0.5, count+1, str(count), verticalalignment="center", horizontalalignment="center")

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
        df1 = df.groupby("task_id").agg({"zfs_D": "max"}).reset_index()["zfs_D"]
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

    def get_zfs_vs_antisite_specie(self, df=None):
        if df is None:
            df = self.get_hse_candidate_df().copy() if df is None else df.copy()
            da = DefectAnalysis(df)
            df_antisite_proj = da.get_antisite_projection()
        else:
            df_antisite_proj = df.copy()
        df_antisite_proj["group_row"] = df_antisite_proj.apply(
            lambda x: f'{x["level_source_group"]}-{x["level_source_row"]}', axis=1)
        x = df_antisite_proj.groupby(["level_source_group"])
        fig, ax = plt.subplots()
        for k, v in x:
            print(k)
            r = None
            if k in [6]:
                r = "average_cationic_radius"
            elif k in [16, 17]:
                r = "Atomic radius calculated"
            else:
                r = "Atomic radius calculated"
            v = v.sort_values("level_source_row")
            display(v.loc[:, ["level_source_specie", r, "zfs_D"]])
            label = {6: "VIB", 13: "IIIA", 16: "VIA", 17: "VIIA", 12: "IIB", 14:"IVA"}
            ax.scatter(v[r], v["zfs_D"], label=label[k])
            # display(v.loc[:, ["group_row", "zfs_D"]])
            # ax.scatter(v["group_row"], v["zfs_D"])
        # set y axis range from 0 to 12
        ax.set_ylim(0, 12)
        # ax.set_xlim(0.8, 2)
        # rotate legend by 90 degree
        ax.legend(loc="upper center", ncol=4, frameon=True)
        ax.set_ylabel("ZFS D (GHz)")
        ax.set_xlabel("Atomic radius (Å)")
        ax.set_title("ZFS D vs. atomic radius of antisite atom")

        #remove subticks in x-axis
        ax.tick_params(axis="x", which="minor", bottom=True, top=False, labelbottom=False, direction="out")
        ax.tick_params(axis="x", which="major", bottom=True, top=False, labelbottom=True, direction="out")
        ax.tick_params(axis="y", which="minor", left=True, right=False, labelbottom=False, direction="out")
        ax.tick_params(axis="y", which="major", left=True, right=False, labelbottom=True, direction="out")

        #rotate x-axis labels by 90 degrees
        # plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        return fig

    def get_zfs_vs_ipr_fig(self, df, turn_on_lr=True):
        from scipy import stats
        plot_df = df.copy()

        x_name = "Spin-up HODL IPR"
        plot_df.rename(columns={"up_homo_ipr": x_name}, inplace=True)
        plot_df.rename(columns={"zfs_D": "ZFS D"}, inplace=True)
        group_df = plot_df.groupby(["chem_group"])
        fig = plt.figure(figsize=(12, 3))
        # use first 3 rows and first 3 columns as a big plot
        ax1 = plt.subplot2grid((3, 26), (0, 0), colspan=10, rowspan=3)
        colors = ["blue", "red", "lime", "purple", "orange", "cyan"]
        group_name = {6: "VIB", 12: "IIB", 13: "IIIA", 14: "IVA", 15: "VA", 16: "VIA", 17: "VIIA"}
        plot_markers = ['(b)', '(c)', '(d)', '(e)', '(f)', '(g)']
        for (name, group), index, color, plot_mark in zip(group_df, range(len(group_df)), colors, plot_markers):
            name = group_name[name]
            ax1.plot(group[x_name], group["ZFS D"], marker="o", linestyle="", label=name,
                     color=color)

            x = group[x_name]
            y = group["ZFS D"]

            # plot scatter plot and regression line for each group at those subplots at last two columns and first 3 rows
            subplot_locs = [(0, 10), (1, 10), (2, 10), (0, 18), (1, 18), (2, 18)]
            ax2 = plt.subplot2grid((3, 26), (subplot_locs[index]), colspan=8, rowspan=1)
            # make text box at middle top to show name
            ax2.text(0.5, 1.03, name, horizontalalignment="center", verticalalignment="center", transform=ax2.transAxes)
            ax2.plot(x, y, marker="o", linestyle="", color=color)
            # add text next to each point
            for i, txt in enumerate(group["level_source_specie"]):
                ax2.annotate(txt, (x.iloc[i], y.iloc[i]))

            if turn_on_lr:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                print(f"r-squared: {r_value**2}, slope: {slope}, intercept: {intercept}")
                ax2.plot(x, intercept + slope * x, color=color, alpha=0.2)
        # set a title for the figure
        # set labels for the big plot
        ax1.set_ylim(0, 12)
        ax1.set_xlabel(x_name)
        ax1.set_ylabel("ZFS D (GHz)")
        ax1.set_title(f"ZFS D vs. {x_name}")
        # turn on legend
        ax1.legend(loc="upper left")
        # set space between ax1 and ax2 larger
        plt.subplots_adjust(wspace=9)
        return fig

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
        # for i in [0, 1]:
        #     count = bin1.count() if i == 0 else bin2.count()
        #     axins.text(i, count-1.8, str(count), verticalalignment="center", horizontalalignment="center")

        df1 = df.loc[:, "TDM_rate"]
        # plot bar chart for the distribution of the data in the main plot, with the range from 0 to 175 with increment of 5
        ax.bar(
            # [i+2.5 for i in range(0, 175, 5)],
            [i + 2.5 for i in range(0, 175, 5)],
            [df1.loc[(df1 >= i) & (df1 < i + 5)].count() for i in range(0, 175, 5)],
            width=5,
            edgecolor="black",
        )
        ax.set_xlim(0, 175)

        # set text of count for each bin at the top of the bar and if the count is 0, then remove the text
        # for i in range(len(ax.patches)):
        #     if ax.patches[i].get_height() == 0:
        #         ax.patches[i].set_visible(False)
        #     else:
        #         ax.text(
        #             ax.patches[i].get_x() + ax.patches[i].get_width() / 2 + 2 if i < 2 else ax.patches[i].get_x() +
        #                                                                                     ax.patches[i].get_width() / 2,
        #             ax.patches[i].get_height(),
        #             str(int(ax.patches[i].get_height())),
        #             ha="center",
        #             va="bottom",
        #         )

        # ax.set_xlim(0, 175)
        # ax.set_ylim(0, 30)
        ax.set_xlabel("Dipole transition rate (MHz)")
        ax.set_ylabel("Frequency")
        ax.set_title("Dipole transition rate distribution")
        ax.tick_params(axis="x", which="minor", bottom=True, top=False, labelbottom=False, direction="out")
        ax.tick_params(axis="x", which="major", bottom=True, top=False, direction="out")
        ax.tick_params(axis="y", which="minor", left=False, right=False, labelbottom=True, direction="out")
        ax.tick_params(axis="y", which="major", left=True, right=False, labelbottom=True, direction="out")

        self.figures.append(fig)
        return fig

    def get_tdm_fig_new(self, df=None):
        df = self.get_hse_candidate_df().copy() if df is None else df.copy()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        df1 = df.loc[:, "TDM_rate"]

        # plot histogram of df1
        ax.hist(df1, bins=35, edgecolor="black")

        # add a vertical line at x=6
        ax.axvline(x=6, color="red", linestyle="--")

        ax.set_xlabel("Dipole transition rate (MHz)")
        ax.set_ylabel("Frequency")
        ax.set_title("Dipole transition rate distribution")
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
                                                  "diff_zpl_vertical_transition", "AB", "BC", "CD", "DA",
                                                  "transition_from", "up_tran_cdft_occ",
                                                  "dn_tran_cdft_occ", "dn_in_gap_deg", "up_in_gap_deg"]]


    @staticmethod
    def get_diff_btw_two_dfs(df1, df2):
        a = df1.copy().drop_duplicates(subset=["uid", "defect_name", "charge"], keep="last")
        b = df2.copy().drop_duplicates(subset=["uid", "defect_name", "charge"], keep="last")

        # get intersection of two dataframes and drop duplicates
        c = pd.merge(a, b, on=["uid", "defect_name", "charge"], how="inner") \
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

    def get_mag_ipr_fig(self, df, mag):
        # set up_in_gap_ipr_avg, dn_in_gap_ipr_avg > 1e-5

        nomag_df = df.loc[(df.mag.round(0) == 0) & (df.up_in_gap_ipr_avg > 1e-5) & (df.dn_in_gap_ipr_avg > 1e-5)]
        mag_df = df.loc[(df.mag.round(0) == mag) & (df.up_in_gap_ipr_avg > 1e-5) & (df.dn_in_gap_ipr_avg > 1e-5)]

        fig, ax = plt.subplots(2, 3, figsize=(12, 8))
        # plot hist plot of nomag_df[up_in_gap_ipr_avg] and mag_df[up_in_gap_ipr_avg] and set name of x-axis to
        # "up_in_gap_ipr_avg"
        import seaborn as sns
        ax[0, 0].set_title("Avg. IPRs of spin-up defect levels for vacancies and antisites")
        sns.histplot(nomag_df["up_in_gap_ipr_avg"], ax=ax[0, 0], color="blue", label="Spin-0 defects")
        sns.histplot(mag_df["up_in_gap_ipr_avg"], ax=ax[0, 0], color="red", label=f"Spin-{mag/2} defects")

        # plot lines of average of nomag_df[up_in_gap_ipr_avg] and mag_df[up_in_gap_ipr_avg] on the plot
        ax[0, 0].axvline(nomag_df["up_in_gap_ipr_avg"].mean(), color="blue", linestyle="--", label="Mean IPR of "
                                                                                                   "spin-0 defects")
        ax[0, 0].axvline(mag_df["up_in_gap_ipr_avg"].mean(), color="red", linestyle="--", label=f"Mean IPR of spin-"
                                                                                                f"{mag/2} defects")
        ax[0, 0].set_xlabel("up in gap ipr avg")
        ax[0, 0].set_ylabel("count")
        ax[0, 0].legend()

        # print average of nomag_df[up_in_gap_ipr_avg] and mag_df[up_in_gap_ipr_avg]
        print(
            "up_nomag: {:.3e}, up_mag: {:.3e}".format(
                nomag_df.loc[:, "up_in_gap_ipr_avg"].mean(), mag_df.loc[:,
                                                             "up_in_gap_ipr_avg"].mean()
                )
            )

        # plot hat plot of nomag_df[dn_in_gap_ipr_avg] and mag_df[dn_in_gap_ipr_avg] and set name of x-axis to
        # "dn_in_gap_ipr_avg"
        ax[1, 0].set_title("Avg. IPRs of spin-down defect levels for vacancies and antisites")
        sns.histplot(nomag_df["dn_in_gap_ipr_avg"], ax=ax[1, 0], color="blue", label="Spin-0 defects")
        sns.histplot(mag_df["dn_in_gap_ipr_avg"], ax=ax[1, 0], color="red", label=f"Spin-{mag/2} defects")

        # plot lines of average of nomag_df[dn_in_gap_ipr_avg] and mag_df[dn_in_gap_ipr_avg] on the plot
        ax[1, 0].axvline(nomag_df["dn_in_gap_ipr_avg"].mean(), color="blue", linestyle="--", label=f"Mean IPR of "
                                                                                                   f"spin-0 "
                                                                                                   f"defects")
        ax[1, 0].axvline(mag_df["dn_in_gap_ipr_avg"].mean(), color="red", linestyle="--", label=f"Mean IPR of "
                                                                                                f"spin-{mag/2} "
                                                                                                f"defects")
        ax[1, 0].set_xlabel("dn in gap ipr avg")
        ax[1, 0].set_ylabel("count")
        ax[1, 0].legend()

        print(
            "dn_nomag: {:.3e}, dn_mag: {:.3e}".format(
                nomag_df.loc[:, "dn_in_gap_ipr_avg"].mean(), mag_df.loc[:,
                                                             "dn_in_gap_ipr_avg"].mean()
                )
            )


        # plot hat plot of nomag_df[up_in_gap_ipr_avg] and mag_df[up_in_gap_ipr_avg] and set name of x-axis to "up_in_gap_ipr_avg"
        ax[0, 1].set_title("Avg. IPRs of spin-up defect levels for vacancies")
        sns.histplot(nomag_df.loc[nomag_df.defect_type == "vacancy", "up_in_gap_ipr_avg"], ax=ax[0, 1], color="blue",
                     label=f"Spin-0 vacancies")
        sns.histplot(mag_df.loc[mag_df.defect_type == "vacancy", "up_in_gap_ipr_avg"], ax=ax[0, 1], color="red",
                     label=f"Spin-{mag/2:.0f} vacancies")
        # plot lines of average of nomag_df[up_in_gap_ipr_avg] and mag_df[up_in_gap_ipr_avg] on the plot
        ax[0, 1].axvline(nomag_df.loc[nomag_df.defect_type == "vacancy", "up_in_gap_ipr_avg"].mean(), color="blue",
                         linestyle="--", label="Mean IPR of spin-0 vacancies")
        ax[0, 1].axvline(mag_df.loc[mag_df.defect_type == "vacancy", "up_in_gap_ipr_avg"].mean(), color="red",
                         linestyle="--", label=f"Mean IPR of spin-{mag/2:.0f} vacancies")
        ax[0, 1].set_xlabel("up in gap ipr avg")
        ax[0, 1].set_ylabel("count")
        ax[0, 1].legend()

        print(
            "up_nomag: {:.3e}, up_mag: {:.3e}".format(
                nomag_df.loc[nomag_df.defect_type == "vacancy", "up_in_gap_ipr_avg"].mean
                (), mag_df.loc[mag_df.defect_type == "vacancy", "up_in_gap_ipr_avg"].mean()
                )
            )

        # plot hat plot of nomag_df[dn_in_gap_ipr_avg] and mag_df[dn_in_gap_ipr_avg] and set name of x-axis to "dn_in_gap_ipr_avg"
        ax[1, 1].set_title("Avg. IPRs of spin-down defect levels for vacancies")
        sns.histplot(mag_df.loc[mag_df.defect_type == "vacancy", "dn_in_gap_ipr_avg"], ax=ax[1, 1], color="red",
                     label=f"Spin-{mag/2} vacancies")
        sns.histplot(nomag_df.loc[nomag_df.defect_type == "vacancy", "dn_in_gap_ipr_avg"], ax=ax[1, 1], color="blue",
                     label="Spin-0 vacancies")

        # plot lines of average of nomag_df[dn_in_gap_ipr_avg] and mag_df[dn_in_gap_ipr_avg] on the plot
        ax[1, 1].axvline(nomag_df.loc[nomag_df.defect_type == "vacancy", "dn_in_gap_ipr_avg"].mean(), color="blue",
                         linestyle="--", label="Mean IPR of spin-0 vacancies")
        ax[1, 1].axvline(mag_df.loc[mag_df.defect_type == "vacancy", "dn_in_gap_ipr_avg"].mean(), color="red",
                         linestyle="--", label=f"Mean IPR of spin-{mag/2} vacancies")
        ax[1, 1].set_xlabel("dn in gap ipr avg")
        ax[1, 1].set_ylabel("count")
        ax[1, 1].legend()

        print(
            "dn_nomag: {:.3e}, dn_mag: {:.3e}".format(
                nomag_df.loc[nomag_df.defect_type == "vacancy", "dn_in_gap_ipr_avg"].mean
                (), mag_df.loc[mag_df.defect_type == "vacancy", "dn_in_gap_ipr_avg"].mean()
                )
            )

        # plot hat plot of nomag_df[up_in_gap_ipr_avg] and mag_df[up_in_gap_ipr_avg] and set name of x-axis to "up_in_gap_ipr_avg"
        ax[0, 2].set_title("Avg. IPRs of spin-up defect levels for antisites")
        sns.histplot(nomag_df.loc[nomag_df.defect_type == "antisite", "up_in_gap_ipr_avg"], ax=ax[0, 2],
                     color="blue", label="Spin-0 antisites")
        sns.histplot(mag_df.loc[mag_df.defect_type == "antisite", "up_in_gap_ipr_avg"], ax=ax[0, 2], color="red",
                     label=f"Spin-{mag/2} antisites")

        # plot lines of average of nomag_df[up_in_gap_ipr_avg] and mag_df[up_in_gap_ipr_avg] on the plot
        ax[0, 2].axvline(nomag_df.loc[nomag_df.defect_type == "antisite", "up_in_gap_ipr_avg"].mean(), color="blue",
                            linestyle="--", label="Mean IPR of spin-0 antisites")
        ax[0, 2].axvline(mag_df.loc[mag_df.defect_type == "antisite", "up_in_gap_ipr_avg"].mean(), color="red",
                            linestyle="--", label=f"Mean IPR of spin-{mag/2} antisites")
        ax[0, 2].set_xlabel("up in gap ipr avg")
        ax[0, 2].set_ylabel("count")
        ax[0, 2].legend()

        print(
            "up_nomag: {:.3e}, up_mag: {:.3e}".format(
                nomag_df.loc[nomag_df.defect_type == "antisite", "up_in_gap_ipr_avg"].mean
                (), mag_df.loc[mag_df.defect_type == "antisite", "up_in_gap_ipr_avg"].mean()
                )
            )

        # plot hat plot of nomag_df[dn_in_gap_ipr_avg] and mag_df[dn_in_gap_ipr_avg] and set name of x-axis to "dn_in_gap_ipr_avg"
        ax[1, 2].set_title("Avg. IPRs of spin-down defect levels for antisites")
        sns.histplot(nomag_df.loc[nomag_df.defect_type == "antisite", "dn_in_gap_ipr_avg"], ax=ax[1, 2],
                     color="blue", label="Spin-0 antisites")
        sns.histplot(mag_df.loc[mag_df.defect_type == "antisite", "dn_in_gap_ipr_avg"], ax=ax[1, 2], color="red",
                     label=f"Spin-{mag/2} antisites")
        # plot lines of average of nomag_df[dn_in_gap_ipr_avg] and mag_df[dn_in_gap_ipr_avg] on the plot
        ax[1, 2].axvline(nomag_df.loc[nomag_df.defect_type == "antisite", "dn_in_gap_ipr_avg"].mean(), color="blue",
                            linestyle="--", label=f"Mean IPR of spin-0 antisites")
        ax[1, 2].axvline(mag_df.loc[mag_df.defect_type == "antisite", "dn_in_gap_ipr_avg"].mean(), color="red",
                            linestyle="--", label=f"Mean IPR of spin-{mag/2} antisites")

        ax[1, 2].set_xlabel("dn in gap ipr avg")
        ax[1, 2].set_ylabel("count")
        ax[1, 2].legend()


        print(
            "dn_nomag: {:.3e}, dn_mag: {:.3e}".format(
                nomag_df.loc[nomag_df.defect_type == "antisite", "dn_in_gap_ipr_avg"].mean
                (), mag_df.loc[mag_df.defect_type == "antisite", "dn_in_gap_ipr_avg"].mean()
                )
            )

        for i in range(2):
            for j in range(3):
                #mark subfigures with letters (a), (b), (c), (d), (e), (f) in the upper left corner and outside the
                # plot area and set font as Arial, size 12, weight bolder
                import string
                ax[i, j].text(-0.1, 1.1, string.ascii_lowercase[i * 3 + j], transform=ax[i, j].transAxes,
                                size=14, weight='bold', fontname="Arial")

                ax[i, j].tick_params(axis="x", which="minor", bottom=True, top=False, labelbottom=False,
                                     direction="out")
                ax[i, j].tick_params(axis="x", which="major", bottom=True, top=False, labelbottom=True,
                                        direction="out")
                ax[i, j].tick_params(axis="y", which="minor", left=True, right=False, labelleft=False,
                                        direction="out")
                ax[i, j].tick_params(axis="y", which="major", left=True, right=False, labelleft=True,
                                        direction="out")
                # present x ticks in scientific notation with 3 decimal places
                ax[i, j].ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

                # set x label
                ax[i, j].set_xlabel("IPR (arb. units)")
                # set y label
                ax[i, j].set_ylabel("Count")
                ax[i, j].legend()

        self.figures.append(fig)
        return fig


    def get_ST_en_diff_vs_ipr_fig(self, df, turn_on_lr=True, turn_on_inset=False, exclude_VII=False):
        from scipy import stats
        plot_df = df.copy()

        if exclude_VII:
            plot_df = plot_df.loc[plot_df.chem_group != 17]
        # plot ZPL vs abs_pot_mean and color by groups which are grouped by defect_symmetry and chem_group
        x_value, x_name = "ST_in_gap_occ_ipr", "Tot. occ. defect level IPR"
        # x_value, x_name = "zpl_ipr_avg", "Avg. IPR of LUDL and HODL"

        plot_df.rename(columns={x_value: x_name}, inplace=True)
        plot_df.rename(columns={"singlet_triplet_energy_diff": "Energy diff. between spin singlet and triplet states"},
        inplace=True)
        group_df = plot_df.groupby(["chem_group"])
        fig = plt.figure(figsize=(12, 6))
        # use first 3 rows and first 3 columns as a big plot
        ax1 = plt.subplot2grid((3, 26), (0, 0), colspan=10, rowspan=3)
        colors = ["blue", "red", "lime", "purple", "orange", "cyan"]
        group_name = {6: "VIB", 12: "IIB", 13: "IIIA", 14: "IVA", 15: "VA", 16: "VIA", 17: "VIIA"}

        plot_markers = ['(b)', '(c)', '(d)', '(e)', '(f)', '(g)']
        for (name, group), index, color, plot_mark in zip(group_df, range(len(group_df)), colors, plot_markers):
            name = group_name[name]
            x = group[x_name]
            y = group["Energy diff. between spin singlet and triplet states"]

            ax1.plot(group[x_name], group["Energy diff. between spin singlet and triplet states"],
                     marker="o", linestyle="", label=name, color=color)
            # make values of x axis in scientific notation with 3 decimal places
            ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

            if name != "VIIA" and turn_on_inset:
                # plot an inset figure at the left center of the ax1
                ax1_inset = plt.axes([0.2, 0.45, 0.2, 0.2])
                ax1_inset.plot(x, y, marker="o", linestyle="", label=name, color=color)
                ax1_inset.set_xlim(x.min() * 0.1, x.max() * 1.2)
                ax1_inset.set_ylim(y.min() * 0.1, y.max() * 2)



            # plot scatter plot and regression line for each group at those subplots at last two columns and first 3 rows
            subplot_locs = [(0, 10), (1, 10), (2, 10), (0, 18), (1, 18), (2, 18)]
            ax2 = plt.subplot2grid((3, 26), (subplot_locs[index]), colspan=8, rowspan=1)
            # make ax2 x and y limit increase +- 5% of the max and min values of x and y
            if name == "VIB":
                ax2.set_xlim(x.min() * 0.80, x.max() * 1.1)
                ax2.set_ylim(y.min() * 0.80, y.max() * 1.1)
            elif name == "IIB":
                ax2.set_xlim(x.min() * 0.80, x.max() * 1.1)
                ax2.set_ylim(y.min() * 0.85, y.max() * 1.1)
            elif name == "IIIA":
                ax2.set_xlim(x.min() * 0.95, x.max() * 1.1)
                ax2.set_ylim(y.min() * 0.50, y.max() * 1.1)
            elif name == "IVA":
                ax2.set_xlim(x.min() * 0.95, x.max() * 1.05)
                ax2.set_ylim(y.min() * 0.5, y.max() * 1.1)
            elif name == "VIA":
                ax2.set_xlim(x.min() * 0.95, x.max() * 1.05)
                ax2.set_ylim(y.min() * 0.95, y.max() * 1.1)
            elif name == "VIIA":
                ax2.set_xlim(x.min() * 0.85, x.max() * 1.1)
                ax2.set_ylim(y.min() * 0.95, y.max() * 1.05)

            # make text box at middle top to show name
            ax2.text(0.5, 1.03, name, horizontalalignment="center", verticalalignment="center", transform=ax2.transAxes)
            ax2.plot(x, y, marker="o", linestyle="", color=color)
            ax2.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)


            #make line very light color
            # add text next to each point
            # for i, txt in enumerate(group["level_source_specie"]):
            #     ax2.annotate(txt, (x.iloc[i], y.iloc[i]))

            if turn_on_lr:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                ax2.plot(x, intercept + slope * x, color=color, alpha=0.2)

        # set a title for the figure
        # set labels for the big plot
        ax1.set_xlabel(x_name)
        ax1.set_ylabel("Energy diff. between spin singlet and triplet states (eV)")
        ax1.set_title(f"Energy diff. between spin singlet and triplet states vs. {x_name}")
        # turn on legend
        if exclude_VII:
            ax1.legend(loc="upper left")
        else:
            ax1.legend(loc="upper right")
        # set space between ax1 and ax2 larger
        plt.subplots_adjust(wspace=9)
        return fig




    def get_zpl_vs_pot(self, df):
        import re
        test_df = df.copy()
        def add_nn_info(df):
            entries = HSEQubitDefect.collection.find({"task_id": {"$in": df.task_id.tolist()}})
            # get a dataframe that contains task_id and NN
            d = {
                "task_id": [],
                "NN": [],
                "NN_dist": [],
                "NN_part_charge": [],
                "NN_bader_charge": [],
                "Ant_part_charge": [],
                "Ant_bader_charge": [],
                "NN_atom_vol": [],
                "Ant_atom_vol": [],
                "NN_pot": [],
                "NN_dip": [],
                "NN_part_chg_diff": [],
                "the_other_el_vol": [],
                "NN_dist_inverse": [],
                "NN_part_chg_mult_ant_chg": [],
            }
            for e in entries:
                d["task_id"].append(e["task_id"])
                d["NN"].append(e["NN"])
                # loop over NN except the last one
                nn_dist = []
                nn_dist_inverse = []
                nn_part_charge = []
                nn_charge = []
                nn_atom_vol = []
                nn_part_charge_diff = []
                nn_part_chg_mult_ant_chg = []
                the_other_el_vol = []
                st = Structure.from_dict(e["output"]["structure"])
                for i in e["NN"][:-1]:
                    nn_dist.append(st.get_distance(i, e["NN"][-1]))
                    nn_part_charge.append(-1 * e["Bader"]["charge_transfer"][i])
                    nn_charge.append(e["Bader"]["charge"][i])
                    nn_atom_vol.append(e["Bader"]["atomic_volume"][i])
                    nn_dist_inverse = 1 / np.array(nn_dist)
                    nn_part_chg_mult_ant_chg.append(
                        e["Bader"]["charge_transfer"][i] * e["Bader"]["charge_transfer"][e["NN"][-1]]
                    )
                    nn_part_charge_diff.append(
                        e["Bader"]["charge_transfer"][i] - e["Bader"]["charge_transfer"][e["NN"][-1]]
                        )

                ant_el = st.sites[e["NN"][-1]].specie.symbol
                the_other_el_index = [index for index, i in enumerate(st.sites) if i.specie.symbol != ant_el]
                for el_index in the_other_el_index:
                    the_other_el_vol.append(e["Bader"]["atomic_volume"][el_index])

                d["NN_dist"].append(nn_dist)
                d["NN_dist_inverse"].append(nn_dist_inverse)
                d["NN_part_charge"].append(nn_part_charge)
                d["NN_bader_charge"].append(nn_charge)
                d["NN_atom_vol"].append(nn_atom_vol)
                d["NN_part_chg_mult_ant_chg"].append(nn_part_chg_mult_ant_chg)
                d["NN_pot"].append(np.array(nn_part_chg_mult_ant_chg) * np.array(nn_dist_inverse))
                d["NN_dip"].append(np.array(nn_part_charge_diff) * np.array(nn_dist_inverse) ** 3)
                d["Ant_part_charge"].append(-1 * e["Bader"]["charge_transfer"][e["NN"][-1]])
                d["Ant_bader_charge"].append(e["Bader"]["charge"][e["NN"][-1]])
                d["Ant_atom_vol"].append(e["Bader"]["atomic_volume"][e["NN"][-1]])
                d["NN_part_chg_diff"].append(nn_part_charge_diff)
                d["the_other_el_vol"].append(the_other_el_vol)
            nn_dist_df = pd.DataFrame(d)
            nn_dist_df["the_other_el_vol"] = nn_dist_df["the_other_el_vol"].apply(lambda x: np.mean(x))
            nn_dist_df["Ant_the_other_el_vol_diff"] = nn_dist_df["Ant_atom_vol"] - nn_dist_df["the_other_el_vol"]
            return nn_dist_df

        def get_data_df(df):
            data_df = df[[
                "task_id",
                # "tr",
                # "vertical_transition_up",
                # "vertical_transition_down",
                # "atomic_radius",
                "NN_dist",
                "NN_dist_inverse",
                # "chem_group",
                "NN_part_charge",
                "NN_bader_charge",
                "Ant_part_charge",
                "Ant_bader_charge",
                "NN_atom_vol",
                "Ant_atom_vol",
                "NN_pot",
                "NN_dip",
                "NN_part_chg_diff",
                "Ant_the_other_el_vol_diff",
                # "lattice_a",
                # "lattice_b",
                "NN_part_chg_mult_ant_chg"
            ]]

            # separate the NN_dist into columns named nn_dist_1, nn_dist_2, nn_dist_3, nn_dist_4 ... etc
            data_df = data_df.join(pd.DataFrame(data_df["NN_dist"].tolist(), index=data_df.index))
            # rename columns 0, 1, 2, 3 ... etc to nn_dist_1, nn_dist_2, nn_dist_3, nn_dist_4 ... etc
            data_df = data_df.rename(columns=lambda x: "nn_dist_" + str(x) if isinstance(x, int) else x)
            data_df = data_df.drop(columns=["NN_dist"])

            data_df = data_df.join(pd.DataFrame(data_df["NN_dist_inverse"].tolist(), index=data_df.index))
            data_df = data_df.rename(columns=lambda x: "nn_dist_inverse_" + str(x) if isinstance(x, int) else x)
            data_df = data_df.drop(columns=["NN_dist_inverse"])

            # separate the part_charge into columns named part_charge_1, part_charge_2, part_charge_3, part_charge_4 ... etc
            data_df = data_df.join(pd.DataFrame(data_df["NN_part_charge"].tolist(), index=data_df.index))
            # rename columns 0, 1, 2, 3 ... etc to part_charge_1, part_charge_2, part_charge_3, part_charge_4 ... etc
            data_df = data_df.rename(columns=lambda x: "NN_part_charge_" + str(x) if isinstance(x, int) else x)
            data_df = data_df.drop(columns=["NN_part_charge"])
            # separate the bader_charge into columns named bader_charge_1, bader_charge_2, bader_charge_3, bader_charge_4
            # ... etc
            data_df = data_df.join(pd.DataFrame(data_df["NN_bader_charge"].tolist(), index=data_df.index))
            # rename columns 0, 1, 2, 3 ... etc to bader_charge_1, bader_charge_2, bader_charge_3, bader_charge_4 ... etc
            data_df = data_df.rename(columns=lambda x: "NN_bader_charge_" + str(x) if isinstance(x, int) else x)
            data_df = data_df.drop(columns=["NN_bader_charge"])
            # separate the atom_vol into columns named atom_vol_1, atom_vol_2, atom_vol_3, atom_vol_4 ... etc
            data_df = data_df.join(pd.DataFrame(data_df["NN_atom_vol"].tolist(), index=data_df.index))
            # rename columns 0, 1, 2, 3 ... etc to atom_vol_1, atom_vol_2, atom_vol_3, atom_vol_4 ... etc
            data_df = data_df.rename(columns=lambda x: "NN_atom_vol_" + str(x) if isinstance(x, int) else x)
            data_df = data_df.drop(columns=["NN_atom_vol"])
            # separate the pot into columns named pot_1, pot_2, pot_3, pot_4 ... etc
            data_df = data_df.join(pd.DataFrame(data_df["NN_pot"].tolist(), index=data_df.index))
            # rename columns 0, 1, 2, 3 ... etc to pot_1, pot_2, pot_3, pot_4 ... etca
            data_df = data_df.rename(columns=lambda x: "NN_pot_" + str(x) if isinstance(x, int) else x)
            data_df = data_df.drop(columns=["NN_pot"])
            # separate the dip into columns named dip_1, dip_2, dip_3, dip_4 ... etc
            data_df = data_df.join(pd.DataFrame(data_df["NN_dip"].tolist(), index=data_df.index))
            # rename columns 0, 1, 2, 3 ... etc to dip_1, dip_2, dip_3, dip_4 ... etc
            data_df = data_df.rename(columns=lambda x: "NN_dip_" + str(x) if isinstance(x, int) else x)
            data_df = data_df.drop(columns=["NN_dip"])
            # separate the part_chg_diff into columns named part_chg_diff_1, part_chg_diff_2, part_chg_diff_3,
            # part_chg_diff_4 ... etc
            data_df = data_df.join(pd.DataFrame(data_df["NN_part_chg_diff"].tolist(), index=data_df.index))
            # rename columns 0, 1, 2, 3 ... etc to part_chg_diff_1, part_chg_diff_2, part_chg_diff_3, part_chg_diff_4 ...
            # etc
            data_df = data_df.rename(columns=lambda x: "NN_part_chg_diff_" + str(x) if isinstance(x, int) else x)
            data_df = data_df.drop(columns=["NN_part_chg_diff"])

            data_df = data_df.join(pd.DataFrame(data_df["NN_part_chg_mult_ant_chg"].tolist(), index=data_df.index))
            data_df = data_df.rename(columns=lambda x: "NN_part_chg_mult_ant_chg_" + str(x) if isinstance(x, int) else x)
            data_df = data_df.drop(columns=["NN_part_chg_mult_ant_chg"])
            return data_df



        def replace_nan_with_zero(data_df):
            # find all columns containing nn_dist
            # regex nn_dist_n where n is a number
            nn_dist_cols = [col for col in data_df.columns if re.search("nn_dist_\d+", col)]
            max_len_of_nn_dist = len(nn_dist_cols)
            # calculate mean of nn_dist from 0 to 7 not including any nan value in the list of nn_dist from 0 to
            # max_len_of_nn_dist-1
            data_df["nn_dist_mean"] = data_df.apply(lambda row: np.nanmean(row[nn_dist_cols[0:max_len_of_nn_dist]]), axis=1)
            for i in range(1, max_len_of_nn_dist):
                data_df["nn_dist_" + str(i)] = data_df["nn_dist_" + str(i)].fillna(1000)

            # drop the column "nn_dist_mean"
            # data_df = data_df.drop(columns=["nn_dist_mean"])

            nn_dist_inverse_cols = [col for col in data_df.columns if "nn_dist_inverse" in col]
            max_len_of_nn_dist_inverse = len(nn_dist_inverse_cols)
            data_df["nn_dist_inverse_mean"] = data_df.apply(lambda row: np.nanmean(row[nn_dist_inverse_cols[0:max_len_of_nn_dist_inverse]]), axis=1)
            # data_df["nn_dist_inverse_mean"] = data_df.apply(
            #     lambda x: np.nanmean(
            #         x[["nn_dist_inverse_0", "nn_dist_inverse_1", "nn_dist_inverse_2", "nn_dist_inverse_3",
            #            "nn_dist_inverse_4", "nn_dist_inverse_5", "nn_dist_inverse_6", "nn_dist_inverse_7"]]
            #         ), axis=1
            #     )
            # fill nan value in nn_dist_inverse_1 to nn_dist_inverse_7 with data_df["nn_dist_inverse_mean"]
            for i in range(1, max_len_of_nn_dist_inverse):
                data_df["nn_dist_inverse_" + str(i)] = data_df["nn_dist_inverse_" + str(i)].fillna(0)
            # calculate mean of part_charge from 0 to 7 not including any nan value in the list of part_charge from 0
            # to 7

            part_charge_cols = [col for col in data_df.columns if "NN_part_charge" in col]
            max_len_of_part_charge = len(part_charge_cols)
            data_df["part_charge_mean"] = data_df.apply(lambda row: np.nanmean(row[part_charge_cols[0:max_len_of_part_charge]]), axis=1)
            for i in range(1, max_len_of_part_charge):
                data_df["NN_part_charge_" + str(i)] = data_df["NN_part_charge_" + str(i)].fillna(0)
            # drop the column "part_charge_mean"
            # data_df = data_df.drop(columns=["part_charge_mean"])


            # calculate mean of bader_charge from 0 to 7 not including any nan value in the list of bader_charge from
            # 0 to 7
            bader_charge_cols = [col for col in data_df.columns if "NN_bader_charge" in col]
            max_len_of_bader_charge = len(bader_charge_cols)
            data_df["bader_charge_mean"] = data_df.apply(lambda row: np.nanmean(row[bader_charge_cols[0:max_len_of_bader_charge]]), axis=1)
            for i in range(1, max_len_of_bader_charge):
                data_df["NN_bader_charge_" + str(i)] = data_df["NN_bader_charge_" + str(i)].fillna(0)
            # drop the column "bader_charge_mean"
            # data_df = data_df.drop(columns=["bader_charge_mean"]


            # calculate mean of atom_vol from 0 to 7 not including any nan value in the list of atom_vol from 0 to 7
            atom_vol_cols = [col for col in data_df.columns if "NN_atom_vol" in col]
            max_len_of_atom_vol = len(atom_vol_cols)
            data_df["atom_vol_mean"] = data_df.apply(lambda row: np.nanmean(row[atom_vol_cols[0:max_len_of_atom_vol]]), axis=1)
            for i in range(1, max_len_of_atom_vol):
                data_df["NN_atom_vol_" + str(i)] = data_df["NN_atom_vol_" + str(i)].fillna(0)
            # drop the column "atom_vol_mean"
            # data_df = data_df.drop(columns=["atom_vol_mean"])


            # calculate mean of pot from 0 to 7 not including any nan value in the list of pot from 0 to 7
            pot_cols = [col for col in data_df.columns if "NN_pot" in col]
            max_len_of_pot = len(pot_cols)
            data_df["pot_mean"] = data_df.apply(lambda row: np.nanmean(row[pot_cols[0:max_len_of_pot]]), axis=1)
            for i in range(1, max_len_of_pot):
                data_df["NN_pot_" + str(i)] = data_df["NN_pot_" + str(i)].fillna(0)
            # drop the column "pot_mean"
            # data_df = data_df.drop(columns=["pot_mean"])


            # calculate mean of dip from 0 to 7 not including any nan value in the list of dip from 0 to 7
            dip_cols = [col for col in data_df.columns if "NN_dip" in col]
            max_len_of_dip = len(dip_cols)
            data_df["dip_mean"] = data_df.apply(lambda row: np.nanmean(row[dip_cols[0:max_len_of_dip]]), axis=1)
            for i in range(1, max_len_of_dip):
                data_df["NN_dip_" + str(i)] = data_df["NN_dip_" + str(i)].fillna(0)
            # drop the column "dip_mean"
            # data_df = data_df.drop(columns=["dip_mean"])


            part_charge_cols = [col for col in data_df.columns if "NN_part_chg_diff" in col]
            max_len_of_part_charge = len(part_charge_cols)
            data_df["NN_part_chg_diff_mean"] = data_df.apply(lambda row: np.nanmean(row[part_charge_cols[0:max_len_of_part_charge]]), axis=1)
            for i in range(1, max_len_of_part_charge):
                data_df["NN_part_chg_diff_" + str(i)] = data_df["NN_part_chg_diff_" + str(i)].fillna(0)


            part_chg_mult_cols = [col for col in data_df.columns if "NN_part_chg_mult_ant_chg" in col]
            max_len_of_part_chg_mult = len(part_chg_mult_cols)
            data_df["NN_part_chg_mult_ant_chg_mean"] = data_df.apply(lambda row: np.nanmean(row[part_chg_mult_cols[0:max_len_of_part_chg_mult]]), axis=1)
            for i in range(1, max_len_of_part_chg_mult):
                data_df["NN_part_chg_mult_ant_chg_" + str(i)] = data_df["NN_part_chg_mult_ant_chg_" + str(i)].fillna(0)
            return data_df

        def add_addtional_info():
            df = pd.DataFrame()
            for i, row in test_df.iterrows():
                lattice = Structure.from_dict(
                    C2DB_IR_calc_data.collection.find_one(
                        {"c2db_uid": row["uid"], "task_label": "hse line"}
                    )["output"]["structure"]
                    ).lattice
                ant_specie = Element(row["level_source_specie"])
                # add row of (task_id, lattice_a, lattice_b, lattice_c, is_lattice_a_eq_b)
                lumo_homo_deg = row["dn_tran_lumo_homo_deg"] if row["transition_from"] == 'dn' else row[
                                "up_tran_lumo_homo_deg"]
                df = df.append(
                    pd.DataFrame(
                        {   "task_id": row["task_id"],
                            "defect_type": row["defect_type"],
                            "C2DB_uid": row["C2DB_uid"],
                            "defect_name": row["defect_name"],
                            "charge": row["charge"],
                            "lattice_a": round(lattice.a, 3),
                            "lattice_b": round(lattice.b, 3),
                            "lattice_c": round(lattice.c, 3),
                            "is_lattice_a_eq_b": round(lattice.a, 3) == round(lattice.b, 3),
                            "tr": max(row["vertical_transition_up"], row["vertical_transition_down"]),
                            "chem_row": ant_specie.row,
                            "chem_group": ant_specie.group,
                            "atomic_radius": ant_specie.data["Atomic radius calculated"] if
                            ant_specie.group !=6 else ant_specie.average_cationic_radius,
                            "antisite_ion_type": "cation" if ant_specie.common_oxidation_states[0] > 0 else "anion",
                            "is_transition_metal": ant_specie.is_transition_metal,
                            "ant_elec_neg": ant_specie.X,
                            "level_source_specie": row["level_source_specie"],
                            "defect_symmetry": row["defect_symmetry"],
                            "ZPL": row["ZPL"],
                            "prototype": row["prototype"],
                            "spacegroup": row["spacegroup"],
                            "lumo_homo_deg": str(lumo_homo_deg)

                        }, index=[0]
                    )
                )
            # if ant_specie == "Sn", set antisite_ion_type to "cation"
            df.loc[df["level_source_specie"].isin(["Ge", "Sn"]), "antisite_ion_type"] = "cation"
            display(df.head())
            return df

        def merge_all(df):
            df = add_nn_info(df)
            df = get_data_df(df)
            df = replace_nan_with_zero(df)
            df = add_addtional_info().merge(df, on="task_id")

            df["atom_vol_diff"] = df["Ant_atom_vol"] - df["atom_vol_mean"]
            df["abs_NN_part_chg_diff_mean"] = df["NN_part_chg_diff_mean"].abs()
            df["sign_NN_part_chg_diff_mean"] = df["NN_part_chg_diff_mean"].apply(
                lambda x: "positive" if x > 0 else "negative"
                )
            df["sign_ant_part_charge"] = df["Ant_part_charge"].apply(
                lambda x: "positive" if x > 0 else "negative"
                )
            df["sign_part_charge"] = df["part_charge_mean"].apply(
                lambda x: "positive" if x > 0 else "negative"
                )
            df["abs_pot_mean"] = df["pot_mean"].abs()

            # merge data_df with the column ZPL in sheet.hse_candidate_df
            # test_df = pd.merge(df, test_df.loc[:, ["task_id", "ZPL"]], on="task_id", how="left")
            # test_df = test_df.loc[test_df.ZPL.notnull(), :]
            # # test_df = test_df.loc[test_df.antisite_ion_type=="anions", :]

            return df
        return merge_all(df=test_df)

    def get_zpl_vs_pot_fig(self, df, turn_on_lr=True):
        from scipy import stats
        plot_df = df.copy()
        # plot ZPL vs abs_pot_mean and color by groups which are grouped by defect_symmetry and chem_group
        plot_df.rename(columns={"abs_pot_mean": "Abs. Coulombic potential"}, inplace=True)
        group_df = plot_df.groupby(["chem_group"])
        fig = plt.figure(figsize=(12, 6))
        # use first 3 rows and first 3 columns as a big plot
        ax1 = plt.subplot2grid((3, 26), (0, 0), colspan=10, rowspan=3)
        colors = ["blue", "red", "lime", "purple", "orange", "cyan"]
        group_name = {6: "VIB", 12: "IIB", 13: "IIIA", 14: "IVA", 15: "VA", 16: "VIA", 17: "VIIA"}
        plot_markers = ['(b)', '(c)', '(d)', '(e)', '(f)', '(g)']
        for (name, group), index, color, plot_mark in zip(group_df, range(len(group_df)), colors, plot_markers):
            name = group_name[name]
            ax1.plot(group["Abs. Coulombic potential"], group["ZPL"], marker="o", linestyle="", label=name, color=color)
            x = group["Abs. Coulombic potential"]
            y = group["ZPL"]
            # plot scatter plot and regression line for each group at those subplots at last two columns and first 3 rows
            subplot_locs = [(0, 10), (1, 10), (2, 10), (0, 18), (1, 18), (2, 18)]
            ax2 = plt.subplot2grid((3, 26), (subplot_locs[index]), colspan=8, rowspan=1)
            # make text box at middle top to show name
            ax2.text(0.5, 1.03, name, horizontalalignment="center", verticalalignment="center", transform=ax2.transAxes)
            ax2.plot(x, y, marker="o", linestyle="", color=color)
            # add text in plot_markers at the top left of the small plots and make them bold
            # ax2.text(-0.08, 1.1, plot_mark, horizontalalignment="left", verticalalignment="top",
            #          transform=ax2.transAxes, fontweight="bold")

            #make line very light color
            if turn_on_lr:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                ax2.plot(x, intercept + slope * x, color=color, alpha=0.2)
            # add text next to each point
            # for i, txt in enumerate(group["level_source_specie"]):
            #     ax2.annotate(txt, (x.iloc[i], y.iloc[i]))

            if name == "VIB":
                ax2.set_xlim(x.min() * 0.60, x.max() * 1.1)
                ax2.set_ylim(y.min() * 0.80, y.max() * 1.1)
            elif name == "IIB":
                ax2.set_xlim(x.min() * 0.95, x.max() * 1.05)
                ax2.set_ylim(y.min() * 0.8, y.max() * 1.2)
            elif name == "IIIA":
                ax2.set_xlim(x.min() * 0.8, x.max() * 1.2)
                ax2.set_ylim(y.min() * 0.9, y.max() * 1.1)
            elif name == "IVA":
                ax2.set_xlim(-0.01, x.max() * 1.1)
                ax2.set_ylim(y.min() * 0.9, y.max() * 1.1)
            elif name == "VIA":
                ax2.set_xlim(0, x.max() * 1.1)
                ax2.set_ylim(y.min() * 0.95, y.max() * 1.05)
            elif name == "VIIA":
                ax2.set_xlim(x.min() * 0.85, x.max() * 1.1)
                ax2.set_ylim(y.min() * 0.95, y.max() * 1.05)

        # set a title for the figure
        # set labels for the big plot
        ax1.set_xlabel("Abs. Coulombic potential energy (a.u.)")
        ax1.set_ylabel("ZPL (eV)")
        ax1.set_title("ZPL vs. Abs. Coulombic potential energy")
        # turn on legend
        ax1.legend(loc="upper right")

        # set space between ax1 and ax2 larger
        plt.subplots_adjust(wspace=9)


        return fig


    def get_zpl_vs_prototype(self, df):
        plot_df = df.copy()
        plot_df.rename(columns={"abs_pot_mean": "Abs. Coulombic potential"}, inplace=True)
        # plot sub-plots of Abs. Coulombic potential vs ZPL for each prototype. The number of sub-plots depende on
        # the number of unique prototypes in the df.
        import plotly.express as px
        import plotly.graph_objects as go
        from scipy import stats

        df = plot_df.groupby(
            [ "defect_symmetry", "chem_group"],
            as_index=False
            )

        for name, group in df:
            # do linear fitting for each group
            x = group["Abs. Coulombic potential"]
            y = group["ZPL"]
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            fig = px.scatter(
                group, x="Abs. Coulombic potential", y="ZPL", color="chem_row",
                title="ZPL vs Abs. Coulombic potential for {}-lr: slope:{:.3f} R^2:{:.2%}".format(name, slope,
                                                                                                r_value**2),
                hover_data=["defect_name", "C2DB_uid", "task_id"]
                )
            # plot linear fitting line
            fig.add_trace(
                go.Scatter(
                    x=x, y=slope * x + intercept, mode="lines", name="linear fitting"
                    )
                )
            fig.show()

        return df

            # lambda x: px.scatter(
            #     x, x="Abs. Coulombic potential", y="ZPL", color="antisite_ion_type", title="{} {} {} {}".format(
            #         x["anion_specie_chem_group"].iloc[0], x["defect_symmetry"].iloc[0], x["spacegroup"].iloc[0],
            #         x["prototype"].iloc[0]
            #         ), hover_data=["task_id", "defect_name", "chem_row"]
            #     ).show()
            # )








class DefectAnalysis(SheetCollection):
    def __init__(self, df, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hse_candidate_df = df if df is not None else self.get_hse_candidate_df()
        self.ptmcs = self.hse_candidate_df.loc[
            (self.hse_candidate_df.prototype.isin(["GaS", "GaSe"])) &
            (self.hse_candidate_df.chemsys.isin(["Ga-S", "Ga-Se", "Ga-Te", "In-S", "In-Se", "In-Te"]))]

        self.xm = self.ptmcs.loc[self.ptmcs.charge == -1]
        self.xm_GaS = self.xm.loc[self.xm.prototype == "GaS"]
        self.xm_GaSe = self.xm.loc[self.xm.prototype == "GaSe"]

        self.mx = self.ptmcs.loc[self.ptmcs.charge == 1]
        self.mx_GaS = self.mx.loc[self.mx.prototype == "GaS"]
        self.mx_GaSe = self.mx.loc[self.mx.prototype == "GaSe"]


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

    def get_proj_analysis(self, taskids):
        geo_df = self.get_geo_NN()
        self.xm_GaS = self.xm_GaS.merge(geo_df.loc[:, ["task_id", "nn"]], on="task_id", how="left") \
            if "nn" not in self.xm_GaS.columns else self.xm_GaS

        def get_proj_hor_vs_vert(taskids):
            hor_vs_vert = {"taskid": [], "band":[], "hor-vs-vert": [], "hor-proj": [], "vert-proj": [],
                           "BandName": [], "level_source_specie": [], "chemsys": [], }
            for taskid in taskids:
                e_df = self.xm_GaS.loc[self.xm_GaS["task_id"]==taskid]
                a, b, c, d = Tools.extract_defect_levels_v2_hse(int(e_df["task_id"].iloc[0]), localisation=float(
                    e_df["localisation_threshold"].iloc[0]), edge_tol=(0.5, 0.5), fig_on=None)

                a1_band = []
                a1_band.append(e_df["dn_in_gap_band_index"].iloc[0][-1])
                a1_band = tuple(a1_band)
                e_band = e_df["dn_in_gap_band_index"].iloc[0][:2]

                for bands, band_name in zip([a1_band, e_band], ["a1", "e"]):
                    for band in bands:
                        print(band)
                        hor_vs_vert["BandName"].append(band_name)
                        hor_vs_vert["band"].append(band)
                        hor_vs_vert["level_source_specie"].append(e_df["level_source_specie"].iloc[0])
                        hor_vs_vert["chemsys"].append(e_df["chemsys"].iloc[0])
                        d["vertical_species_proj"] = d.apply(lambda x: x[e_df["nn"].iloc[0][0]], axis=1)
                        d["horizontal_species_proj"] = d.apply(lambda x: sum([x[i] for i in e_df["nn"].iloc[0][1:5]]), axis=1)
                        df = d.loc[(d["orbital"].isin(["s", "pz", "px", "py"])) & (d["spin"]=="-1") &
                                   (d.index==band-1)]
                        display(df)

                        hor_vs_vert["taskid"].append(taskid)
                        hor_vs_vert["hor-vs-vert"].append(df["horizontal_species_proj"].sum()/
                                                          df["vertical_species_proj"].sum())
                        hor_vs_vert["hor-proj"].append(df["horizontal_species_proj"].sum())
                        hor_vs_vert["vert-proj"].append(df["vertical_species_proj"].sum())

            df = pd.DataFrame(hor_vs_vert)
            # sum columns of hor-vs-vert, hor-proj, and vert-proj for the rows with the same taskid and BandName

            df = df.groupby(["taskid", "BandName"]).agg({"hor-vs-vert": "sum", "hor-proj": "sum", "vert-proj": "sum",
                                                         "band": "first", "level_source_specie": "first", "chemsys": "first"})
            df = df.reset_index()
            df = df.sort_values(by=["level_source_specie"])
            df = df.round(3)
            #rear

            display(df)

            fig, axes = plt.subplots(2, 1)
            ax = axes[0]
            #plot taskid vs hor-proj with labels of BandName
            for band in df["BandName"].unique():
                df_band = df.loc[df["BandName"]==band]
                # rearrange the order of chemsys to [In-S, In-Se, In-Te]

                df_band = df_band.sort_values(by=["chemsys"])
                ax.plot(df_band["chemsys"], df_band["hor-proj"], label=band, marker="o")

            ax = axes[1]
            #plot taskid vs vert-proj with labels of BandName
            for band in df["BandName"].unique():
                df_band = df.loc[df["BandName"]==band]
                ax.plot(df_band["chemsys"], df_band["vert-proj"], label=band, marker="o")

            # share xticks and xticklabels
            axes[0].get_shared_x_axes().join(axes[0], axes[1])
            axes[0].set_xticklabels(axes[1].get_xticklabels())
            # set no minor xticks
            from matplotlib import ticker
            axes[0].xaxis.set_minor_locator(ticker.NullLocator())
            axes[1].xaxis.set_minor_locator(ticker.NullLocator())

            #set ylabel for both axes
            axes[0].set_ylabel("Arb. unit")
            axes[1].set_ylabel("Arb. unit")

            #set title for both axes
            axes[0].set_title("Projection of anions")
            axes[1].set_title("Projection of cation")

            # set legend for both axes
            axes[0].legend()
            axes[1].legend()


            fig.tight_layout()
            fig.savefig("anion_cation_proj.pdf", format="pdf")
            return df


        def get_proj_low3_vs_antisite(taskids):
            hor_vs_vert = {"taskid": [], "band":[], "hor-vs-vert": [], "hor-proj": [], "vert-proj": [], "tran-en": [],
                           "BandName": [],  "level_source_specie": [], "chemsys": []}
            for taskid in taskids:
                e_df = self.xm_GaS.loc[self.xm_GaS["task_id"]==taskid]
                a, b, c, d = Tools.extract_defect_levels_v2_hse(
                    int(e_df["task_id"].iloc[0]),
                    localisation=float(e_df["localisation_threshold"].iloc[0]),
                    edge_tol=(0.5, 0.5),
                    fig_on=None
                )
                a1_band = []
                a1_band.append(e_df["dn_in_gap_band_index"].iloc[0][-1])
                a1_band = tuple(a1_band)
                e_band = e_df["dn_in_gap_band_index"].iloc[0][:2]

                for bands, band_name in zip([a1_band, e_band], ["a1", "e"]):
                    for band in bands:
                        print(band)
                        hor_vs_vert["BandName"].append(band_name)
                        hor_vs_vert["band"].append(band)
                        hor_vs_vert["taskid"].append(taskid)
                        hor_vs_vert["level_source_specie"].append(e_df["level_source_specie"].iloc[0])
                        hor_vs_vert["chemsys"].append(e_df["chemsys"].iloc[0])

                        d["vertical_species_proj"] = d.apply(lambda x: x[e_df["nn"].iloc[0][-1]], axis=1)
                        d["horizontal_species_proj"] = d.apply(lambda x: sum([x[i] for i in e_df["nn"].iloc[0][1:4]]), axis=1)
                        df = d.loc[(d["orbital"].isin(["s", "pz", "px", "py"])) & (d["spin"]=="-1") &
                                   (d.index==band-1)]
                        display(df)

                        hor_vs_vert["hor-vs-vert"].append(df["horizontal_species_proj"].sum()/df[
                            "vertical_species_proj"].sum())
                        hor_vs_vert["hor-proj"].append(df["horizontal_species_proj"].sum())
                        hor_vs_vert["vert-proj"].append(df["vertical_species_proj"].sum())
                        hor_vs_vert["tran-en"].append(e_df["dn_tran_en"].iloc[0])

            df = pd.DataFrame(hor_vs_vert)
            df = df.groupby(["taskid", "BandName"]).agg({"hor-vs-vert": "sum", "hor-proj": "sum", "vert-proj": "sum",
                                                         "band": "first", "level_source_specie": "first", "chemsys": "first"})
            df = df.reset_index()
            df = df.sort_values(by=["level_source_specie"])
            df = df.round(3)

            fig, axes = plt.subplots(2, 1)
            ax = axes[0]
            #plot taskid vs hor-proj with labels of BandName
            for band in df["BandName"].unique():
                df_band = df.loc[df["BandName"]==band]
                ax.plot(df_band["chemsys"], df_band["hor-proj"], label=band, marker="o")

            ax = axes[1]
            #plot taskid vs vert-proj with labels of BandName
            for band in df["BandName"].unique():
                df_band = df.loc[df["BandName"]==band]
                ax.plot(df_band["chemsys"], df_band["vert-proj"], label=band, marker="o")

            # share xticks and xticklabels
            axes[0].get_shared_x_axes().join(axes[0], axes[1])
            axes[0].set_xticklabels(axes[1].get_xticklabels())
            # set no minor xticks
            from matplotlib import ticker
            axes[0].xaxis.set_minor_locator(ticker.NullLocator())
            axes[1].xaxis.set_minor_locator(ticker.NullLocator())

            #set ylabel for both axes
            axes[0].set_ylabel("Arb. unit")
            axes[1].set_ylabel("Arb. unit")

            #set title for both axes
            axes[0].set_title("Projection of bottom anions")
            axes[1].set_title("Projection of antisite")

            # set legend for both axes
            axes[0].legend()
            axes[1].legend()

            fig.tight_layout()
            fig.savefig("antisite_bottom_anions_proj.pdf", format="pdf")
            return df

        return get_proj_hor_vs_vert(taskids), get_proj_low3_vs_antisite(taskids)


    def get_geo_NN(self, df=None):
        df = df if df is not None else self.ptmcs
        from pymatgen import PeriodicSite
        dist_df = {"task_id": [], "a":[], "b":[], "gamma":[], "NN_dist": [], "NN_index": [], "distinct_NN_dist": [],
                   "host_NN_dist":[], "host_NN_index":[], "distinct_host_NN_dist":[], "gap_hse_line":[],
                   "d1 in host": [], "d2 in host": [], "d1 in defect": [], "d2 in defect": [], "C2DB_uid": [],
                   "defect_name": [], "charge": [], "nn": []}
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
            dist_df["nn"].append(nn)

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

    def get_host_materials_info(self):
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


    def get_singlet_triplet_en_diff(self, df=None):
        df = df if df is not None else self.hse_candidate_df
        triplets = df["task_id"].unique()
        singlets = {"task_id": [], "singlet_taskid": [], "C2DB_uid": [], "defect_name": [], "charge":[],
                    "singlet_energy": [], "triplet_energy": [], "singlet_triplet_energy_diff": [],}
        for t in triplets:
            triplet_e = HSEQubitDefect.collection.find_one({"task_id": int(t)})

            singlet_e = HSEQubitDefect.collection.find_one({"pc_from": triplet_e["pc_from"], "nupdown_set": 0,
                                                            "defect_entry.name": triplet_e["defect_entry"]["name"],
                                                            "charge_state": triplet_e["charge_state"], "task_label":
                                                                "HSE_scf"})
            if singlet_e is None:
                print(f"Processing triplet taskid:{t}, No singlet found!")
            else:
                print(f"Processing triplet taskid:{t}, Found singlet taskid:{singlet_e['task_id']}!")

            singlets["task_id"].append(int(t))

            if singlet_e:
                singlets["singlet_taskid"].append(singlet_e["task_id"])
            else:
                singlets["singlet_taskid"].append("None")

            if singlet_e:
                singlets["singlet_energy"].append(singlet_e["output"]["energy"])
            else:
                singlets["singlet_energy"].append("None")

            singlets["C2DB_uid"].append(df.loc[df["task_id"] == t]["C2DB_uid"].tolist()[0])
            singlets["defect_name"].append(df.loc[df["task_id"] == t]["defect_name"].tolist()[0])
            singlets["charge"].append(df.loc[df["task_id"] == t]["charge"].tolist()[0])
            singlets["triplet_energy"].append(triplet_e["output"]["energy"])
            if singlet_e:
                singlets["singlet_triplet_energy_diff"].append(
                    singlet_e["output"]["energy"] - triplet_e["output"]["energy"])
            else:
                singlets["singlet_triplet_energy_diff"].append("None")
        self.singlet_df = pd.DataFrame(singlets).round(3)

    def get_defect_symmetry(self, df=None):
        def convert_character_table(t):
            d = {}
            import sympy as sy
            for row in range(len(t)):
                if row == 0:
                    d["class"] = [i for i in t[0].split(" ") if i !=""]
                    mult = []
                    for i in t[row].split(" "):
                        if i != "":
                            try:
                                mult.append(int(i[0]))
                            except Exception:
                                mult.append(1)
                    d["mult"] = mult
                else:
                    v = [i for i in t[row].split(" ") if i !=""]
                    ir = v[0]
                    ir2 = v[1]
                    character = []
                    e = np.exp(2/3*np.pi*1j)
                    for scalar in v[2:]:
                        if scalar == "/3":
                            scalar = 1/3
                        elif scalar == "-/3":
                            scalar = -1/3
                        elif scalar == "/2":
                            scalar = 1/2
                        elif scalar == "-/2":
                            scalar = -1/2
                        elif scalar == "i":
                            scalar = 1j
                        elif scalar == "-i":
                            scalar = -1j
                        elif scalar == "e":
                            scalar = e
                        elif scalar == "e*":
                            scalar = np.conj(e)
                        elif scalar == "-e":
                            scalar = -1*e
                        elif scalar == "-e*":
                            scalar = -1*np.conj(e)
                        elif scalar == "ie":
                            scalar = 1j*e
                        elif scalar == "ie*":
                            scalar = 1j*np.conj(e)
                        elif scalar == "-ie":
                            scalar = -1j*e
                        elif scalar == "-ie*":
                            scalar = -1j*np.conj(e)
                        else:
                            scalar = float(sy.Rational(scalar))
                        character.append(scalar)

                    d[ir] = character
                    d[ir2] = character
            return d

        df = df if df is not None else self.hse_candidate_df
        taskids = df["task_id"].unique()
        defect_symmetry_df = {"task_id": [], "defect_symmetry": [], "defect_sym_character_table": [],}
        for taskid in taskids:
            defect_ir_e = HSEQubitIR.collection.find_one({"prev_fw_taskid": int(taskid)})
            defect_symmetry_df["task_id"].append(taskid)
            if defect_ir_e is not None:
                print(f"{taskid} has defect symmetry")
                defect_symmetry = defect_ir_e["irvsp"]["parity_eigenvals"]["single_kpt"]["(0.0, 0.0, 0.0)"]["up"][
                    "point_group"]
                defect_sym_charac_table = defect_ir_e["irvsp"]["parity_eigenvals"]["single_kpt"]["(0.0, 0.0, 0.0)"][
                    "up"]["pg_character_table"]
                defect_sym_charac_table = convert_character_table(defect_sym_charac_table)
                defect_symmetry_df["defect_symmetry"].append(defect_symmetry)
                defect_symmetry_df["defect_sym_character_table"].append(defect_sym_charac_table)
            else:
                print(f"No defect symmetry found for {taskid}")
                defect_symmetry_df["defec_symmetry"].append("None")
                defect_symmetry_df["defect_sym_character_table"].append("None")

        self.defect_symmetry_df = pd.DataFrame(defect_symmetry_df)

    def get_electronic_config(self, pg, df=None):
        df = df if df is not None else self.hse_candidate_df.copy()
        df = df.loc[df["defect_symmetry"] == pg]

        # if transition_from is "up", df["focused_homo_lumo"] = df[up_tran_lumo_homo_ir
        df["focused_ir"] = df.apply(lambda x: x["up_tran_ir"] if x["transition_from"] == "up" else
        x["dn_tran_ir"], axis=1)

        df["focused_occ"] = df.apply(lambda x: x["up_tran_occ"] if x["transition_from"] == "up" else
        x["dn_tran_occ"], axis=1)

        df["focused_deg"] = df.apply(lambda x: x["up_tran_deg"] if x["transition_from"] == "up" else
        x["dn_tran_deg"], axis=1)

        df["focused_lumo_homo_ir"] = df.apply(lambda x: x["up_tran_lumo_homo_ir"] if x["transition_from"] == "up"
        else x["dn_tran_lumo_homo_ir"], axis=1)

        df["focused_lumo_homo_deg"] = df.apply(lambda x: x["up_tran_lumo_homo_deg"] if x["transition_from"] == "up"
        else x["dn_tran_lumo_homo_deg"], axis=1)

        analysis_df = df.loc[:, ["defect_symmetry", "task_id", "prototype", "defect_name", "focused_ir", "focused_occ",
                                 "focused_deg", "focused_lumo_homo_ir", "focused_lumo_homo_deg"]]
        grouping = ["focused_lumo_homo_ir", "focused_lumo_homo_deg", "prototype",
                    "defect_name", "focused_ir", "focused_deg", "focused_occ"]
        analysis_df_gp = analysis_df.groupby(grouping).agg({"task_id": "unique"})

        return analysis_df, analysis_df_gp

    def get_e_config_summary_df(self):
        get_electronic_config = lambda pg: self.get_electronic_config(pg)
        def D3_df():
            d3_df, d3_df_gp = get_electronic_config("D3")
            d3_df.loc[d3_df["task_id"]==124, "e_config_G"] = "unknown"
            d3_df.loc[d3_df["task_id"]==654, "e_config_G"] = "unknown"
            d3_df.loc[d3_df["task_id"]==124, "e_config"] = "unknown"
            d3_df.loc[d3_df["task_id"]==654, "e_config"] = "unknown"

            d3_df.loc[d3_df["task_id"]==107, "e_config_G"] = "G3.2/G2.0"
            d3_df.loc[d3_df["task_id"]==107, "h_config_G"] = "G2.2/G3.2"
            d3_df.loc[d3_df["task_id"]==107, "e_config"] = "E.2/A2.0"
            d3_df.loc[d3_df["task_id"]==107, "h_config"] = "A2.2/E.2"

            d3_df.loc[d3_df["task_id"]==107, "e2_from"] = "e"
            return d3_df

        def D3h_df():
            d3h_df, d3h_df_gp = get_electronic_config("D3h")
            d3h_df.loc[d3h_df["task_id"]==109, "e_config_G"] = "G5.4/G1.2/G6.2"
            d3h_df.loc[d3h_df["task_id"]==109, "h_config_G"] = "G6.2/G1.0/G5.0"
            d3h_df.loc[d3h_df["task_id"]==109, "e_config"] = 'E".4/A1`.2/E`.2'
            d3h_df.loc[d3h_df["task_id"]==109, "h_config"] = 'E`.2/A1`.0/E".0'
            d3h_df.loc[d3h_df["task_id"]==109, "e2_from"] = "h"

            d3h_df.loc[d3h_df["task_id"]==129, "e_config_G"] = "G4.2/G1.2/G6.2"
            d3h_df.loc[d3h_df["task_id"]==129, "h_config_G"] = "G6.2/G1.0/G4.0"
            d3h_df.loc[d3h_df["task_id"]==129, "e_config"] = 'A2".2/A1`.2/E`.2'
            d3h_df.loc[d3h_df["task_id"]==129, "h_config"] = 'E`.2/A1`.0/A2".0'
            d3h_df.loc[d3h_df["task_id"]==129, "e2_from"] = "h"
            return d3h_df

        def D2d_df():
            d2d_df, d2d_df_gp = get_electronic_config("D2d")
            d2d_df.loc[:, "e_config_G"] = "G4.2/G5.2"
            d2d_df.loc[:, "h_config_G"] = "G5.2/G4.0"
            d2d_df.loc[:, "e_config"] = "B2.2/E.2"
            d2d_df.loc[:, "h_config"] = "E.2/B2.0"
            d2d_df.loc[:, "e2_from"] = "h"
            return d2d_df

        def C3v_df():
            c3v_df, c3v_df_gp = get_electronic_config("C3v")
            c3v_df.loc[c3v_df["task_id"].isin([383, 499]), "e_config_G"] = "G1.2/G3.2/G1.0"
            c3v_df.loc[c3v_df["task_id"].isin([383, 499]), "h_config_G"] = "G1.2/G3.2/G1.0"
            c3v_df.loc[c3v_df["task_id"].isin([383, 499]), "e_config"] = "A1.2/E.2/A1.0"
            c3v_df.loc[c3v_df["task_id"].isin([383, 499]), "h_config"] = "A1.2/E.2/A1.0"
            c3v_df.loc[c3v_df["task_id"].isin([383, 499]), "e2_from"] = "h"

            c3v_df.loc[c3v_df["task_id"].isin(
                [1088, 585, 1007, 220, 2688, 2657, 924, 306, 1234, 1423, 2658, 2671, 870, 13, 785, 868]), "e_config_G"]\
                = "G3.2/G1.0"
            c3v_df.loc[c3v_df["task_id"].isin(
                [1088, 585, 1007, 220, 2688, 2657, 924, 306, 1234, 1423, 2658, 2671, 870, 13, 785, 868]), "h_config_G"] = "G1.2/G3.2"
            c3v_df.loc[c3v_df["task_id"].isin(
                [1088, 585, 1007, 220, 2688, 2657, 924, 306, 1234, 1423, 2658, 2671, 870, 13, 785, 868]), "e_config"] = "E.2/A1.0"
            c3v_df.loc[c3v_df["task_id"].isin(
                [1088, 585, 1007, 220, 2688, 2657, 924, 306, 1234, 1423, 2658, 2671, 870, 13, 785, 868]), "h_config"] = "A1.2/E.2"
            c3v_df.loc[c3v_df["task_id"].isin(
                [1088, 585, 1007, 220, 2688, 2657, 924, 306, 1234, 1423, 2658, 2671, 870, 13, 785, 868]), "e2_from"] = "e"

            c3v_df.loc[c3v_df["task_id"].isin([441, 722, 774, 727, 1033, 929,
                                               2590, 151, 135, 162, 182, 208, 575, 171,
                                               970, 238, 253, 196, 258, 291, 2563,
                                               70
                                               ]), "e_config_G"] = "G1.2/G3.2"

            c3v_df.loc[c3v_df["task_id"].isin([441, 722, 774, 727, 1033, 929,
                                               2590, 151, 135, 162, 182, 208, 575, 171,
                                               970, 238, 253, 196, 258, 291, 2563,
                                               70
                                               ]), "h_config_G"] = "G3.2/G1.0"

            c3v_df.loc[c3v_df["task_id"].isin([441, 722, 774, 727, 1033, 929,
                                               2590, 151, 135, 162, 182, 208, 575, 171,
                                               970, 238, 253, 196, 258, 291, 2563,
                                               70
                                               ]), "e_config"] = "A1.2/E.2"

            c3v_df.loc[c3v_df["task_id"].isin([441, 722, 774, 727, 1033, 929,
                                               2590, 151, 135, 162, 182, 208, 575, 171,
                                               970, 238, 253, 196, 258, 291, 2563,
                                               70
                                               ]), "h_config"] = "E.2/A1.0"


            c3v_df.loc[c3v_df["task_id"].isin([441, 722, 774, 727, 1033, 929,
                                               2590, 151, 135, 162, 182, 208, 575, 171,
                                               970, 238, 253, 196, 258, 291, 2563,
                                               70
                                               ]), "e2_from"] = "h"
            return c3v_df

        self.e_config_summary_df = pd.concat([D3_df(), D3h_df(), D2d_df(), C3v_df()])
        self.e_config_summary_df = self.e_config_summary_df.fillna("unknown")
        self.e_config_summary_df_gp = self.e_config_summary_df.groupby(
            ["defect_symmetry", "focused_lumo_homo_ir", "focused_lumo_homo_deg", "prototype",
             "defect_name", "focused_ir", "focused_deg", "focused_occ", "e_config", "h_config", "e2_from"]).agg(
            {'task_id': "unique"})


        df = self.e_config_summary_df.copy()
        df["LUMO_HOMO_deg"] = df["focused_lumo_homo_deg"]
        df.set_index(["task_id", "LUMO_HOMO_deg"], inplace=True)
        df = df.loc[:, ["e_config", "h_config", "e2_from"]]
        # merge with self.hse_candidate_df with df on task_id
        updated_hse_candidate_df = pd.merge(self.hse_candidate_df, df, how="left", on=["task_id", "LUMO_HOMO_deg"])

        return updated_hse_candidate_df


    def get_antisite_projection(self):
        df = self.hse_candidate_df.copy()
        def add_lumo_homo_projections(x):
            task_id = x["task_id"]
            transition_from = x["transition_from"]
            transition_lumo_homo = x["up_tran_lumo_homo_band_id"] if transition_from == "up" else x[
                "dn_tran_lumo_homo_band_id"]
            lumo_band_id, homo_band_id  = transition_lumo_homo[0], transition_lumo_homo[1]
            a, b, c, d, e = Tools.extract_defect_levels_v2_hse(task_id, localisation=0.1, edge_tol=[1, 1],
                                                               fig_on=None)
            lumo_antisite_proj = a.loc[(a.band_id==lumo_band_id) & (a.spin=="1" if transition_from == "up" else
                                                                    a.spin=="-1"),
                                       "antisite_proj"].iloc[0].round(0)
            homo_antisite_proj = a.loc[(a.band_id==homo_band_id) & (a.spin=="1" if transition_from == "up" else
                                                                    a.spin=="-1"),
                                       "antisite_proj"].iloc[0].round(0)

            return (lumo_antisite_proj, homo_antisite_proj)

        def add_lumo_homo_max_antisite_proj(x):
            lumo_antisite_proj = x["lumo_homo_antisite_proj"][0]
            homo_antisite_proj = x["lumo_homo_antisite_proj"][1]
            lumo_homo_deg = x["LUMO_HOMO_deg"]
            lumo_homo_deg_proj = dict(zip(lumo_homo_deg, [lumo_antisite_proj, homo_antisite_proj]))
            lumo_homo_name_proj = dict(zip(["lumo", "homo"], [lumo_antisite_proj, homo_antisite_proj]))
            # find the key with largest value for lumo_homo_deg_proj
            level_deg_with_max_antisite_proj = max(lumo_homo_deg_proj, key=lumo_homo_deg_proj.get)
            level_name_with_max_antisite_proj = max(lumo_homo_name_proj, key=lumo_homo_name_proj.get)
            return (level_deg_with_max_antisite_proj, level_name_with_max_antisite_proj)

        def add_corr_antisite_proj_and_level_cat(x):
            level_cat = x["level_cat"]
            lumo_antisite_proj, homo_antisite_proj =  x["lumo_homo_antisite_proj"][0], x["lumo_homo_antisite_proj"][1]
            if level_cat == 2:
                if lumo_antisite_proj >= homo_antisite_proj:
                    farest_level_antisite_proj_dominates = True
                else:
                    farest_level_antisite_proj_dominates = False
            elif level_cat == 3:
                if homo_antisite_proj >= lumo_antisite_proj:
                    farest_level_antisite_proj_dominates = True
                else:
                    farest_level_antisite_proj_dominates = False
            else:
                farest_level_antisite_proj_dominates = "unknown"
            return farest_level_antisite_proj_dominates

        df["lumo_homo_antisite_proj"] = df.apply(lambda x: add_lumo_homo_projections(x), axis=1)
        df["lumo_antisite_proj"] = df.apply(lambda x: x["lumo_homo_antisite_proj"][0], axis=1)
        df["homo_antisite_proj"] = df.apply(lambda x: x["lumo_homo_antisite_proj"][1], axis=1)
        df["lumo_homo_max_antisite_proj"] = df.apply(lambda x: add_lumo_homo_max_antisite_proj(x), axis=1)
        df["farest_level_antisite_proj_dominates"] = df.apply(lambda x: add_corr_antisite_proj_and_level_cat(x), axis=1)

        def get_table9(df=df):
            df = df.loc[:, ["prototype", "spacegroup", "uid", "task_id", "defect_name", "charge",
                                          "level_cat","chemsys", "site_oxi_state", "defect_symmetry",
                                          "transition_from", "LUMO-HOMO_down", "LUMO-HOMO_up",
                                          "lumo_homo_antisite_proj", "lumo_homo_max_antisite_proj",
                                          "farest_level_antisite_proj_dominates", "lumo_antisite_proj", "homo_antisite_proj",
                                          "LUMO_HOMO_deg","LUMO_HOMO_splitting_type",
                                          "level_source_specie", "zfs_D", "zfs_E"]].sort_values(["prototype"])

            df["level_source_gs"] = df.apply(lambda x: Element(x["level_source_specie"]).electronic_structure, axis=1)
            df["level_source_group"] = df.apply(lambda x: Element(x["level_source_specie"]).group, axis=1)
            df["level_source_row"] = df.apply(lambda x: Element(x["level_source_specie"]).row, axis=1)

            df.loc[df["task_id"].isin([171, 2563, 196, 182, 349, 135, 208, 258]), "level_cat"] = 2
            df.loc[df["task_id"].isin([1033]), "level_cat"] = 3

            def fun(x):
                level_source_specie = x["level_source_specie"]
                el = Element(level_source_specie)
                gs_term = None
                try:
                    gs_term = el.ground_state_term_symbol
                except Exception:
                    gs_term = "unknown"
                return gs_term

            df["level_source_gs_term"] = df.apply(lambda x: fun(x), axis=1)

            def fun(x):
                transition_from = x["transition_from"]
                lumo_homo_dn = x["LUMO-HOMO_down"]
                lumo_homo_up = x["LUMO-HOMO_up"]
                lumo_homo_gap = lumo_homo_up if transition_from == "up" else lumo_homo_dn
                return lumo_homo_gap

            df["lumo_homo_gap"] = df.apply(lambda x: fun(x), axis=1)
            df["lumo_homo_deg"] = df["LUMO_HOMO_deg"]

            df.loc[df.task_id.isin([2688, 2657]), "level_source_gs_term"] = "5D0.0"

            df = df.loc[~(df["task_id"].isin([383, 499]) & (df["LUMO_HOMO_deg"] == (2,1)))]
            return df

        return get_table9()

    def get_IPR_SCAN(self, df):
        ipr_df = {"task_id": [], "dn_in_gap_ipr": [], "up_in_gap_ipr": []}
        for taskid in df.task_id:
            print("=="*10, taskid, "=="*10)
            ipr_df["task_id"].append(taskid)
            entry = SCAN2dDefect.collection.find_one({"task_id": taskid})
            up_ipr = entry["IPR"]["up"]["ipr"]
            dn_ipr = entry["IPR"]["down"]["ipr"]
            level_dn = df.loc[df.task_id==taskid, "dn_in_gap_band_index"].iloc[0]
            level_up = df.loc[df.task_id==taskid, "up_in_gap_band_index"].iloc[0]
            if level_dn != ():
                print("dn_in_gap_band_index", level_dn)
                level_dn_ipr = []
                for level in level_dn:
                    level_dn_ipr.append(dn_ipr[level-1])
                    # add level_dn_ipr to "dn_ipr" of df
                ipr_df["dn_in_gap_ipr"].append(tuple(level_dn_ipr))
            else:
                ipr_df["dn_in_gap_ipr"].append(())
            if level_up != ():
                print("up_in_gap_band_index", level_up)
                level_up_ipr = []
                for level in level_up:
                    level_up_ipr.append(up_ipr[level-1])
                ipr_df["up_in_gap_ipr"].append(tuple(level_up_ipr))
            else:
                ipr_df["up_in_gap_ipr"].append(())
        self.ipr_df = pd.DataFrame(ipr_df)

    def get_IPR_HSE(self, df, include_singlet_in_gap_occ_ipr=True):
        def get_singlet_in_gap_occ_ipr(singlet_taskid):
            _, _, in_gap_levels, _ = GenerateDefectTable.extract_defect_levels_v2_hse(
                singlet_taskid,
                localisation=0.2,
                selected_bands=None,
                edge_tol=(.5, .5)
            )

            entry = HSEQubitDefect.collection.find_one({"task_id": singlet_taskid})
            up_ipr = entry["IPR"]["up"]["ipr"]
            dn_ipr = entry["IPR"]["down"]["ipr"]


            in_gap_dn, in_gap_dn_occ = in_gap_levels["dn_in_gap_band_index"], in_gap_levels["dn_in_gap_occ"]
            dn_in_gap_df = pd.DataFrame({"in_gap_dn": in_gap_dn, "in_gap_dn_occ": in_gap_dn_occ})
            in_gap_up, in_gap_up_occ = in_gap_levels["up_in_gap_band_index"], in_gap_levels["up_in_gap_occ"]
            up_in_gap_df = pd.DataFrame({"in_gap_up": in_gap_up, "in_gap_up_occ": in_gap_up_occ})

            # keep those rows with occ = 1
            dn_in_gap_df = dn_in_gap_df.loc[dn_in_gap_df.in_gap_dn_occ > 0]
            up_in_gap_df = up_in_gap_df.loc[up_in_gap_df.in_gap_up_occ > 0]

            # get the ipr of dn_in_gap_df and up_in_gap_df
            dn_in_gap_df["in_gap_dn_ipr"] = dn_in_gap_df.in_gap_dn.apply(lambda x: dn_ipr[x-1])
            up_in_gap_df["in_gap_up_ipr"] = up_in_gap_df.in_gap_up.apply(lambda x: up_ipr[x-1])

            # get the sum of ipr of dn_in_gap_df and up_in_gap_df and add to ipr_df
            singlet_in_gap_occ_ipr = dn_in_gap_df.in_gap_dn_ipr.sum() + up_in_gap_df.in_gap_up_ipr.sum()
            return singlet_in_gap_occ_ipr


        ipr_df = {"task_id": [], "dn_lumo_homo_ipr": [], "up_lumo_homo_ipr": [], "transition_from": [],
                  "chem_group": [], "chem_row": [],"up_homo_ipr":[], "triplet_in_gap_occ_ipr": [],
                  "singlet_in_gap_occ_ipr": [], "ST_in_gap_occ_ipr": []}
        # for taskid in df.task_id:
        #loop over rows of ipr_df
        for row in df.itertuples():
            print("=="*10, row.task_id, "=="*10)
            taskid = row.task_id
            ipr_df["task_id"].append(taskid)

            entry = HSEQubitDefect.collection.find_one({"task_id": taskid})

            up_ipr = entry["IPR"]["up"]["ipr"]
            dn_ipr = entry["IPR"]["down"]["ipr"]

            in_gap_dn, in_gap_dn_occ = row.dn_in_gap_band_index, row.dn_in_gap_occ
            dn_in_gap_df = pd.DataFrame({"in_gap_dn": in_gap_dn, "in_gap_dn_occ": in_gap_dn_occ})
            in_gap_up, in_gap_up_occ = row.up_in_gap_band_index, row.up_in_gap_occ
            up_in_gap_df = pd.DataFrame({"in_gap_up": in_gap_up, "in_gap_up_occ": in_gap_up_occ})

            # keep those rows with occ = 1
            dn_in_gap_df = dn_in_gap_df.loc[dn_in_gap_df.in_gap_dn_occ > 0]
            up_in_gap_df = up_in_gap_df.loc[up_in_gap_df.in_gap_up_occ > 0]

            # get the ipr of dn_in_gap_df and up_in_gap_df
            dn_in_gap_df["in_gap_dn_ipr"] = dn_in_gap_df.in_gap_dn.apply(lambda x: dn_ipr[x-1])
            up_in_gap_df["in_gap_up_ipr"] = up_in_gap_df.in_gap_up.apply(lambda x: up_ipr[x-1])

            # get the sum of ipr of dn_in_gap_df and up_in_gap_df and add to ipr_df
            triplet_in_gap_occ_ipr = dn_in_gap_df.in_gap_dn_ipr.sum() + up_in_gap_df.in_gap_up_ipr.sum()
            ipr_df["triplet_in_gap_occ_ipr"].append(triplet_in_gap_occ_ipr)

            if include_singlet_in_gap_occ_ipr:
                singlet_taskid = row.singlet_taskid
                singlet_in_gap_occ_ipr = get_singlet_in_gap_occ_ipr(singlet_taskid)
                singlet_triplet_in_gap_occ_sum_ipr = triplet_in_gap_occ_ipr + singlet_in_gap_occ_ipr
                ipr_df["singlet_in_gap_occ_ipr"].append(singlet_in_gap_occ_ipr)
                ipr_df["ST_in_gap_occ_ipr"].append(singlet_triplet_in_gap_occ_sum_ipr)

            # get "up_in_gap_band_index" and "up_in_gap_occ" of row
            up_in_gap_band_index = getattr(row, "up_in_gap_band_index")
            up_in_gap_occ = getattr(row, "up_in_gap_occ")
            # find the max element of "up_in_gap_band_index" that has "up_in_gap_occ" > 0
            up_homo_band_index = pd.DataFrame(
                list(zip(up_in_gap_band_index, up_in_gap_occ)), columns=["band_index", "occ"]
                ).query("occ > 0").band_index.max()

            ipr_df["up_homo_ipr"].append(up_ipr[up_homo_band_index-1])


            transition_from = row.transition_from
            ipr_df["transition_from"].append(transition_from)

            ant_specie = Element(row.level_source_specie)
            chem_group = ant_specie.group
            chem_row = ant_specie.row
            ipr_df["chem_group"].append(chem_group)
            ipr_df["chem_row"].append(chem_row)

            level_dn = row.dn_tran_lumo_homo_band_index
            level_up = row.up_tran_lumo_homo_band_index
            print(transition_from)
            if level_dn != ():
                print("dn_tran_lumo_homo_band_index", level_dn)
                level_dn_ipr = []
                for level in level_dn:
                    level_dn_ipr.append(dn_ipr[level-1])
                    # add level_dn_ipr to "dn_ipr" of df
                ipr_df["dn_lumo_homo_ipr"].append(tuple(level_dn_ipr))
            else:
                ipr_df["dn_lumo_homo_ipr"].append(())
            if level_up != ():
                print("up_tran_lumo_homo_band_index", level_up)
                level_up_ipr = []
                for level in level_up:
                    level_up_ipr.append(up_ipr[level-1])
                ipr_df["up_lumo_homo_ipr"].append(tuple(level_up_ipr))
            else:
                ipr_df["up_lumo_homo_ipr"].append(())

        self.ipr_df = pd.DataFrame(ipr_df)


class HostAnaylsis:
    def __init__(self, host_db):
        self.host_db = host_db

    def get_element_e_a1_resolved_dos(self, dos):
        from pymatgen import Orbital
        import functools
        from pymatgen.electronic_structure.dos import Dos, add_densities

        z_orbs = {"P": (Orbital.s, Orbital.pz), "D": (Orbital.s, Orbital.dz2)}
        xy_orbs = {"P": (Orbital.px, Orbital.py), "D": (Orbital.dxy, Orbital.dx2)}

        e_dos = {}
        a1_dos = {}

        for site, atom_dos in dos.pdos.items():
            el = site.specie
            ground_state_term = el.ground_state_term_symbol
            el = el.name
            if "S" in ground_state_term or "P" in ground_state_term:
                a1_orbs = z_orbs["P"]
                e_orbs = xy_orbs["P"]
            elif "D" in ground_state_term:
                a1_orbs = z_orbs["D"]
                e_orbs = xy_orbs["D"]
            for orb, pdos in atom_dos.items():
                if orb in e_orbs:
                    if el not in e_dos:
                        e_dos[el] = pdos
                    else:
                        e_dos[el] = add_densities(e_dos[el], pdos)
                if orb in a1_orbs:
                    if el not in a1_dos:
                        a1_dos[el] = pdos
                    else:
                        a1_dos[el] = add_densities(a1_dos[el], pdos)
        print(a1_dos.keys(), e_dos.keys())
        out_dict = {el+":a1": Dos(dos.efermi, dos.energies, densities) for el, densities in a1_dos.items()}
        out_dict.update({el+":e": Dos(dos.efermi, dos.energies, densities) for el, densities in e_dos.items()})
        return out_dict

    def get_element_orb_pdos(self, dos):
        from pymatgen import Orbital
        import functools
        from pymatgen.electronic_structure.dos import Dos, add_densities

        orbs = {"P": (Orbital.s, Orbital.pz, Orbital.px, Orbital.py),
                "D": (Orbital.s, Orbital.dz2, Orbital.dxy, Orbital.dx2, Orbital.dxz, Orbital.dyz)}

        orb_dos = {}
        for site, atom_dos in dos.pdos.items():
            el = site.specie
            ground_state_term = el.ground_state_term_symbol
            el = el.name
            if "S" in ground_state_term or "P" in ground_state_term:
                orbitals = orbs["P"]
            elif "D" in ground_state_term:
                orbitals = orbs["D"]
            for orb, pdos in atom_dos.items():
                if orb in orbitals:
                    if el not in orb_dos:
                        orb_dos[el+f":{orb.name}"] = pdos
                    else:
                        orb_dos[el+f":{orb.name}"]  = add_densities( orb_dos[el+f":{orb.name}"], pdos)
        out_dict = {el: Dos(dos.efermi, dos.energies, densities) for el, densities in orb_dos.items()}
        return out_dict

    def get_host_dos(self, task_id):
        from pymatgen.electronic_structure.plotter import DosPlotter
        from pymatgen import Orbital

        host_db = self.host_db
        dos = self.host_db.get_dos(task_id)
        plotter = DosPlotter(sigma=0.07, stack=True)
        plotter.add_dos_dict(dos.get_element_dos())
        el_fig = plotter.get_plot(xlim=[-5, 8])


        plotter = DosPlotter(sigma=0.07, stack=True)
        element_orb_dos = self.get_element_orb_pdos(dos)
        plotter.add_dos_dict(element_orb_dos)
        el_orb_fig = plotter.get_plot(xlim=[-5, 8])


        plotter = DosPlotter(sigma=0.07, stack=True)
        e_a1_dos = self.get_element_e_a1_resolved_dos(dos)
        plotter.add_dos_dict(e_a1_dos)
        e_a1_fig = plotter.get_plot(xlim=[-5, 8])
        return el_orb_fig, el_fig, e_a1_fig

    def get_host_bs_proj(self, df):
        from pymatgen.electronic_structure.plotter import BSPlotterProjected
        must_group_by = ["C2DB_uid", "defect_name", "level_source_specie", "vbm_max_el", "cbm_max_el", "host_taskid"]
        df = df.loc[:, must_group_by]
        display(df)
        for idx, row in df.iterrows():
            uid = row[0]
            defect_name = row[1]
            level_source_specie = row[2]
            vbm_el = row[3]
            cbm_el = row[4]
            task_id = row[5].split("-")[-1]

            vbm_orbitals = ["s", "pz", "px", "py"] # ["s", "dz2", "dx2", "dxy", "dxz", "dyz"] if Element(vbm_el).group in range(3, 13) else ["s", "pz", "px", "py"]
            cbm_orbitals = ["s", "pz", "px", "py"] #["s", "dz2", "dx2", "dxy", "dxz", "dyz"] if Element(cbm_el).group in range(3, 13) else ["s", "pz", "px", "py"]

            bs = self.host_db.get_band_structure(int(task_id))
            from pymatgen import Spin
            vbm_en, vbm_band = bs.get_vbm()["energy"], bs.get_vbm()["band_index"][Spin.up][0]
            plotter = BSPlotterProjected(bs)
            print(f"vbm: {vbm_en}, en@gamma: {bs.bands[Spin.up][vbm_band][0]}, diff: {bs.bands[Spin.up][vbm_band][0]-vbm_en}")

            print("=="*30, row[-1], uid, defect_name, level_source_specie)

            print("***"*20, f"vbm: {vbm_el}")
            fig_vbm = plotter.get_projected_plots_dots_patom_pmorb(
                {vbm_el: vbm_orbitals},
                dictpa={vbm_el:["all"]},
                sum_atoms={vbm_el: ["all"]},
                sum_morbs={vbm_el: vbm_orbitals[-2:]},
                selected_branches=[1,2,3])


            print("***"*20, f"cbm: {cbm_el}")
            fig_cbm = plotter.get_projected_plots_dots_patom_pmorb(
                {cbm_el: cbm_orbitals},
                dictpa={cbm_el:["all"]},
                sum_atoms={cbm_el: ["all"]},
                sum_morbs={cbm_el: cbm_orbitals[-2:]},
                selected_branches=[1,2,3])

            fig_vbm.show()
            fig_cbm.show()


class DecomposeSymRep:
    def _character_tables(self, pg):
        # revising the character table from phonopy to phonopy.phonon.character_table
        charater_table = {
            '-42m':
                [
                    {'rotation_list': ('E', 'S4', 'C2z', 'C2\'', 'sgd'),
                     'character_table': {'A1': ( 1, 1, 1, 1, 1 ),
                                         'A2': ( 1, 1, 1,-1,-1 ),
                                         'B1': ( 1,-1, 1, 1,-1 ),
                                         'B2': ( 1,-1, 1,-1, 1 ),
                                         'E':  ( 2, 0,-2, 0, 0 )},
                     'mapping_table': {'E'   : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0, 1 )),),
                                       'S4'  : ((( 0, 1, 0 ),
                                                 (-1, 0, 0 ),
                                                 ( 0, 0,-1 )),
                                                (( 0,-1, 0 ),
                                                 ( 1, 0, 0 ),
                                                 ( 0, 0,-1 )),),
                                       'C2z' : (((-1, 0, 0 ),
                                                 ( 0,-1, 0 ),
                                                 ( 0, 0, 1 )),),
                                       'C2\'': (((-1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0,-1 )),
                                                (( 1, 0, 0 ),
                                                 ( 0,-1, 0 ),
                                                 ( 0, 0,-1 )),),
                                       'sgd':  ((( 0,-1, 0 ),
                                                 (-1, 0, 0 ),
                                                 ( 0, 0, 1 )),
                                                (( 0, 1, 0 ),
                                                 ( 1, 0, 0 ),
                                                 ( 0, 0, 1 )),)}},

                ],
            'D2d':
                [
                    {'rotation_list': ('E', '2IC4', 'C2', '2C2`', '2IC2"'),
                     "multiplicity": (1, 2, 1, 2, 2),
                     'character_table': {'A1': ( 1, 1, 1, 1, 1 ),
                                         'A2': ( 1, 1, 1,-1,-1 ),
                                         'B1': ( 1,-1, 1, 1,-1 ),
                                         'B2': ( 1,-1, 1,-1, 1 ),
                                         'E':  ( 2, 0,-2, 0, 0 ),
                                         'G1': ( 1, 1, 1, 1, 1 ),
                                         'G2': ( 1, 1, 1,-1,-1 ),
                                         'G3': ( 1,-1, 1, 1,-1 ),
                                         'G4': ( 1,-1, 1,-1, 1 ),
                                         'G5':  ( 2, 0,-2, 0, 0 )},
                     'mapping_table': {'E': ((( 1, 0, 0 ),
                                              ( 0, 1, 0 ),
                                              ( 0, 0, 1 )),),
                                       '2IC4'  : ((( 0, 1, 0 ),
                                                   (-1, 0, 0 ),
                                                   ( 0, 0,-1 )),
                                                  (( 0,-1, 0 ),
                                                   ( 1, 0, 0 ),
                                                   ( 0, 0,-1 )),),
                                       'C2' : (((-1, 0, 0 ),
                                                ( 0,-1, 0 ),
                                                ( 0, 0, 1 )),),
                                       '2C2`': (((-1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0,-1 )),
                                                (( 1, 0, 0 ),
                                                 ( 0,-1, 0 ),
                                                 ( 0, 0,-1 )),),
                                       '2IC2"':  ((( 0,-1, 0 ),
                                                   (-1, 0, 0 ),
                                                   ( 0, 0, 1 )),
                                                  (( 0, 1, 0 ),
                                                   ( 1, 0, 0 ),
                                                   ( 0, 0, 1 )),)},
                     "function": {
                         'A2': ("Rz"),
                         "B2": ("z"),
                         "E":("(x,y)", "(Rx,Ry)"),
                         "G2": ("Rz"),
                         "G4": ("z"),
                         "G5": ("(x,y)", "(Rx,Ry)")}},

                ],
            '32':
                [
                    {'rotation_list': ('E', 'C3', 'C2\''),
                     'character_table': {'A1': ( 1, 1, 1 ),
                                         'A2': ( 1, 1,-1 ),
                                         'E' : ( 2,-1, 0 )},
                     'mapping_table': {'E'   : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0, 1 )),),
                                       'C3'  : ((( 0,-1, 0 ),
                                                 ( 1,-1, 0 ),
                                                 ( 0, 0, 1 )),
                                                ((-1, 1, 0 ),
                                                 (-1, 0, 0 ),
                                                 ( 0, 0, 1 )),),
                                       'C2\'': ((( 0, 1, 0 ),
                                                 ( 1, 0, 0 ),
                                                 ( 0, 0,-1 )),
                                                (( 1,-1, 0 ),
                                                 ( 0,-1, 0 ),
                                                 ( 0, 0,-1 )),
                                                ((-1, 0, 0 ),
                                                 (-1, 1, 0 ),
                                                 ( 0, 0,-1 )),)}}
                ],
            'D3':
                [
                    {'rotation_list': ('E', '2C3', '3C2'),
                     "multiplicity": (1, 2, 3),
                     'character_table': {'A1': ( 1, 1, 1 ),
                                         'A2': ( 1, 1,-1 ),
                                         'E' : ( 2,-1, 0 ),
                                         'G1': ( 1, 1, 1 ),
                                         'G2': ( 1, 1,-1 ),
                                         'G3' : ( 2,-1, 0 )},
                     'mapping_table': {'E'   : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0, 1 )),),
                                       '2C3'  : ((( 0,-1, 0 ),
                                                  ( 1,-1, 0 ),
                                                  ( 0, 0, 1 )),
                                                 ((-1, 1, 0 ),
                                                  (-1, 0, 0 ),
                                                  ( 0, 0, 1 )),),
                                       '3C2': ((( 0, 1, 0 ),
                                                ( 1, 0, 0 ),
                                                ( 0, 0,-1 )),
                                               (( 1,-1, 0 ),
                                                ( 0,-1, 0 ),
                                                ( 0, 0,-1 )),
                                               ((-1, 0, 0 ),
                                                (-1, 1, 0 ),
                                                ( 0, 0,-1 )),)},
                     "function": {"A2": ("z", "Rz"), "E": ("(x,y)", "(Rx,Ry)"), "G2": ("z", "Rz"), "G3": ("z", "Rz")}},
                ],
            'D3_cartesian':
                [
                    {'rotation_list': ('E', '2C3', '3C2'),
                     "multiplicity": (1, 2, 3),
                     'character_table': {'A1': ( 1, 1, 1 ),
                                         'A2': ( 1, 1,-1 ),
                                         'E' : ( 2,-1, 0 ),
                                         'G1': ( 1, 1, 1 ),
                                         'G2': ( 1, 1,-1 ),
                                         'G3' : ( 2,-1, 0 )},
                     'mapping_table': {'E'   : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0, 1 )),),
                                       '2C3'  : ((( -0.5,-np.sqrt(3)/2, 0 ),
                                                  ( np.sqrt(3)/2,-0.5, 0 ),
                                                  ( 0, 0, 1 )),
                                                 ((-0.5, np.sqrt(3)/2, 0 ),
                                                  (-np.sqrt(3)/2, -0.5, 0 ),
                                                  ( 0, 0, 1 )),),
                                       '3C2': ((( 1, 0, 0 ),
                                                ( 0, -1, 0 ),
                                                ( 0, 0,-1 )),
                                               (( -0.5,np.sqrt(3)/2, 0 ),
                                                ( np.sqrt(3)/2,0.5, 0 ),
                                                ( 0, 0,-1 )),
                                               ((-0.5, -np.sqrt(3)/2, 0 ),
                                                (-np.sqrt(3)/2, 0.5, 0 ),
                                                ( 0, 0,-1 )),)},
                     "function": {"A2": ("z", "Rz"), "E": ("(x,y)", "(Rx,Ry)"), "G2": ("z", "Rz"), "G3": ("z", "Rz")}},
                ],
            '3m':
                [
                    {'rotation_list': ('E', 'C3', 'sgv'),
                     'character_table': {'A1': (1, 1, 1),
                                         'A2': (1, 1,-1),
                                         'E' : (2,-1, 0)},
                     'mapping_table': {'E'   : ((( 1, 0, 0),
                                                 ( 0, 1, 0),
                                                 ( 0, 0, 1)),),
                                       'C3'  : ((( 0,-1, 0),
                                                 ( 1,-1, 0),
                                                 ( 0, 0, 1)),
                                                ((-1, 1, 0),
                                                 (-1, 0, 0),
                                                 ( 0, 0, 1)),),
                                       'sgv' : ((( 0,-1, 0),
                                                 (-1, 0, 0),
                                                 ( 0, 0, 1)),
                                                ((-1, 1, 0),
                                                 ( 0, 1, 0),
                                                 ( 0, 0, 1)),
                                                (( 1, 0, 0),
                                                 ( 1,-1, 0),
                                                 ( 0, 0, 1)),)}}
                ],
            'C3v':
                [
                    {'rotation_list': ('E', '2C3', '3IC2'),
                     "multiplicity": (1, 2, 3),
                     'character_table': {'A1': (1, 1, 1),
                                         'A2': (1, 1,-1),
                                         'E' : (2,-1, 0),
                                         'G1': (1, 1, 1),
                                         'G2': (1, 1,-1),
                                         'G3' : (2,-1, 0)},
                     'mapping_table': {'E'   : ((( 1, 0, 0),
                                                 ( 0, 1, 0),
                                                 ( 0, 0, 1)),),
                                       '2C3'  : ((( 0,-1, 0),
                                                  ( 1,-1, 0),
                                                  ( 0, 0, 1)),
                                                 ((-1, 1, 0),
                                                  (-1, 0, 0),
                                                  ( 0, 0, 1)),),
                                       '3IC2' : ((( 0,1, 0),
                                                  (1, 0, 0),
                                                  ( 0, 0, 1)),
                                                 ((1, -1, 0),
                                                  ( 0, -1, 0),
                                                  ( 0, 0, 1)),
                                                 (( -1, 0, 0),
                                                  ( -1,1, 0),
                                                  ( 0, 0, 1)),)},
                     "function": {
                         "A1": ("z"),
                         "A2": ("Rz"),
                         "E": ("(x,y)", "(Rx,Ry)"),
                         "G1": ("z"),
                         "G2": ("Rz"),
                         "G3": ("(x,y)", "(Rx,Ry)")}
                    }
                ],
            'C3v_diag':
                [
                    {'rotation_list': ('E', '2C3', '3IC2'),
                     "multiplicity": (1, 2, 3),
                     'character_table': {'A1': (1, 1, 1),
                                         'A2': (1, 1,-1),
                                         'E' : (2,-1, 0),
                                         'G1': (1, 1, 1),
                                         'G2': (1, 1,-1),
                                         'G3' : (2,-1, 0)},
                     'mapping_table': {'E'   : ((( 1, 0, 0),
                                                 ( 0, 1, 0),
                                                 ( 0, 0, 1)),),
                                       '2C3'  : ((( np.exp(2*np.pi*1j/3), 0, 0),
                                                  ( 0,  np.exp(-2*np.pi*1j/3), 0),
                                                  ( 0, 0, 1)),
                                                 (( np.exp(-2*np.pi*1j/3), 0, 0),
                                                  (0,  np.exp(2*np.pi*1j/3), 0),
                                                  ( 0, 0, 1)),),
                                       '3IC2' : ((( 0,1, 0),
                                                  (1, 0, 0),
                                                  ( 0, 0, 1)),
                                                 ((0,  np.exp(-2*np.pi*1j/3), 0),
                                                  (  np.exp(2*np.pi*1j/3), 0, 0),
                                                  ( 0, 0, 1)),
                                                 (( 0,  np.exp(2*np.pi*1j/3), 0),
                                                  (  np.exp(-2*np.pi*1j/3),0, 0),
                                                  ( 0, 0, 1)),)},
                     "function": {
                         "A1": ("z"),
                         "A2": ("Rz"),
                         "E": ("(x,y)", "(Rx,Ry)"),
                         "G1": ("z"),
                         "G2": ("Rz"),
                         "G3": ("(x,y)", "(Rx,Ry)")}
                    }
                ],
            'C3v_cartesian':
                [
                    {'rotation_list': ('E', '2C3', '3IC2'),
                     "multiplicity": (1, 2, 3),
                     'character_table': {'A1': (1, 1, 1),
                                         'A2': (1, 1,-1),
                                         'E' : (2,-1, 0),
                                         'G1': (1, 1, 1),
                                         'G2': (1, 1,-1),
                                         'G3' : (2,-1, 0)},
                     'mapping_table': {'E'   : ((( 1, 0, 0),
                                                 ( 0, 1, 0),
                                                 ( 0, 0, 1)),),
                                       '2C3'  : ((( -1/2,-np.sqrt(3)/2, 0),
                                                  ( np.sqrt(3)/2,-1/2, 0),
                                                  ( 0, 0, 1)),
                                                 ((-1/2, np.sqrt(3)/2, 0),
                                                  (-np.sqrt(3)/2, -1/2, 0),
                                                  ( 0, 0, 1)),),
                                       '3IC2' : ((( 1,0, 0),
                                                  (0, -1, 0),
                                                  ( 0, 0, 1)),
                                                 ((-1/2, np.sqrt(3)/2, 0),
                                                  ( np.sqrt(3)/2, 1/2, 0),
                                                  ( 0, 0, 1)),
                                                 (( -1/2, -np.sqrt(3)/2, 0),
                                                  ( -np.sqrt(3)/2,1/2, 0),
                                                  ( 0, 0, 1)),)},
                                       # '3IC2' : ((( -1,0, 0),
                                       #            (0, 1, 0),
                                       #            ( 0, 0, 1)),
                                       #           ((1/2, -np.sqrt(3)/2, 0),
                                       #            ( -np.sqrt(3)/2, -1/2, 0),
                                       #            ( 0, 0, 1)),
                                       #           (( 1/2, np.sqrt(3)/2, 0),
                                       #            ( np.sqrt(3)/2,-1/2, 0),
                                       #            ( 0, 0, 1)),)},
                     "function": {
                         "A1": ("z"),
                         "A2": ("Rz"),
                         "E": ("(x,y)", "(Rx,Ry)"),
                         "G1": ("z"),
                         "G2": ("Rz"),
                         "G3": ("(x,y)", "(Rx,Ry)")}
                    }
                ],
            '-6m2':
                [
                    {'rotation_list': ('E', 'C3', 'C\'2', 'sgh', 'S3', 'sgv'),
                     'character_table': {'A1\''  : ( 1, 1, 1, 1, 1, 1 ),
                                         'A2\''  : ( 1, 1,-1, 1, 1,-1 ),
                                         'E\''   : ( 2,-1, 0, 2,-1, 0 ),
                                         'A1\'\'': ( 1, 1, 1,-1,-1,-1 ),
                                         'A2\'\'': ( 1, 1,-1,-1,-1, 1 ),
                                         'E\'\'' : ( 2,-1, 0,-2, 1, 0 )},
                     'mapping_table': {'E'   : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0, 1 )),),
                                       'C3'  : ((( 0,-1, 0 ),
                                                 ( 1,-1, 0 ),
                                                 ( 0, 0, 1 )),
                                                ((-1, 1, 0 ),
                                                 (-1, 0, 0 ),
                                                 ( 0, 0, 1 )),),
                                       'C\'2': ((( 0,-1, 0 ),
                                                 (-1, 0, 0 ),
                                                 ( 0, 0,-1 )),
                                                (( 1, 0, 0 ),
                                                 ( 1,-1, 0 ),
                                                 ( 0, 0,-1 )),
                                                ((-1, 1, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0,-1 )),),
                                       'sgh' : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0,-1 )),),
                                       'S3'  : ((( 0,-1, 0 ),
                                                 ( 1,-1, 0 ),
                                                 ( 0, 0,-1 )),
                                                ((-1, 1, 0 ),
                                                 (-1, 0, 0 ),
                                                 ( 0, 0,-1 )),),
                                       'sgv' : ((( 0,-1, 0 ),
                                                 (-1, 0, 0 ),
                                                 ( 0, 0, 1 )),
                                                (( 1, 0, 0 ),
                                                 ( 1,-1, 0 ),
                                                 ( 0, 0, 1 )),
                                                ((-1, 1, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0, 1 )),)}}
                ],
            'D3h':
                [
                    {'rotation_list': ('E', 'IC2', '2C3', '2IC6', '3C2`', '3IC2"'),
                     "multiplicity": (1, 1, 2, 2, 3, 3),
                     'character_table': {'A1`'  : ( 1, 1, 1, 1, 1, 1 ),
                                         'A2`'  : ( 1, 1, 1, 1, -1, -1),
                                         'E`'   : (2, 2, -1, -1, 0, 0),
                                         'A1"': ( 1, -1, 1, -1, 1, -1),
                                         'A2"': ( 1, -1, 1, -1, -1, 1),
                                         'E"' : ( 2, -2, -1, 1, 0, 0),
                                         'G1'  : ( 1, 1, 1, 1, 1, 1 ),
                                         'G2'  : ( 1, 1, 1, 1, -1, -1),
                                         'G6'   : (2, 2, -1, -1, 0, 0),
                                         'G3': ( 1, -1, 1, -1, 1, -1),
                                         'G4': ( 1, -1, 1, -1, -1, 1),
                                         'G5' : ( 2, -2, -1, 1, 0, 0)},
                     'mapping_table': {'E'   : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0, 1 )),),
                                       '2C3'  : ((( 0,-1, 0 ),
                                                  ( 1,-1, 0 ),
                                                  ( 0, 0, 1 )),
                                                 ((-1, 1, 0 ),
                                                  (-1, 0, 0 ),
                                                  ( 0, 0, 1 )),),
                                       '3C2`': ((( 0,-1, 0 ),
                                                 (-1, 0, 0 ),
                                                 ( 0, 0,-1 )),
                                                (( 1, 0, 0 ),
                                                 ( 1,-1, 0 ),
                                                 ( 0, 0,-1 )),
                                                ((-1, 1, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0,-1 )),),
                                       'IC2' : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0,-1 )),),
                                       '2IC6'  : ((( 0,-1, 0 ),
                                                   ( 1,-1, 0 ),
                                                   ( 0, 0,-1 )),
                                                  ((-1, 1, 0 ),
                                                   (-1, 0, 0 ),
                                                   ( 0, 0,-1 )),),
                                       '3IC2"' : ((( 0,-1, 0 ),
                                                   (-1, 0, 0 ),
                                                   ( 0, 0, 1 )),
                                                  (( 1, 0, 0 ),
                                                   ( 1,-1, 0 ),
                                                   ( 0, 0, 1 )),
                                                  ((-1, 1, 0 ),
                                                   ( 0, 1, 0 ),
                                                   ( 0, 0, 1 )),)},
                     "function": {"A2`": ("Rz"),
                                  "E`": ("(x,y)"),
                                  'A2"': ("z"),
                                  'E"': ("(Rx,Ry)"),
                                  "G2": ("Rz"),
                                  "G3": ("(x,y)"),
                                  'G4': ("z"),
                                  'G5': ("(Rx,Ry)")}},
                ],
            'D3h_cartesian':
                [
                    {'rotation_list': ('E', 'IC2', '2C3', '2IC6', '3C2`', '3IC2"'),
                     "multiplicity": (1, 1, 2, 2, 3, 3),
                     'character_table': {'A1`'  : ( 1, 1, 1, 1, 1, 1 ),
                                         'A2`'  : ( 1, 1, 1, 1, -1, -1),
                                         'E`'   : (2, 2, -1, -1, 0, 0),
                                         'A1"': ( 1, -1, 1, -1, 1, -1),
                                         'A2"': ( 1, -1, 1, -1, -1, 1),
                                         'E"' : ( 2, -2, -1, 1, 0, 0),
                                         'G1'  : ( 1, 1, 1, 1, 1, 1 ),
                                         'G2'  : ( 1, 1, 1, 1, -1, -1),
                                         'G6'   : (2, 2, -1, -1, 0, 0),
                                         'G3': ( 1, -1, 1, -1, 1, -1),
                                         'G4': ( 1, -1, 1, -1, -1, 1),
                                         'G5' : ( 2, -2, -1, 1, 0, 0)},
                     'mapping_table': {'E'   : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0, 1 )),),
                                       'IC2' : ((( 1, 0, 0 ),
                                                 ( 0, 1, 0 ),
                                                 ( 0, 0,-1 )),),

                                       '2C3'  : ((( -0.5,-np.sqrt(3)/2, 0 ),
                                                  ( np.sqrt(3)/2,-0.5, 0 ),
                                                  ( 0, 0, 1 )),
                                                 ((-0.5, np.sqrt(3)/2, 0 ),
                                                  (-np.sqrt(3)/2, -0.5, 0 ),
                                                  ( 0, 0, 1 )),),
                                       '2IC6'  : ((( -0.5,-np.sqrt(3)/2, 0 ),
                                                   ( np.sqrt(3)/2,-0.5, 0 ),
                                                   ( 0, 0,-1 )),
                                                  ((-0.5, np.sqrt(3)/2, 0 ),
                                                   (-np.sqrt(3)/2, -0.5, 0 ),
                                                   ( 0, 0,-1 )),),

                                       '3C2`': ((( 1,0, 0 ),
                                                 (0, -1, 0 ),
                                                 ( 0, 0,-1 )),
                                                (( -0.5, np.sqrt(3)/2, 0 ),
                                                 ( np.sqrt(3)/2,0.5, 0 ),
                                                 ( 0, 0,-1 )),
                                                ((-0.5, -np.sqrt(3)/2, 0 ),
                                                 ( -np.sqrt(3)/2, 0.5, 0 ),
                                                 ( 0, 0,-1 )),),

                                       '3IC2"' : ((( 1,0, 0 ),
                                                   (0, -1, 0 ),
                                                   ( 0, 0, 1 )),
                                                  (( -0.5, np.sqrt(3)/2, 0 ),
                                                   ( np.sqrt(3)/2,0.5, 0 ),
                                                   ( 0, 0, 1 )),
                                                  ((-0.5, -np.sqrt(3)/2, 0 ),
                                                   ( -np.sqrt(3)/2, 0.5, 0 ),
                                                   ( 0, 0, 1 )),)},
                     "function": {"A2`": ("Rz"),
                                  "E`": ("(x,y)"),
                                  'A2"': ("z"),
                                  'E"': ("(Rx,Ry)"),
                                  "G2": ("Rz"),
                                  "G3": ("(x,y)"),
                                  'G4': ("z"),
                                  'G5': ("(Rx,Ry)")}},
                ],
        }
        return charater_table[pg][0]

    def _get_spin_conjugate_reps(self, n_electrons, pg, ir):
        if n_electrons > 4:
            raise ValueError("n_electrons must be less than 5")
        table_pg = self._character_tables(pg)
        op_to_char = dict(zip(table_pg["rotation_list"], table_pg["character_table"][ir]))
        origin_ir_vec = np.array(list(op_to_char.values()))

        eqivalent_ops_dict =lambda x: self._get_power_n_operation(pg, x)
        n_power_rep_vec = lambda x: np.array([op_to_char[i] for i in eqivalent_ops_dict(x).values()])

        rep_vec = {}
        if n_electrons == 1:
            doublet = n_power_rep_vec(1)
            rep_vec[n_electrons] = {"1/2": doublet} #Sn, spin, rep_vec
            return rep_vec
        elif n_electrons == 2:
            singlet = 1/2*origin_ir_vec**2 + 1/2*n_power_rep_vec(2)
            triplet = 1/2*origin_ir_vec**2 - 1/2*n_power_rep_vec(2)
            rep_vec[n_electrons] = {"0": singlet, "1": triplet}
            return rep_vec
        elif n_electrons == 3:
            doublet = 1/3*origin_ir_vec**3 - 1/3*n_power_rep_vec(3)
            quartet = 1/6*origin_ir_vec**3 - 1/2*origin_ir_vec*n_power_rep_vec(2) + 1/3*n_power_rep_vec(3)
            rep_vec[n_electrons] = {"1/2": doublet, "3/2": quartet}
            return rep_vec
        elif n_electrons == 4:
            singlet = 1/12*origin_ir_vec**4 - 1/3*origin_ir_vec*n_power_rep_vec(3) + 1/4*n_power_rep_vec(2)**2
            triplet = 1/8*origin_ir_vec**4 - 1/4*origin_ir_vec**2*n_power_rep_vec(2) + 1/4*n_power_rep_vec(4) - \
                      1/8*n_power_rep_vec(2)**2
            quintet = 1/24*origin_ir_vec**4 - 1/4*origin_ir_vec**2*n_power_rep_vec(2) + \
                      1/3*origin_ir_vec*n_power_rep_vec(3) - 1/4*n_power_rep_vec(4) + 1/8*n_power_rep_vec(2)**2
            rep_vec[n_electrons] = {"0": singlet, "1": triplet, "2": quintet}
            return rep_vec


    def _get_power_n_operation(self, pg, power_n):
        table_pg = self._character_tables(pg)
        operations = table_pg["mapping_table"]

        eqivalent_ops_dict = {}
        for src_op_name, src_ops in operations.items():
            src_op = src_ops[0] # pick first op for each class
            power_n_operation = np.linalg.matrix_power(src_op, power_n)
            for op_name, ops in operations.items():
                for op in ops:
                    if np.array_equal(power_n_operation, np.array(op)):
                        eqivalent_ops_dict[src_op_name] = op_name
                        break

        eqivalent_ops_dict = {op: eqivalent_ops_dict[op] for op in table_pg["rotation_list"]}
        print(f"Power of {power_n} Ops:\n{eqivalent_ops_dict}")

        return eqivalent_ops_dict



    def _decomposition_engine(self, pg, rep):
        table_pg = self._character_tables(pg)
        rep = np.repeat(rep, table_pg["multiplicity"])
        irs, irs_G = {}, {}

        for k,v in table_pg["character_table"].items():
            if "G" in k:
                irs_G[k] = np.repeat(v, table_pg["multiplicity"])
            else:
                irs[k] = np.repeat(v, table_pg["multiplicity"])

        decomp_coeff, decomp_coeff_G = {}, {}
        for k, v in irs.items():
            decomp_coeff[k] = np.dot(rep, v)/np.sum(table_pg["multiplicity"])

        for k, v in irs_G.items():
            decomp_coeff_G[k] = np.dot(rep, v)/np.sum(table_pg["multiplicity"])
        return decomp_coeff, decomp_coeff_G

    def _get_tensor_operators(self, config='xy*xy', pg='C3v'):
        table_pg = self._character_tables(pg)
        tensor_table = {}
        for op_name, ops in table_pg["mapping_table"].items():
            print(op_name, ops)
            tensor_ops = []
            for idx, op in enumerate(ops):
                op = np.array(op)[:2, :2]
                if config == 'xy*xy':
                    tensor_op = np.kron(op, op)
                    tensor_ops.append(tensor_op)
                elif config == 'xy*z':
                    for ir, function in table_pg["function"].items():
                        if "z" in function:
                            tensor_op = np.kron(op, table_pg["character_table"][ir][idx])
                            tensor_ops.append(tensor_op)
                            break
                    tensor_ops.append(tensor_op)
            tensor_table.update({op_name: tensor_ops})
        return tensor_table

    def decomp_product(self, pg, *irs):
        table_pg = self._character_tables(pg)
        for idx, i in enumerate(irs):
            if idx == 0:
                rep = table_pg["character_table"][i]
                continue
            rep *= np.array(table_pg["character_table"][i])
        print(f"Point_group:: {pg}, Rep: {rep}")
        return self._decomposition_engine(pg, rep)

    def decomp_any_rep(self, pg, rep):
        table_pg = self._character_tables(pg)
        print(f"Point_group: {pg}, Rep: {rep}")
        return self._decomposition_engine(pg, rep)

    def get_term_symbol(self, n_electrons, pg, ir):
        from collections import defaultdict
        term_symbol = {}
        term_symbol["n_electrons"] = n_electrons
        term_symbol["pg"] = pg
        term_symbol["ir"] = ir
        spin_conjugate_rep = self._get_spin_conjugate_reps(n_electrons, pg, ir)
        for spin, rep in spin_conjugate_rep[n_electrons].items():
            print(f"Spin {spin}, Rep: {rep}")
            term_symbol[spin] =  self.decomp_any_rep(pg, rep)
        return term_symbol

    def get_sublevel_sym(self, pg, term_symbol):
        table_pg = self._character_tables(pg)
        singlet = "G1"
        ang_mom_xy = []
        ang_mom_z = []
        for k, v in table_pg["function"].items():
            if "(Rx,Ry)" in v:
                ang_mom_xy.append(k)
            elif "Rz" in v:
                ang_mom_z.append(k)
        print(f"Sx, Sy: {ang_mom_xy}, Sz: {ang_mom_z}")

        ang_mom_xy = [i for i in ang_mom_xy if "G" not in i]
        ang_mom_z = [i for i in ang_mom_z if "G" not in i]

        triplet_sublevel_sym = {"(Sx,Sy)": self.decomp_product(pg, term_symbol, *ang_mom_xy),
                                "Sz": self.decomp_product(pg, term_symbol, *ang_mom_z)}
        singlet_sublevel_sym = {"S0": self.decomp_product(pg, term_symbol, singlet)}

        return triplet_sublevel_sym, singlet_sublevel_sym

    def get_term_symbol_sublevel_sym(self, n_electrons, pg, ir):
        if n_electrons != 2:
            raise ValueError("n_electrons must be 2")

        term_symbol = self.get_term_symbol(n_electrons, pg, ir)
        singlet_terms = term_symbol["0"][0]
        # get those keys whose value is not 0
        singlet_terms = [k for k, v in singlet_terms.items() if v != 0]
        triplet_terms = term_symbol["1"][0]
        triplet_terms = [k for k, v in triplet_terms.items() if v != 0]

        sublevel_sym = {"0": {}, "1": {}}
        for term in singlet_terms:
            print(f"Singlet term: {term}")
            sublevel_sym["0"][term] = self.get_sublevel_sym(pg, term)[1]
        for term in triplet_terms:
            print(f"Triplet term: {term}")
            sublevel_sym["1"][term] = self.get_sublevel_sym(pg, term)[0]
        return sublevel_sym


    # a numpy vector, each element of it times a set of matrices in the basis
    def get_sym_adopted_basis(self, ir, pg="D2d", config="xy*xy"):
        table_pg = self._character_tables(pg)
        tensor_table = self._get_tensor_operators(config=config, pg=pg)
        ir_vec = table_pg["character_table"][ir]
        ir_vec = np.repeat(ir_vec, table_pg["multiplicity"])
        tensor_ops = np.concatenate(list(tensor_table.values()))
        print(f"\nir_vec: {ir_vec}, ir: {ir}")
        print(f"\ntensor_ops:\n{tensor_ops}")
        proj_op = np.array([m*v for m, v in zip(tensor_ops, ir_vec)])
        print(f"\nproj_op:\n{proj_op}")

        new_basis = sum(proj_op)/len(ir_vec)*ir_vec[0]
        print(f"\nnew_basis:\n{new_basis}")

        import scipy as sp
        orthonormal_basis = sp.linalg.orth(new_basis)
        print(f"\northonormal_basis:\n{orthonormal_basis}")
        return new_basis, orthonormal_basis









if __name__ == '__main__':
    a = {1:2, 2:3}
    #find the key who has the largest value
    print(max(a, key=a.get))






