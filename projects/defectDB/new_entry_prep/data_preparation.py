from JPack_independent.projects.defectDB.analysis.data_analysis import *
from JPack_independent.projects.defectDB.analysis.analysis_api import *
from JPack_independent.projects.defectDB.wf.wf import INPUT_PATH

from qubitPack.tool_box import get_db, get_good_ir_sites, get_unique_sites_from_wy, get_band_edges_characters, IOTools

from pymatgen.io.vasp.inputs import Structure
from monty.json import jsanitize
import pandas as pd

import os

C2DB_IR_calc_data = get_db("C2DB_IR", "calc_data", user="Jeng", password="qimin", port=12345)
C2DB_IR_ir_data = get_db("C2DB_IR", "ir_data", user="Jeng_ro", password="qimin", port=12345)
C2DB = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan", port=12345)

"""
Prepare the data for the new entry in the database, and then to update the database.
Generate tables for the new entry.

"""
class DataPrepHost:
    @classmethod
    def sym_data(cls):
        # 2

        # db = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=1236)
        db = C2DB_IR_calc_data

        for e in db.collection.find({"task_label": "hse line",}):
            symprec = 1e-1
            st = Structure.from_dict(e["output"]["structure"])
            species, site_syms, wyckoffs, pg, site_idx, spg, spg_number = get_unique_sites_from_wy(st, symprec=symprec).values()
            good_ir = get_good_ir_sites(st, symprec=symprec)

            print(e["c2db_uid"])
            d = {}
            wy_fields = {
                "species": species,
                "site_sym": site_syms,
                "wyckoffs": wyckoffs,
                "site_idx": site_idx,
            }

            if good_ir["site_syms"]:
                ir = {
                    "species": good_ir["species"],
                    "site_sym": good_ir["site_syms"],
                    "wyckoffs": good_ir["wys"],
                    "site_idx": good_ir["site_idx"],
                }
            else:
                ir = {
                    "species": [],
                    "site_sym": [],
                    "wyckoffs": [],
                    "site_idx": [],
                }

            d["unique_wyckoff"] = wy_fields
            d["good_ir_info"] = ir

            d = jsanitize(d)
            d.update({"symprec":symprec, "pmg_pg": pg, "pmg_spg": spg, "pmg_spg_number": spg_number})

            db.collection.update_one({"task_id": e["task_id"]}, {"$set":{"sym_data":d}})

    @classmethod
    def band_edge(cls):
        # 3
        db = C2DB_IR_calc_data
        col = db.collection
        filter = {"task_label": "hse line", "task_id": 8565}

        for e in list(col.find(filter)):
            print(e["task_id"], e["c2db_uid"])
            bs = db.get_band_structure(e["task_id"])
            if bs.get_band_gap()["energy"] == 0:
                print("skip")
                continue
            data = get_band_edges_characters(bs, ir_data={"db": C2DB_IR_ir_data,
                                                          "entry_filter": {"c2db_uid": e["c2db_uid"]}})
            d = jsanitize(data)
            db.collection.update_one({"task_id": e["task_id"]}, {"$set":{"band_edges":d}})

    @classmethod
    def site_oxi_state(cls):
        from monty.json import jsanitize
        from qubitPack.tool_box import get_db
        from pymatgen import Structure

        db = C2DB_IR_calc_data
        col = db.collection
        filter = {"task_label": "hse line"}

        for e in list(col.find(filter))[:]:
            print(e["task_id"], e["c2db_uid"])
            structure = Structure.from_dict(e["output"]["structure"]).copy()
            structure.add_oxidation_state_by_guess()
            site_oxi_state = []
            for site in structure.sites:
                specie = site.specie.as_dict()
                site_oxi_state.append((specie["element"], specie["oxidation_state"]))
            data = {"site_oxi_state": site_oxi_state}
            d = jsanitize(data)
            db.collection.update_one({"task_id": e["task_id"]}, {"$set": data})

    @classmethod
    def cp_c2db_info(cls):
        src = C2DB
        tgt = C2DB_IR_calc_data

        col = tgt.collection
        filter = {}

        for e in col.find({"task_label": "hse line", "c2db_info": {"$exists": 0}}):
            try:
                print(e["task_id"], e["c2db_uid"])
                src_entry = src.collection.find_one({"uid": e["c2db_uid"]})
                tgt.collection.update_one({"task_id": e["task_id"]}, {"$set":{
                    "c2db_info.uid": src_entry["uid"],
                    "c2db_info.formula": src_entry["formula"],
                    "c2db_info.spacegroup": src_entry["spacegroup"],
                    "c2db_info.gap_hse": src_entry["gap_hse"],
                    "c2db_info.gap_hse_nosoc": src_entry["gap_hse_nosoc"],
                    "c2db_info.gap_nosoc": src_entry["gap_nosoc"],
                    "c2db_info.soc_intensity": src_entry["gap_hse_nosoc"] - src_entry["gap_hse"],
                    "c2db_info.ehull": src_entry["ehull"],
                    "c2db_info.class": src_entry.get("class"),
                    "c2db_info.prototype": src_entry["prototype"],
                }})
            except Exception as er:
                print(er)


class DataPrepDefect:

    @classmethod
    def cp_symdata_bandedges(cls, filter, read_c2db_uid_key=False):
        # 4
        from qubitPack.tool_box import get_db, get_band_edges_characters
        from monty.json import jsanitize

        src = C2DB_IR_calc_data
        tgt = HSEQubitDefect

        col = tgt.collection
        entries = list(col.find(filter))
        for idx, e in enumerate(entries):
            c2db_uid = e["pc_from"].split("/")[-1] if not read_c2db_uid_key else e["c2db_uid"]
            try:
                print(e["task_id"], e["pc_from"])
                src_entry = src.collection.find_one({"c2db_uid":c2db_uid, "task_label": "hse line"})
                # vacuum = self.c2db.collection.find_one({"uid": c2db_uid})["evacmean"]
                sym_data = src_entry["sym_data"]
                hse_bs_output = src_entry["output"]
                for remove in ["structure", "density", "energy", "energy_per_atom", "forces", "stress", "spacegroup"]:
                    hse_bs_output.pop(remove)
                band_edges = src_entry["band_edges"]
                hse_bs_output.update({"band_edges": band_edges})
                tgt.collection.update_one({"task_id": e["task_id"]}, {"$set":{
                    "host_info.sym_data": sym_data,
                    "host_info.c2db_ir_hse_line": hse_bs_output,
                }})
            except Exception as er:
                print(er)

    @classmethod
    def is_site_sym_uniform(cls, filter):
        # 5
        from qubitPack.tool_box import get_db, get_band_edges_characters
        from monty.json import jsanitize

        tgt = HSEQubitDefect

        col = tgt.collection
        entries = list(col.find(filter))
        for e in entries:
            print(e["task_id"], e["pc_from"])
            try:
                site_syms = e["host_info"]["sym_data"]["unique_wyckoff"]["site_sym"]
                if site_syms[0] == site_syms[1]:
                    site_symmetry_uniform = True
                else:
                    site_symmetry_uniform = False

                tgt.collection.update_one({"task_id": e["task_id"]},
                                          {"$set":{"site_symmetry_uniform": site_symmetry_uniform}})
            except Exception as er:
                print(er)

    @classmethod
    def cp_site_oxi_state(cls, filter, read_c2db_uid_key=False):
        # 8
        from monty.json import jsanitize
        from qubitPack.tool_box import get_db
        from pymatgen import Structure

        host = C2DB_IR_calc_data
        db = HSEQubitDefect

        col = db.collection
        entries = list(col.find(filter))
        for idx, e in enumerate(entries):
            c2db_uid = e["pc_from"].split("/")[-1] if not read_c2db_uid_key else e["c2db_uid"]
            print(e["task_id"], e["pc_from"])
            host_entry = host.collection.find_one({"c2db_uid": c2db_uid, "task_label": "hse line"})
            try:
                site_oxi_state = host_entry["site_oxi_state"]
            except Exception:
                site_oxi_state = {}
            db.collection.update_one({"task_id": e["task_id"]}, {"$set": {"host_info.c2db_ir_hse_line.site_oxi_state":
                                                                              site_oxi_state}})

class GenerateDefectTable(BackProcess):
    def __init__(self, df_filter):
        super(GenerateDefectTable, self).__init__(None)
        self.defect_df = None
        self.df_filter = df_filter


    def extract_defect_levels_v2(self, defect_taskid, localisation=0.05, edge_tol=(0.25, 0.25), selected_bands=None):
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
                "eigen",
                None,
                None,  #(host_db, host_taskid, 0, vbm_dx, cbm_dx),
                localisation,  #0.2
                locpot_c2db=None,  #(c2db, c2db_uid, 0)
                is_vacuum_aligment_on_plot=True,
                edge_tol=edge_tol,
                ir_db=ir_db,
                ir_entry_filter={"pc_from_id": pc_from_id, "defect_name": defect_name, "charge_state": charge_state},
                selected_bands=selected_bands,
            )
            tot, proj, d_df, levels, defect_levels = state
            level_info = d_df.to_dict("records")[0]
        except Exception as er:
            print(er)
            level_info = {}
            levels = {}
            defect_levels = {}
        return level_info, levels, defect_levels

    def get_defect_df_v2(self):
        data = []
        col = SCAN2dDefect.collection
        # condition = {"task_label": "SCAN_scf", "task_id": {"$lte": 6084}}
        for e in list(col.find(self.df_filter))[:]:
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

            localisation, selected_bands, edge_tol = 0.05, None, (0.5, 0.5)
            d_df, levels, in_gpa_levels = self.extract_defect_levels_v2(
                e["task_id"], localisation=localisation,
                selected_bands=selected_bands,
                edge_tol=edge_tol
            )
            info.update({"localisation_threshold": localisation})
            info.update(d_df)
            info.update(levels)
            info.update(in_gpa_levels)
            data.append(info)

        self.defect_df = pd.DataFrame(data)
        self.defect_df.fillna("None", inplace=True)
        self.defect_df.replace({"up_tran_en": "None"}, 0, inplace=True)
        self.defect_df.replace({"dn_tran_en": "None"}, 0, inplace=True)

        self.input_df = self.defect_df

        # IOTools(cwd=save_xlsx_path, pandas_df=self.defect_df).to_excel("defect")



    def extract_defect_levels_v2_hse(self, defect_taskid, localisation=0.2, selected_bands=None, edge_tol=(0.5, 0.5)):
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
                None,
                None,
                None,  #(host_db, host_taskid, 0, vbm_dx, cbm_dx),
                localisation,  #0.2
                locpot_c2db=None,  #(c2db, c2db_uid, 0)
                is_vacuum_aligment_on_plot=True,
                edge_tol=edge_tol, # defect state will be picked only if it's above vbm by
                # 0.025 eV
                # and below
                # cbm by 0.025 eV
                ir_db=ir_db,
                ir_entry_filter={"prev_fw_taskid": defect_taskid},
                selected_bands=selected_bands,
            )
            tot, proj, d_df, levels, defect_levels = state
            level_info = d_df.to_dict("records")[0]
        except Exception as er:
            print(er)
            level_info = {}
            levels = {}
            defect_levels = {}
            proj = None
        return level_info, levels, defect_levels, proj


    def get_defect_df_v2_hse(self, read_c2db_uid_key=False):
        data = []
        col = HSEQubitDefect.collection
        entries = list(col.find(self.df_filter))
        for idx, e in enumerate(entries):
            print(e["task_id"])
            c2db_uid = e["pc_from"].split("/")[-1] if not read_c2db_uid_key else e["c2db_uid"]
            e_from_host = C2DB_IR_calc_data.collection.find_one({"c2db_uid": c2db_uid, "task_label": "hse line"})
            host_c2db_info = e_from_host["c2db_info"]
            for field in ["spacegroup", "pmg_point_gp", "irreps", "formula"]:
                if host_c2db_info.get(field):
                    host_c2db_info.pop(field)

            host_band_edges = e_from_host["band_edges"]
            for field in ["vbm_up_proj_on_el", "vbm_up_orbital_proj_on_el",
                          "vbm_down_proj_on_el", "vbm_down_orbital_proj_on_el",
                          "cbm_up_proj_on_el", "cbm_up_orbital_proj_on_el",
                          "cbm_down_proj_on_el", "cbm_down_orbital_proj_on_el"]:
                if host_band_edges.get(field):
                    host_band_edges.pop(field)

            vbm_ir = host_band_edges.get("vbm_up_band_ir") or host_band_edges.get("vbm_down_band_ir")
            cbm_ir = host_band_edges.get("cbm_up_band_ir") or host_band_edges.get("cbm_down_band_ir")
            level_edge_ir = (vbm_ir, cbm_ir)
            host_band_edges.update({"level_edge_ir": level_edge_ir})

            host_sym_data = e_from_host["sym_data"]
            host_sym_data.update({"reduced_site_sym": tuple(e_from_host["sym_data"]["good_ir_info"]["site_sym"]),
                                  "reduced_site_specie": tuple(e_from_host["sym_data"]["good_ir_info"]["species"])
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
                "host_taskid": "{}-{}-{}".format(C2DB_IR_calc_data.db_name, C2DB_IR_calc_data.collection.name,
                                                 e_from_host["task_id"]),
                "gap_hse_c2db_ir": e_from_host["output"]["bandgap"],
                "defect_name": e["defect_entry"]["name"],
                "defect_type": e["defect_entry"]["defect_type"],
                "charge": e["charge_state"],
                "mag": round(e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]),
                "chemsys": e["chemsys"],
                "is_nelect_even": is_nelect_even,
                "site_oxi_state": tuple([tuple(i) for i in e_from_host["site_oxi_state"]]),
                "number_NN": len(e["NN"]),
                "nbands": e["input"]["parameters"]["NBANDS"]
            }
            for host_info in [host_c2db_info, host_band_edges, host_sym_data]:
                info.update(host_info)

           # if e["task_id"] in [575]:
            #     localisation = 0.2
            # elif e["task_id"] in [1088, 2590, 585, 545, 571, 2563, 603, 644, ]:
            #     localisation = 0.1
            # elif e["task_id"] in [605, 2569]:
            #     localisation = 0.05
            # else:
            #     localisation = 0.2
            def settings_for_36group(task_id):
                localisation, selected_bands, edge_tol = 0.2, None, (0.5, 0.5)
                if task_id in [238]:
                    localisation = 0.25
                elif task_id in [220]:
                    localisation = 0.114
                elif task_id in [1088]:
                    localisation = 0.1
                elif task_id in [585]:
                    localisation = 0.08
                    selected_bands = {"1": [306, 307, 308], "-1": [306, 307, 309]}
                elif task_id in [2576]:
                    localisation = 0.02
                    selected_bands = {"1": [306, 307, 312], "-1": [306, 307, 312]}
                elif task_id in [605]:
                    localisation = 0.11
                    selected_bands = {"1": [306, 307, 309], "-1": [306, 307, 309]}
                elif task_id in [2592]:
                    localisation = 0.05
                    selected_bands = {"1": [306, 307, 313], "-1": [306, 307, 324]}
                    edge_tol = (1, 1)
                elif task_id in [571]:
                    localisation = 0.085
                elif task_id in [545]:
                    localisation = 0.15
                elif task_id in [2569]:
                    localisation = 0.14
                return localisation, selected_bands, edge_tol

            localisation, selected_bands, edge_tol = settings_for_36group(e["task_id"])


            d_df, levels, in_gpa_levels, proj_df = self.extract_defect_levels_v2_hse(
                e["task_id"],
                localisation=localisation,
                selected_bands=selected_bands,
                edge_tol=edge_tol
            )
            info.update({"localisation_threshold": localisation})
            info.update(d_df)
            info.update(levels)
            info.update(in_gpa_levels)
            data.append(info)

        self.defect_df = pd.DataFrame(data)
        self.defect_df.fillna("None", inplace=True)
        self.defect_df.replace({"up_tran_en": "None"}, 0, inplace=True)
        self.defect_df.replace({"dn_tran_en": "None"}, 0, inplace=True)

        self.input_df = self.defect_df
        # IOTools(cwd=save_xlsx_path, pandas_df=self.defect_df).to_excel("taskid_970")

    def backprocess(self, excel_name=None):
        self.add_band_edges_and_defects()
        self.add_level_category()
        self.add_transition_wavevlength()
        self.add_cdft_occupations()
        self.add_hse_fworker()
        if excel_name:
            self.df_to_excel(excel_name=excel_name)

class DataPrepCDFT(CDFT):
    def __init__(self, defect_entry_df):
        super().__init__()
        self.defect_entry_df = defect_entry_df # hse_screened_qubits_df
        self.zpl_df = None

    def get_tgt_df(self):
        self.get_data_sheet({"taskid": {"$in": list(self.defect_entry_df["task_id"])}}, read_c2db_uid_key=True)
        # You need to update the entry in databased for TDM information first!

    def run_TDM(self):
        from wf.wf import transition_dipole_moment
        transition_dipole_moment(self.foundation_df)

    def get_zpl_df(self):
        zpl_df = self.get_zpl_data()
        zpl_df["task_id"] = zpl_df["gs_taskid"]
        zpl_df.drop(columns=["charge", "prototype"], inplace=True)
        self.zpl_df = zpl_df

    def one_shot_calc_TDM(self):
        self.get_tgt_df()
        self.run_TDM()

    def one_shot_zpl_df(self):
        self.get_tgt_df()
        self.get_zpl_df()

def main():
    # define a function that updates c2db_ir_calc_data with the new data
    def update_c2db_ir_calc_data():
        # DataPrepHost.sym_data()
        DataPrepHost.band_edge()
        # DataPrepHost.site_oxi_state()
        # DataPrepHost.cp_c2db_info()

    # define a function to generate a table of defects
    def get_defect_table():
        from analysis.analysis_api import hse_qubit_df
        # antisite_tmd = {"pc_from": {"$regex": "owls"}, "task_label": "HSE_scf",
        #                 "chemsys": {"$in": ["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"]}}

        # hse_qubit = hse_qubit_df.copy()
        # taskid_list = hse_qubit["task_id"].to_list()
        # for i in [173, 275,934]:
        #     taskid_list.remove(i)
        # print(taskid_list)

        tkids = defects_36_groups_df["task_id"].to_list()

        # filter = antisite_tmd
        # DataPrepDefect.cp_symdata_bandedges(filter, read_c2db_uid_key=True)
        # DataPrepDefect.is_site_sym_uniform(filter)
        # DataPrepDefect.cp_site_oxi_state(filter, read_c2db_uid_key=True)

        test = GenerateDefectTable({"task_id": {"$in": tkids}})
        test.get_defect_df_v2_hse(read_c2db_uid_key=False)
        test.backprocess()
        test.df_to_excel(excel_name="defects_36_groups")

    def get_defect_table_scan():
        table1_df_failed = table1_df.loc[(table1_df["up_tran_en"]==0) & (table1_df["dn_tran_en"]==0)]
        tkids = table1_df["task_id"].to_list()[:]

        test = GenerateDefectTable({"task_id": {"$in": tkids}})
        test.get_defect_df_v2()
        test.backprocess()
        test.df_to_excel(excel_name="table1_failed")

    def get_zpl_df():
        # run TDM in db1
        p_path = "/home/qimin/sdb_tsai/site-packages/JPack_independent/projects/defectDB"
        tgt_df = IOTools(excel_file="wte2_mote2",
                         cwd=os.path.join(p_path, INPUT_PATH)).read_excel()
        zpl = DataPrepCDFT(defect_entry_df=tgt_df)
        zpl.one_shot_calc_TDM()
        zpl.one_shot_zpl_df()
        IOTools(pandas_df=zpl.zpl_df, cwd=os.path.join(p_path, "analysis/output/xlsx")).to_excel(
            "test_zpl_df")

    update_c2db_ir_calc_data()

if __name__ == '__main__':
    main()
    # run TDM in db1
    # p_path = "/home/qimin/sdb_tsai/site-packages/JPack_independent/projects/defectDB"
    # tgt_df = IOTools(excel_file="Table_4_df_2022-04-21", cwd=INPUT_PATH).read_excel()
    # zpl = DataPrepCDFT(defect_entry_df=tgt_df.loc[tgt_df["task_id"] == 605])
    # zpl.one_shot_calc_TDM()
    # zpl.one_shot_zpl_df()
    # zpl.get_zpl_df()
    # a = zpl.zpl_df
    # IOTools(pandas_df=zpl.zpl_df, cwd="analysis/output/xlsx").to_excel(
    #     "test_zpl_df")