from analysis.data_analysis import BackProcess, HSEQubitDefect, HSEQubitIR
from qubitPack.tool_box import get_db, get_good_ir_sites, get_unique_sites_from_wy, get_band_edges_characters

from pymatgen.io.vasp.inputs import Structure

from monty.json import jsanitize
import pandas as pd

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
        filter = {"task_label": "hse line"}

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
    def cp_symdata_bandedges(cls, filter, c2db_uids=None):
        # 4
        from qubitPack.tool_box import get_db, get_band_edges_characters
        from monty.json import jsanitize

        src = C2DB_IR_calc_data
        tgt = HSEQubitDefect

        col = tgt.collection
        entries = list(col.find(filter))
        for idx, e in enumerate(entries):
            c2db_uid = e["pc_from"].split("/")[-1] if c2db_uids is None else c2db_uids[idx]
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
    def cp_site_oxi_state(cls, filter, c2db_uids=None):
        # 8
        from monty.json import jsanitize
        from qubitPack.tool_box import get_db
        from pymatgen import Structure

        host = C2DB_IR_calc_data
        db = HSEQubitDefect

        col = db.collection
        entries = list(col.find(filter))
        for idx, e in enumerate(entries):
            c2db_uid = e["pc_from"].split("/")[-1] if c2db_uids is None else c2db_uids[idx]
            print(e["task_id"], e["pc_from"])
            host_entry = host.collection.find_one({"c2db_uid": c2db_uid, "task_label": "hse line"})
            site_oxi_state = host_entry["site_oxi_state"]
            db.collection.update_one({"task_id": e["task_id"]}, {"$set": {"host_info.c2db_ir_hse_line.site_oxi_state":
                                                                              site_oxi_state}})


class GenerateDefectTable(BackProcess):
    def __init__(self, df_filter):
        super(GenerateDefectTable, self).__init__(None)
        self.defect_df = None
        self.df_filter = df_filter


    def extract_defect_levels_v2_hse(self, defect_taskid, localisation=0.2):
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
                None, #"all",
                None,
                None,  #(host_db, host_taskid, 0, vbm_dx, cbm_dx),
                localisation,  #0.2
                locpot_c2db=None,  #(c2db, c2db_uid, 0)
                is_vacuum_aligment_on_plot=True,
                edge_tol= (0, 0), # defect state will be picked only if it's above vbm by
                # 0.025 eV
                # and below
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


    def get_defect_df_v2_hse(self, c2db_uids=None):
        data = []
        col = HSEQubitDefect.collection
        entries = list(col.find(self.df_filter))
        for idx, e in enumerate(entries):
            print(e["task_id"])
            c2db_uid = e["pc_from"].split("/")[-1] if c2db_uids is None else c2db_uids[idx]
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
                "gap_hse_nosoc": e_from_host["output"]["bandgap"],
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

            localisation = None
            if e["task_id"] in [575]:
                localisation = 0.35
            elif e["task_id"] in [1088, 2590, 585, 545, 571, 2563, 603, 644, ]:
                localisation = 0.1
            elif e["task_id"] in [605, 2569]:
                localisation = 0.05
            else:
                localisation = 0.2
            d_df, levels, in_gpa_levels = self.extract_defect_levels_v2_hse(e["task_id"], localisation=localisation)
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
        # IOTools(cwd=save_xlsx_path, pandas_df=self.defect_df).to_excel("defects_36_groups")

    def backprocess(self):
        self.add_band_edges_and_defects()
        self.add_level_category()
        self.add_transition_wavevlength()
        self.add_cdft_occupations()
        self.add_hse_fworker()
        self.df_to_excel(excel_name="test")


def main():
    # define a function that updates c2db_ir_calc_data with the new data
    def update_c2db_ir_calc_data():
        DataPrepHost.sym_data()
        DataPrepHost.band_edge()
        DataPrepHost.site_oxi_state()
        DataPrepHost.cp_c2db_info()

    # define a function to generate a table of defects
    def get_defect_table():
        antisite_tmd = {
            "task_label": "HSE_scf", "nupdown_set": 2,
            "pc_from":
                {
                    "$in": ["owls/mx2_antisite_pc/{}".format(taskid) for taskid in [3102]]
                }
        }
        # DataPrepDefect.cp_symdata_bandedges(antisite_tmd, c2db_uids=["MoTe2-MoS2-NM"])
        # DataPrepDefect.is_site_sym_uniform(antisite_tmd)
        # DataPrepDefect.cp_site_oxi_state(antisite_tmd, c2db_uids=["MoTe2-MoS2-NM"])

        test = GenerateDefectTable(antisite_tmd)
        test.get_defect_df_v2_hse(c2db_uids=["MoTe2-MoS2-NM"])
        test.backprocess()
    get_defect_table()


if __name__ == '__main__':

    main()
