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

SCAN2dMat = get_db("Scan2dMat", "calc_data",  user="Jeng_ro", password="qimin", port=27017)
SCAN2dDefect = get_db("Scan2dDefect", "calc_data",  user="Jeng_ro", password="qimin", port=27017)
SCAN2dIR = get_db("Scan2dDefect", "ir_data", port=27017)

save_xlsx_path = "output/xlsx"
save_plt_path = "output/plt"
#%%
class Host:
    def __init__(self, host_xlsx="host_2021-10-25"):
        if host_xlsx:
            self.host_df = IOTools(cwd=save_xlsx_path, excel_file=host_xlsx).read_excel()
    @classmethod
    def host_df(cls):
        pass

Host().host_df()
#%%
class Defect:
    @classmethod
    def defect_df(cls):
        data = []
        col = SCAN2dDefect.collection
        for e in col.find({"task_label": "SCAN_scf"}): #"host_info.c2db_info.prototype": {"$in": ["AgBr3"]}
            #     print(e["host_info"]["c2db_info"]["uid"], e["task_id"])
            try:
                is_vbm_cbm_same = host_df.loc[host_df["c2db_uid"] == e["host_info"]["c2db_info"]["uid"],
                                              "is_vbm_cbm_from_same_element"].iloc[0]
                hse_gap_nosoc = host_df.loc[host_df["c2db_uid"] == e["host_info"]["c2db_info"]["uid"], "c2db_gap_hse_nosoc"].iloc[0]
                host_spg_number = host_df.loc[host_df["c2db_uid"] == e["host_info"]["c2db_info"]["uid"], "spg_number"].iloc[0]
                host_spg =  host_df.loc[host_df["c2db_uid"] == e["host_info"]["c2db_info"]["uid"], "spg"].iloc[0]
                host_taskid =  host_df.loc[host_df["c2db_uid"] == e["host_info"]["c2db_info"]["uid"], "taskid"].iloc[0]

            except:
                is_vbm_cbm_same = None
                hse_gap_nosoc = None

            is_nelect_even = None
            if e["input"]["incar"]["NELECT"] % 2 == 0:
                is_nelect_even = True
            else:
                is_nelect_even = False

            host_pot_a, defect_pot_a = Potential(e["task_id"], SCAN2dDefect, SCAN2dMat).linear_fit_potential()

            info = {
                "task_id": e["task_id"],
                "is_vbm_cbm_same": is_vbm_cbm_same,
                "c2db_uid": e["host_info"]["c2db_info"]["uid"],
                "prototype": e["host_info"]["c2db_info"]["prototype"],
                "host_spg": host_spg,
                "host_spg_number": host_spg_number,
                "host_taskid": host_taskid,
                "gap_scan": e["host_info"]["scan_bs"]["bandgap"],
                "hse_gap_nosoc":hse_gap_nosoc,
                "defect_name": e["defect_name"],
                "defect_type": e["defect_entry"]["defect_type"],
                "charge": e["charge_state"],
                "mag": e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"],
                "reduced_site_sym": tuple(e["host_info"]["sym_data"]["good_ir_info"]["site_sym"]),
                "reduced_site_specie": tuple(e["host_info"]["sym_data"]["good_ir_info"]["species"]),
                "chemsys": e["chemsys"],
                "is_nelect_even": is_nelect_even,
                "is_host_pot_steep": abs(host_pot_a[0]) > 5e-4,
                "is_defect_pot_steep": abs(defect_pot_a[0]) > 5e-4,
                "host_pot_a": host_pot_a[0],
                "defect_pot_a": defect_pot_a[0]
            }

            level_info, levels = extract_defect_levels(e["task_id"])
            info.update(level_info)
            info.update(levels)
            data.append(info)


        defect_df = pd.DataFrame(data)
        defect_df.fillna("None", inplace=True)
        display(defect_df)
        IOTools(cwd=save_xlsx_path, pandas_df=defect_df).to_excel("defect")