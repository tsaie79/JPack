import os, datetime

from qubitPack.tool_box import get_db
from qubitPack.tool_box import IOTools

import pandas as pd

from pymatgen import Spin

from monty.serialization import loadfn

C2DB = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan", port=1234)
SCAN2dMat = get_db("Scan2dMat", "calc_data",  user="Jeng_ro", password="qimin", port=1236)
SCAN2dDefect = get_db("Scan2dDefect", "calc_data",  user="Jeng_ro", password="qimin", port=1236)
C2DB_IR = get_db("C2DB_IR", "calc_data",  user="Jeng_ro", password="qimin", port=1234)

PATH_EXCEL = "searching_triplet/xlsx/"

class ReadHosts:
    def __init__(self):
        self.host_df = IOTools(excel_file="good_scan_pg_2021-09-02").read_excel()

    def screen_hosts(self):
        condition = self.host_df["is_vbm_cbm_from_same_element"] == 1
        good_hosts = self.host_df.loc[condition, :]
        df_good_host_ids = good_hosts.loc[:, ["uid", "tk_id"]]
        return df_good_host_ids

class ReadDefects:
    def __init__(self):
        self.good_hosts_ids = IOTools(json_file="good_hosts_ids_2021-09-02").read_json()
        self.data = []
    def screen_defects(self):
        col = SCAN2dDefect.collection
        for good_host_id in self.good_hosts_ids:
            c2db_uid = good_host_id["uid"]
            task_id = good_host_id["tk_id"]
            for e in col.find({"host_info.c2db_info.uid": c2db_uid, "task_label": "SCAN_nscf uniform"}):
                print(e["host_info"]["c2db_info"]["uid"], e["task_id"])
                info = {
                    "task_id": e["task_id"],
                    "c2db_uid": e["host_info"]["c2db_info"]["uid"],
                    "task_label": e["task_label"],
                    "defect_name": e["defect_name"],
                    "charge": e["charge_state"],
                    "mag": e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]
                }
                self.data.append(info)
        return pd.DataFrame(self.data)

    @classmethod
    def run_screening(cls):
        d = cls().screen_defects()
        IOTools(output_path=PATH_EXCEL, pandas_df=d).to_excel("good_defects_nscf")
        IOTools(output_path=PATH_EXCEL, pandas_df=d).to_json("good_defects_nscf")


    def update_defect_entry(self):
        db = SCAN2dDefect
        for e in self.data[:1]:
            db.collection.update_one({"task_id": e["task_id"]}, {"$set": {"bandedges_uniform": True}})


if __name__ == '__main__':
    ReadDefects.run_screening()


