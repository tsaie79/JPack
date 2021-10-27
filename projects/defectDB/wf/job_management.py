import json
import os, shutil, subprocess
from glob import glob
import monty.serialization
from datetime import datetime
import pandas as pd

from qubitPack.tool_box import get_db, IOTools

from pymatgen.io.vasp.inputs import Structure, Incar
from pymatgen.io.vasp.outputs import Oszicar
from fireworks import LaunchPad

db = get_db("Scan2dDefect", "calc_data", port=12347)
wf_col = get_db("Scan2dDefect", "workflows", port=12347).collection

class CheckJobs:
    def __init__(self):
        self.lpad = LaunchPad(host=db.host, port=db.port, name=db.db_name, username=db.user, password=db.password)
        self.states = ["COMPLETED", "FIZZLED", "RUNNING", "READY", "RESERVED", "ALL"]
        self.data = []
    @classmethod
    def completed(cls, filter):
        new_filter = filter.copy()
        new_filter.update({"state": "COMPLETED"})
        return new_filter

    @classmethod
    def fizzled(cls, filter):
        new_filter = filter.copy()
        new_filter.update({"state": "FIZZLED"})
        return new_filter

    @classmethod
    def running(cls, filter):
        new_filter = filter.copy()
        new_filter.update({"state": "RUNNING"})
        return new_filter

    @classmethod
    def ready(cls, filter):
        new_filter = filter.copy()
        new_filter.update({"state": "READY"})
        return new_filter

    @classmethod
    def reserved(cls, filter):
        new_filter = filter.copy()
        new_filter.update({"state": "RESERVED"})
        return new_filter

    @classmethod
    def all(cls, filter):
        return filter

    def batch_1(self):
        batch_ids = [13, 14, 15]
        filter = {
                "metadata.tags.0.group_id": {"$exists":0},
                "created_on": {"$gte": datetime(2021,8,1), "$lt": datetime(2021, 8,2)}
            }
        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            print(f(filter))
            wfs = self.lpad.get_wf_ids(f(filter))
            wf_names = [wf["name"] for wf in wf_col.find(f(filter))]
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs,
                              "wf_names": wf_names})

    def batch_2(self):
        batch_ids = [12]
        filter = {
                # "metadata.tags.0.group_id": {"$exists":0},
                "created_on": {"$gte": datetime(2021,7,31), "$lt": datetime(2021, 8,1)}
            }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            print(f(filter))
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})


    def batch_3(self):
        batch_ids = [10]
        filter = {
                "metadata.tags.0.group_id": {"$exists":0},
                "created_on": {"$gte": datetime(2021,7,29), "$lt": datetime(2021, 7,30)}
            }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            print(f(filter))
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_4(self):
        batch_ids = [8]
        filter = {
                "metadata.tags.0.group_id": {"$exists":0},
                "created_on": {"$gte": datetime(2021,7,28), "$lt": datetime(2021, 7,29)}
            }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_5(self):
        charge_states = [0, 1]
        batch_ids = [9, 11]
        for chg in charge_states:
            filter = {"metadata.tags.0.group_id":4 , "metadata.tags.0.charge_state.0":chg}
        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
                wfs = self.lpad.get_wf_ids(f(filter))
                print(wfs)
                print("{}: {}".format(chg, len(wfs)))
                self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "charge": chg, "wfs": wfs})


    def batch_6(self):
        batch_ids = [1]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [0],
            "created_on": {"$gte": datetime(2021,8, 20), "$lt": datetime(2021, 8, 21)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_7(self):
        batch_ids = [5]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [0],
            "created_on": {"$gte": datetime(2021,8, 23), "$lt": datetime(2021, 8, 24)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_8(self):
        batch_ids = [4]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [0],
            "created_on": {"$gte": datetime(2021,9, 3), "$lt": datetime(2021, 9, 4)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_9(self):
        batch_ids = [0]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [0],
            "created_on": {"$gte": datetime(2021,9, 5), "$lt": datetime(2021, 9, 6)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_10(self):
        batch_ids = [0]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [-1],
            "created_on": {"$gte": datetime(2021,10,11), "$lte": datetime(2021, 10, 16)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            wf_names = [wf["name"] for wf in wf_col.find(f(filter))]
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs,
                              "wf_names": wf_names})
        df = pd.DataFrame(self.data)
        return df

    def batch_11(self):
        batch_ids = [0]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [1],
            "created_on": {"$gte": datetime(2021,10,17), "$lt": datetime(2021, 10, 18)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            wf_names = [wf["name"] for wf in wf_col.find(f(filter))]
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs,
                              "wf_names": wf_names})
        df = pd.DataFrame(self.data)
        return df

    def batch_12(self):
        batch_ids = [0]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [-1],
            "created_on": {"$gte": datetime(2021,10,20), "$lt": datetime(2021, 10, 21)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            wf_names = [wf["name"] for wf in wf_col.find(f(filter))]
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs,
                              "wf_names": wf_names})
        df = pd.DataFrame(self.data)
        return df

    def batch_13(self):
        batch_ids = [0]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [1],
            "created_on": {"$gte": datetime(2021,10,23), "$lt": datetime(2021, 10, 24)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            wf_names = [wf["name"] for wf in wf_col.find(f(filter))]
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs,
                              "wf_names": wf_names})
        df = pd.DataFrame(self.data)
        return df

    @classmethod
    def test(cls):
        batch_10 = cls().batch_10()
        batch_11 = cls().batch_11()
        batch_10_list = batch_10.loc[batch_10["state"] == "ALL", "wf_names"].tolist()[0]
        batch_11_list = batch_11.loc[batch_11["state"] == "ALL", "wf_names"].tolist()[0]

        # batch_10_list = [wf.split(":")[0] for wf in batch_10_list]
        # batch_11_list= [wf.split(":")[0] for wf in batch_11_list]

        print(len(dict.fromkeys(batch_10_list)))
        print(len(dict.fromkeys(batch_11_list)))

        # diff = []
        # inter = []
        # for wf_name in batch_11_list:
        #     if wf_name not in batch_10_list:
        #         print(wf_name)
        #         diff.append(wf_name)
        #     if wf_name in batch_10_list:
        #         print(wf_name)
        #         inter.append(wf_name)
        # print(len(diff), len(inter))

    @classmethod
    def run(cls):
        chk = cls()
        return  chk.batch_0()
        # df = pd.DataFrame(chk.data)
        # df.to_clipboard()
        # print(df)
        # return df

class RerunJobs:
    def __init__(self, json_file):
        self.lpad = LaunchPad(host=db.host, port=db.port, name=db.db_name, username=db.user, password=db.password)
        self.fws_df = pd.DataFrame(monty.serialization.loadfn("wf/unfinished_jobs/{}".format(json_file)))

    def rerun_fw(self):
        """
        fworker = "nersc", "owls", "efrc"
        """
        def delet_dir(dir_name):
            shutil.rmtree(dir_name)

        def rerun_opt(fw_id):
            prev_path = self.lpad.get_launchdir(fw_id, 0)
            if prev_path:
                print(fw_id, prev_path)
                os.chdir(prev_path)
                try:
                    for f_gz in glob("*.gz"):
                        subprocess.call("gunzip {}.gz".format(f_gz).split(" "))
                except Exception:
                    pass
                if glob("*") == []:
                    self.lpad.rerun_fw(fw_id)
                elif glob("CONTCAR.relax*") != []:
                    shutil.copy("CONTCAR", "POSCAR")
                    self.lpad.rerun_fw(fw_id, recover_launch="last", recover_mode="prev_dir")
                else:
                    shutil.copy("CONTCAR", "POSCAR")
                    self.lpad.rerun_fw(fw_id, recover_launch="last", recover_mode="prev_dir")
            else:
                self.lpad.rerun_fw(fw_id)

        def rerun_scf(fw_id):
            subprocess.call("gunzip INCAR.gz".split(" "))
            incar = Incar.from_file("INCAR")
            incar.update({"NELM": 200, "ALGO": "Fast"})
            incar.write_file("INCAR")
            lpad.rerun_fw(fw_id, recover_launch="last", recover_mode="prev_dir")

        def run():
            conditions = (self.fws_df["state"].isin(["RUNNING", "FIZZLED"])) & \
                         (self.fws_df["charge_state"] == 0) &\
                         (self.fws_df["name"].str.contains("relax")) & \
                         (self.fws_df["fworker"] == "owls")

            fw_ids = self.fws_df.loc[conditions, ["fw_id", "state"]]
            for idx, fw_status in fw_ids.iterrows():
                if self.lpad.detect_lostruns(query={"fw_id": fw_status["fw_id"]})[1] == [] \
                        and fw_status["state"] == "RUNNING":
                    continue
                try:
                    print(fw_status["fw_id"])
                    # delet_dir(prev_path)
                    rerun_opt(fw_status["fw_id"])
                    # rerun_scf(fw_id)
                except Exception as err:
                    print(err)
                    continue
        run()
class ResetFworker:
    def __init__(self):
        self.lpad = LaunchPad(host=db.host, port=db.port, name=db.db_name, username=db.user, password=db.password)
        self.init_fws = None

    def locate_jobs(self):
        self.init_fws = self.lpad.get_fw_ids({"state": "READY", "spec._fworker":"owls",
                                              "name": {"$regex": "SCAN_relax"}})
        print("Number of wfs: {}".format(len(self.init_fws)))

    def change_fworker(self, fraction, new_fworker):
        total_fws = len(self.init_fws)
        print("Number of changing wfs: {}".format(int(total_fws*fraction)))
        for init_fw in self.init_fws[:int(total_fws*fraction)]:
            wf = self.lpad.get_wf_by_fw_id(init_fw)
            fw_ids = list(wf.id_fw.keys())
            print(fw_ids)
            self.lpad.update_spec(fw_ids, {"_fworker": new_fworker})

    @classmethod
    def run(cls):
        reset_fworker = cls()
        reset_fworker.locate_jobs()
        print(reset_fworker.init_fws)
        reset_fworker.change_fworker(1/2, "efrc")


class DeleteCompleteJobs:
    def __init__(self, fworker):
        self.lpad = LaunchPad(host=db.host, port=db.port, name=db.db_name, username=db.user, password=db.password)
        self.fworker = fworker
        self.init_fws = []
    def locate_completed_wf(self):
        # wf_ids = self.lpad.get_wf_ids({"state": "COMPLETED"})
        self.init_fws = self.lpad.get_fw_ids({"state": "COMPLETED", "spec._fworker":self.fworker,
                                              "name": {"$regex": "irvsp"}})

        print("Number of wfs: {}".format(len(self.init_fws)))

    def delete_dirs_completed_wf(self):
        for init_fw in self.init_fws[120:]:
            wf = self.lpad.db["workflows"].find_one({"nodes": {"$elemMatch": {"$eq": init_fw}}})
            fw_ids = list(wf["nodes"])
            print(fw_ids)
            for fw_id in fw_ids:
                l_dir = self.lpad.db["launches"].find_one({"state": "COMPLETED", "fw_id": fw_id})["launch_dir"]
                print(l_dir)
                try:
                    shutil.rmtree(l_dir)
                except FileNotFoundError:
                    continue
    @classmethod
    def run(cls):
        op = cls("efrc")
        op.locate_completed_wf()
        op.delete_dirs_completed_wf()

class RemainingJobs:
    def __init__(self):
        self.lpad = LaunchPad(host=db.host, port=db.port, name=db.db_name, username=db.user, password=db.password)

        self.save_xlsx_path = "/Users/jeng-yuantsai/Research/project/Scan2dDefect/latest_results/xlsx"
        self.save_plt_path = "/Users/jeng-yuantsai/Research/project/Scan2dDefect/latest_results/plt"

        self.defect_df = IOTools(cwd=self.save_xlsx_path, excel_file="defect_0_2_2021-10-22").read_excel()
        self.host_df = IOTools(cwd=self.save_xlsx_path, excel_file="host_2021-10-25").read_excel()

        self.uni = self.host_df["reduced_site_sym"].isin([("-6m2", "-6m2"), ("3m", "3m"), ("4mm", "4mm")])
        self.uniform_host_df = self.host_df.loc[self.uni, :]

        self.nonuni = ~self.host_df["reduced_site_sym"].isin([("-6m2", "-6m2"), ("3m", "3m"), ("4mm", "4mm")])
        self.nonuniform_host_df = self.host_df.loc[self.nonuni, :]

        self.uniform_host_gp = self.uniform_host_df.groupby(["spg_number", "reduced_site_sym"]).agg(
            dict(taskid=["count", "unique"], c2db_uid=["count", "unique"], spg_number=["count", "unique"]))
        self.nonuniform_host_gp = self.nonuniform_host_df.groupby(["spg_number", "reduced_site_sym"]).agg(
            dict(taskid=["count", "unique"], c2db_uid=["count", "unique"], spg_number=["count", "unique"]))


    def nonuni_not_calculated(self):
        c = self.nonuniform_host_df["c2db_uid"].isin(self.defect_df["c2db_uid"])
        result_gp = self.nonuniform_host_df.loc[~c, ["spg_number", "taskid", "c2db_uid"]].groupby(["spg_number"]).agg(dict(
            taskid=["count", "unique"], c2db_uid=["count", "unique"]))
        return  result_gp


    def uni_not_calculated(self):
        c = self.uniform_host_df["c2db_uid"].isin(self.defect_df["c2db_uid"])
        result_gp = self.uniform_host_df.loc[~c, ["spg_number", "taskid", "c2db_uid"]].groupby(["spg_number"]).agg(dict(
            taskid=["count", "unique"], c2db_uid=["count", "unique"]))
        return  result_gp

    def existed_fw_status(self):
        job_states = ["FIZZLED", "RESERVED", "RUNNING"]
        def get_existed_fw_status(c2db_uids, excel_name):
            data = []
            for c2db_uid in c2db_uids[:]:
                nscf_fw_ids = self.lpad.get_fw_ids({"spec._tasks.10.additional_fields.host_info.c2db_info.uid": c2db_uid})
                for nscf_fw_id in nscf_fw_ids:
                    wf = self.lpad.db["workflows"].find_one({"nodes": {"$in": [nscf_fw_id]}})
                    for fw_id in wf["nodes"]:
                        print(fw_id)
                        entry = {}
                        fw = self.lpad.db["fireworks"].find_one({"fw_id": fw_id})

                        todb_firetask = None
                        for ft in fw["spec"]["_tasks"]:
                            if "ToDb" in ft["_fw_name"]:
                                todb_firetask = ft

                        additional_fields = todb_firetask["additional_fields"]
                        charge_state = additional_fields["charge_state"]
                        c2db_uid = additional_fields["host_info"]["c2db_info"]["uid"]
                        host_taskid = additional_fields["pc_from_id"]
                        if "irvsp" in fw["name"]:
                            defect_name = additional_fields["defect_name"]
                            try:
                                spg_number = additional_fields["host_info"]["sym_data"]["number"]
                            except KeyError:
                                spg_number = additional_fields["host_info"]["sym_data"]["pmg_spg_number"]

                        else:
                            defect_name = additional_fields["defect_entry"]["name"]
                            try:
                                spg_number = additional_fields["host_info"]["sym_data"]["number"]
                            except KeyError:
                                spg_number = additional_fields["host_info"]["sym_data"]["pmg_spg_number"]

                        entry["c2db_uid"] = c2db_uid
                        entry["host_taskid"] = host_taskid
                        entry["defect_name"] = defect_name
                        entry["charge_state"] = charge_state
                        entry["spg_number"] = spg_number
                        entry["fw_id"] = fw["fw_id"]
                        entry["fworker"] = fw["spec"]["_fworker"]
                        entry["state"] = fw["state"]
                        entry["name"] = fw["name"]
                        data.append(entry)

            df = pd.DataFrame(data).drop_duplicates()
            df = df.loc[(df["state"].isin(job_states)), :]
            # IOTools(cwd=os.path.join(self.save_xlsx_path, "jobs/existed_fws_status"), pandas_df=df).to_excel(
            #     excel_name, index=True)
            IOTools(cwd="wf/unfinished_jobs", pandas_df=df).to_json(excel_name, index=True)
            return df

        def group_df(df, excel_name):
            df = df.groupby(["spg_number", "host_taskid", "c2db_uid", "defect_name", "charge_state", "name",
                             "fworker"]).agg(
                dict(state=["unique"], fw_id=["unique", "count"]))
            output = df.style\
                .applymap(lambda x: "background-color:red" if x[0] == "FIZZLED" else "", subset=["state"])\
                .applymap(lambda x: "background-color:pink" if x[0] == "RESERVED" else "", subset=["state"])\
                .applymap(lambda x: "background-color:yellow" if x[0] == "RUNNING" else "", subset=["state"])

            IOTools(cwd=os.path.join(self.save_xlsx_path, "jobs/existed_fws_status"),
                    pandas_df=output).to_excel(excel_name, index=True)

        def get_existed_fw_status_run():
            uni_host_uids, nonuni_host_uids = self.uniform_host_df["c2db_uid"].unique(), self.nonuniform_host_df[
                "c2db_uid"].unique()
            for host_sym_type, name in zip([uni_host_uids, nonuni_host_uids], ["uni", "nonuni"]):
                get_existed_fw_status(host_sym_type, name)

        def group_df_run():
            for i, o in zip(
                    ["uni_fizzled_2021-10-27", "nonuni_fizzled_2021-10-27"],
                    ["uni_fizzled_gp", "nonuni_fizzled_gp"]):
                df = IOTools(cwd=os.path.join(self.save_xlsx_path, "jobs/existed_fws_status"),
                             excel_file=i).read_excel()
                group_df(df, excel_name=o)

        def run():
            uni_host_uids, nonuni_host_uids = self.uniform_host_df["c2db_uid"].unique(), self.nonuniform_host_df[
                "c2db_uid"].unique()
            for host_sym_type, name in zip([uni_host_uids, nonuni_host_uids], ["uni", "nonuni"]):
                df = get_existed_fw_status(host_sym_type, name+"_owls")
                group_df(df, excel_name=name+"_owls_gp")

        get_existed_fw_status_run()

def local_main():
    RemainingJobs().existed_fw_status()

def remote_main():
    rerun_list = RerunJobs("nonuni_owls_2021-10-27.json")
    rerun_list.rerun_fw()


if __name__ == '__main__':
    local_main()
    # d = CheckJobs.run()
    # d = d.groupby(["c2db_uid", "defect_name", "charge_state"])
    # a = RerunJobs("nonuni_owls_2021-10-27.json")
    # a.rerun_fw()
    # a.get_fw_ids_from_json()
    # a.rerun_fw()
    # RemainingJobs().existed_fw_status()
    # uni = RemainginJobs().uni_not_calculated()
    # nonuni = RemainginJobs().nonuni_not_calculated()



