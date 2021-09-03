from qubitPack.tool_box import get_db
from fireworks import LaunchPad
from datetime import datetime
import pandas as pd


class CheckJobs:
    def __init__(self):
        db = get_db("Scan2dDefect", "calc_data", port=1236)
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
        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.all], self.states):
            print(f(filter))
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_2(self):
        batch_ids = [12]
        filter = {
                "metadata.tags.0.group_id": {"$exists":0},
                "created_on": {"$gte": datetime(2021,7,31), "$lt": datetime(2021, 8,1)}
            }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.all], self.states):
            print(f(filter))
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})


    def batch_3(self):
        batch_ids = [10]
        filter = {
                "metadata.tags.0.group_id": {"$exists":0},
                "created_on": {"$gte": datetime(2021,7,29), "$lt": datetime(2021, 7,30)}
            }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.all], self.states):
            print(f(filter))
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_4(self):
        batch_ids = [8]
        filter = {
                "metadata.tags.0.group_id": {"$exists":0},
                "created_on": {"$gte": datetime(2021,7,28), "$lt": datetime(2021, 7,29)}
            }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_5(self):
        charge_states = [0, 1]
        batch_ids = [9, 11]
        for chg in charge_states:
            filter = {"metadata.tags.0.group_id":4 , "metadata.tags.0.charge_state.0":chg}
            for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.all], self.states):
                wfs = self.lpad.get_wf_ids(f(filter))
                print(wfs)
                print("{}: {}".format(chg, len(wfs)))
                self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "charge": chg, "wfs": wfs})


    def batch_6(self):
        batch_ids = [0, 1]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [0],
            "created_on": {"$gte": datetime(2021,8, 20), "$lt": datetime(2021, 8, 21)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

    def batch_7(self):
        batch_ids = [4, 5]
        filter = {
            "metadata.tags.0.group_id": {"$exists":1},
            "metadata.tags.0.charge_state": [0],
            "created_on": {"$gte": datetime(2021,8, 23), "$lt": datetime(2021, 8, 24)},
        }

        for f, state in zip([self.completed, self.fizzled, self.running, self.ready, self.reserved, self.all], self.states):
            wfs = self.lpad.get_wf_ids(f(filter))
            self.data.append({"state": state, "number": len(wfs), "batch_ids": batch_ids, "wfs": wfs})

class RerunJobs:
    def rerun_fw(fworker):
        """
        fworker = "nersc", "owls", "efrc"
        """
        from fireworks import LaunchPad
        import os, subprocess, shutil, datetime
        from glob import glob
        from pymatgen.io.vasp.inputs import Structure, Incar
        from pymatgen.io.vasp.outputs import Oszicar
        from qubitPack.tool_box import get_db

        db = get_db("mgb2", "workflows", port=12345, user="Jeng")
        lpad = LaunchPad(host=db.host, port=db.port, name=db.db_name, username=db.user, password=db.password)

        def delet_dir(dir_name):
            shutil.rmtree(dir_name)
        def rerun_opt(fw_id): #421
            # for f in glob("POSCAR*"):
            #     os.remove(f)
            lpad.rerun_fw(fw_id, recover_launch="last", recover_mode="prev_dir")
            fw = lpad.get_fw_by_id(fw_id)
            # fw.spec['_queueadapter'] = {'walltime': '06:00:00'}
            # fw.tasks[1]["other_params"]["user_incar_settings"].update({"ALGO":"All"})
            # fw.tasks[1]["incar_update"].update({"ALGO": "All"})
            # lpad.update_spec([fw.fw_id], fw.as_dict()["spec"])
            shutil.copy("CONTCAR", "POSCAR")
            # subprocess.call("qlaunch -c {} -r singleshot -f {}".format(
            #     os.path.join(os.path.expanduser("~"), "config/project/Scan2dDefect/calc_data/"), fw.fw_id).split(" "))
        def rerun_scf(fw_id):
            subprocess.call("gunzip INCAR.gz".split(" "))
            incar = Incar.from_file("INCAR")
            incar.update({"NELM": 200, "ALGO": "Fast"})
            incar.write_file("INCAR")
            lpad.rerun_fw(fw_id, recover_launch="last", recover_mode="prev_dir")




        fws = lpad.get_fw_ids(
            {
                "state": {"$in": ["READY"]}, #First running, then fizzled
                "name": {"$regex": "static"},
                "spec._fworker": fworker,
                "created_on": {"$gte": str(datetime.datetime(2021,8,12))}
                # "fw_id": {"$gte": 1797}
                # "fw_id": 2550
            }
        )
        print(fws)
        a = []
        for fw_id in fws[::-1][:]:
            prev_fw = lpad.get_fw_by_id(fw_id)
            fworker = prev_fw.spec["_fworker"]
            prev_path = lpad.get_launchdir(fw_id, 0)
            os.chdir(prev_path)
            # a.append(fw_id)
            # if glob("CONTCAR.relax*") != []:
            try:
                print(fw_id, fworker, prev_path, prev_fw.name)
                # delet_dir(prev_path)
                # rerun_opt(fw_id)
                # rerun_scf(fw_id)
            except Exception as err:
                print(err)
                continue
            # try:
            #     oszicar = Oszicar("OSZICAR")
            # except Exception as err:
            #     print(fw_id, err)
            #     continue
            # if len(Oszicar("OSZICAR").ionic_steps) > 10:
            #     a.append(fw_id)
            #     print(fw_id, fworker, path, fw.name)
            #     # rerun_opt(fw)
        print(a, len(a))

#%%
class ResetFworker:
    def __init__(self):
        db = get_db("Scan2dDefect", "calc_data", port=1236, user="Jeng")
        self.lpad = LaunchPad(host=db.host, port=db.port, name=db.db_name, username=db.user, password=db.password)
        self.init_fws = None

    def locate_jobs(self):
        self.init_fws = self.lpad.get_fw_ids({"state": "READY", "spec._fworker":"owls", "name": {"$regex": "SCAN_relax"}})
        print("Number of wfs: {}".format(len(self.init_fws)))

    def change_fworker(self):
        for init_fw in self.init_fws[:]:
            wf = self.lpad.get_wf_by_fw_id(init_fw)
            fw_ids = list(wf.id_fw.keys())
            print(fw_ids)
            self.lpad.update_spec(fw_ids, {"_fworker": "efrc"})



a = ResetFworker()
a.locate_jobs()
a.change_fworker()
#%%
def main():
    chk = CheckJobs()
    # chk.batch_1()
    # chk.batch_2()
    # chk.batch_3()
    # chk.batch_4()
    chk.batch_7()
    df = pd.DataFrame(chk.data)
    df.to_clipboard()
    return df

if __name__ == '__main__':
    df = main()
    print(df)



