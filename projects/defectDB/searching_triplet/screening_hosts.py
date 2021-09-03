#%% looking for candidate hosts
import datetime

from qubitPack.tool_box import *

from mpinterfaces.utils import *

import pandas as pd

from pymatgen import Spin

from JPack.projects.defectDB.searching_triplet.screening_defects import IOTools

C2DB = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan", port=1234)
SCAN2dMat = get_db("Scan2dMat", "calc_data",  user="Jeng_ro", password="qimin", port=1236)
SCAN2dDefect = get_db("Scan2dDefect", "calc_data",  user="Jeng_ro", password="qimin", port=1236)
C2DB_IR = get_db("C2DB_IR", "calc_data",  user="Jeng_ro", password="qimin", port=1234)

class DefectStats:
    @classmethod
    def host_candidate_stats(cls):
        col = SCAN2dMat.collection
        filter = {
            "c2db_info.gap_hse_nosoc":{"$gte":1.5},
            "task_label": "SCAN_nscf line"
            # "ehull":{"$lt":0.4},
            # "nkinds":2,
            # "magstate": "NM",
            # "uid": "Ti3ZrS8-CrW3S8-NM"
        }
        show = {
            "formula":1,
            "class":1,
            "ehull":1,
            "spacegroup":1,
            "gap_hse_nosoc":1
        }

        good_pg = []
        bad_pg = []
        host_candidates = list(col.find(filter))
        for e in host_candidates[:]:
            print(e["task_id"])
            print(e["task_id"], e["c2db_info"]["uid"])
            bs = SCAN2dMat.get_band_structure(e["task_id"])
            data = {}
            if bs.get_band_gap()["energy"] != 0:
                data = get_band_edges_characters(bs)
                data.update({"tk_id": e["task_id"]})


            st = Structure.from_dict(e["input"]["structure"])
            species, site_syms, wyckoffs, pg, site_idx, spg, spg_number = get_unique_sites_from_wy(st, symprec=1e-2).values()
            good_ir = get_good_ir_sites(st, symprec=1e-2)

            print(e["c2db_info"]["uid"])
            if good_ir["site_syms"]:
                d = {
                    "species": species,
                    "reduced_species": good_ir["species"],
                    "site_sym": site_syms,
                    "reduced_site_sym": good_ir["site_syms"],
                    "wyckoffs": wyckoffs,
                    "reduced_wy": good_ir["wys"],
                    "site_idx": site_idx,
                    "reduced_site_idx": good_ir["site_idx"],
                    "pmg_point_gp": pg,
                    "formula": e["c2db_info"]["formula"],
                    "spacegroup": e["c2db_info"]["spacegroup"],
                    "gap_hse": e["c2db_info"]["gap_hse"],
                    "gap_hse_nosoc": e["c2db_info"]["gap_hse_nosoc"],
                    "soc_intensity": e["c2db_info"]["gap_hse_nosoc"] - e["c2db_info"]["gap_hse"],
                    "ehull": e["c2db_info"]["ehull"],
                    "uid":e["c2db_info"]["uid"],
                }
                d.update(data)
                good_pg.append(d)
            else:
                d = {
                    "species": species,
                    "site_sym": site_syms,
                    "wyckoffs": wyckoffs,
                    "site_idx": site_idx,
                    "pmg_point_gp": pg,
                    "formula": e["c2db_info"]["formula"],
                    "spacegroup": e["c2db_info"]["spacegroup"],
                    "gap_hse": e["c2db_info"]["gap_hse"],
                    "gap_hse_nosoc": e["c2db_info"]["gap_hse_nosoc"],
                    "soc_intensity": e["c2db_info"]["gap_hse_nosoc"] - e["c2db_info"]["gap_hse"],
                    "ehull": e["c2db_info"]["ehull"],
                    "uid":e["c2db_info"]["uid"],
                }
                d.update(data)
                bad_pg.append(d)
        # print(bad_pg)
        good = pd.DataFrame(good_pg)
        bad = pd.DataFrame(bad_pg)
        #
        # # print(df)
        for name, df in zip(["good_scan_pg", "bad_scan_pg"], [good, bad]):
            df = df.round(3).sort_values(["soc_intensity", "ehull"], ascending=True)
            # df.set_index("uid", inplace=True)
        # df.to_json("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/gap_gt1-binary-NM-full.json",
        #            orient="index", indent=4, index=True)
            IOTools(pandas_df=df).to_excel(name)

    @classmethod
    def plot(cls):
        import pandas as pd
        from matplotlib import pyplot as plt

        DF = pd.read_excel("bad_pg.xlsx")
        def pmg_point_gp():
            fig = plt.figure(dpi=280)
            df = DF["pmg_point_gp"].value_counts()
            df.plot(kind="barh")
            # df["pmg_point_gp"].apply(pd.value_counts).plot(kind='bar', subplots=False)
            plt.show()

        def gap_hse_nosoc():
            fig = plt.figure(dpi=280)
            data = DF["gap_hse_nosoc"]
            max = C2DBStats.sheet()[2]
            min = C2DBStats.sheet()[1]

            data.plot.hist(bins=np.arange(min, max + 0.2, 0.2))
            plt.xlim([0, 10])
            plt.ylim([0, 70])
            plt.xticks(np.linspace(0,10,11))
            plt.grid()
            plt.show()

        gap_hse_nosoc()


class C2DBStats:
    @classmethod
    def sheet(cls):
        db = C2DB
        filter = {
            "magstate": "NM",
            "gap_hse_nosoc": {"$gte": 0}
        }

        es = db.collection.aggregate(
            [
                {"$match": filter},
                {"$group": {
                    "_id": {"gap": "$gap_hse_nosoc"},
                    "count": {"$sum":1},
                    "uid": {"$push": "$uid"}
                }},
                {"$project":
                    {
                        "_id": 0,
                        "gap_hse_nosoc": "$_id.gap",
                    }
                },
                {"$sort": {"gap":1}}
            ]
        )

        df = pd.DataFrame(es)
        print(df)
        return df, min(df["gap_hse_nosoc"]), max(df["gap_hse_nosoc"])


    @classmethod
    def plot(cls):
        from matplotlib import pyplot as plt
        data = cls.sheet()[0]["gap_hse_nosoc"]
        fig = plt.figure(dpi=280)
        cls.sheet()[0].plot.hist(bins=np.arange(min(data), max(data) + 0.2, 0.2))
        plt.xlim([0, 10])
        plt.ylim([0, 70])
        plt.xticks(np.linspace(0,10,11))
        plt.grid()
        plt.show()

class Scan2dMat:

    @classmethod
    def host_candidate(cls):
        db = C2DB_IR
        db.collection.find_one()
        col = db.collection
        filter = {
            "task_label": "hse line",
            # "c2db_info.gap_hse_nosoc": {"$gte": 1.5}
        }

        d = []
        for e in list(col.find(filter))[:]:
            print(e["task_id"], e["c2db_uid"])
            bs = db.get_band_structure(e["task_id"])
            if bs.get_band_gap()["energy"] == 0:
                continue
            data = get_band_edges_characters(bs)
            data.update({"tk_id": e["task_id"], "c2db_uid": e["c2db_uid"]})
            d.append(data)

        df = pd.DataFrame(d)
        print(df)
        df.to_excel("test.xlsx", index=False)

    @classmethod
    def update_entry_with_c2db(cls):
        for e in C2DB.collection.find({"gap_hse_nosoc": {"$gte": 0}, "magstate": "NM", "nkinds":2}):
            uid = e["uid"]
            print(uid)
            try:
                SCAN2dMat.collection.update_one(
                    {"task_label": "SCAN_nscf line", "c2db_info.uid": uid},
                    {
                        "$set": {
                            "c2db_info.formula": e["formula"],
                            "c2db_info.spacegroup": e["spacegroup"],
                            "c2db_info.gap_hse": e["gap_hse"],
                            "c2db_info.gap_hse_nosoc": e["gap_hse_nosoc"],
                            "c2db_info.soc_intensity": e["gap_hse_nosoc"] - e["gap_hse"],
                            "c2db_info.ehull": e["ehull"],
                            "c2db_info.uid":e["uid"],
                        }
                     }
                )
            except Exception as er:
                print(er)


class Scan2dDefect:
    @classmethod
    def run_time_stats(cls):
        col = SCAN2dDefect.collection
        data = []
        for e in list(col.find({"task_label": "SCAN_nscf uniform"}))[:]:
            print(e["formula_pretty"])
            run_time = e["run_stats"]["overall"]["Total CPU time used (sec)"]
            data.append({"task_id": e["task_id"], "run_time": round(run_time/60, 2)})

        df = pd.DataFrame(data)
        df["run_time"].mean()
        print(df["run_time"].mean(), df["run_time"].min(), df["run_time"].max())
        return df

    @classmethod
    def plot(cls):
        from matplotlib import pyplot as plt
        data = cls.run_time_stats()["run_time"]
        fig = plt.figure(dpi=280)
        data.plot.hist(bins=np.arange(min(data), max(data) + 0.5, 0.5))
        # plt.xlim([0, 15])
        # plt.ylim([0, 70])
        plt.xticks(np.linspace(0,15,21))
        plt.grid()
        plt.show()


def main():
    DefectStats.host_candidate_stats()




if __name__ == '__main__':
    main()


