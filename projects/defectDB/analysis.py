class Pc:
    @classmethod
    def hist_nsites(cls):
        import pandas as pd
        from matplotlib import pyplot as plt
        from qubitPack.tool_box import get_db

        db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan")
        filter = {
            "magstate": "NM",
            "gap_hse_nosoc": {"$gt":0},
            "nkinds": 3
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
                    "gap": "$_id.gap",
                    "count": "$count",
                    }
                },
                {"$sort": {"gap":1}}
            ]
        )

        df = pd.DataFrame(es)
        print(df)

        df.plot(x="gap", y="count", kind="bar", rot=0, figsize=(10,8), title="Hist nsites under filter")
        plt.ylabel("counts")

        plt.show()


class Defect:
    @classmethod
    def find_defect_state(cls):
        from qubitPack.qc_searching.analysis.main import get_defect_state
        from qubitPack.tool_box import get_db

        defect_db = get_db("defect_db", "binary_defect")
        host_db = get_db("Scan2dMat", "calc_data")


        tk_id = 17
        # pc_from_id = defect_db.collection.find_one({"task_id": tk_id})["pc_from_id"]
        # c2db_uid = host_db.collection.find_one({"task_id": pc_from_id})["c2db_info"]["uid"]

        tot, proj, d_df = get_defect_state(
            defect_db,
            {"task_id": tk_id},
            1, -5,
            None,
            True,
            "proj",
            (host_db, 545, 0), #(get_db("antisiteQubit", "W_S_Ef"), 312, 0.),
            0.1,
            locpot_c2db=None #(c2db, c2db_uid, 0)
        )


if __name__ == '__main__':
    Defect.find_defect_state()
    # Pc.hist_nsites()