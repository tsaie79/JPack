class Pc:
    @classmethod
    def hist_nsites(cls):
        import pandas as pd
        from matplotlib import pyplot as plt
        from qubitPack.tool_box import get_db

        db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan")
        filter = {
            "gap_hse_nosoc":{"$gt":0},
            "nkinds":2,
            "magstate": "NM"
        }

        es = db.collection.aggregate(
            [
                {"$match": filter},
                {"$group": {
                    "_id": {"nsites": "$nsites"},
                    "count": {"$sum":1},
                    "uid": {"$push": "$uid"}
                }},
                {"$project":
                    {
                    "_id": 0,
                    "nsites": "$_id.nsites",
                    "count": "$count",
                    }
                },
                {"$sort": {"nsites":1}}
            ]
        )

        df = pd.DataFrame(es)
        print(df)

        df.plot(x="nsites", y="count", kind="bar", rot=0, figsize=(10,8), title="Hist nsites under filter")
        plt.ylabel("counts")

        plt.show()





if __name__ == '__main__':
    Pc.hist_nsites()