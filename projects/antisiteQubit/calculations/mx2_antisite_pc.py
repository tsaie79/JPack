#%% mx2_antisite_pc
"""
HSE after relaxation
"""
from atomate.vasp.database import VaspCalcDb
from pymatgen import Structure
import pandas as pd
from matplotlib import pyplot as plt

col = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_pc/db.json").collection
data = {}
data1 = {}
for i in col.find({
    "chemsys": {"$in":["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"]},
    "formula_anonymous": "AB2",
    "calcs_reversed.run_type": "HSE06",
    "nsites": 3,
    "task_type":"JHSEStaticFW",
    "input.parameters.AEXX": 0.25,
    "input.structure.lattice.c": {"$nin":[20]}
}):
    # st = Structure.from_dict(i["calcs_reversed"][0]["output"]["structure"])
    st = Structure.from_dict(i["input"]["structure"])

# st.to("poscar", "/Users/jeng-yuantsai/Research/project/qubit/calculations/"
    #                 "mx2_antisite_pc/structures/{}_{}.vasp".format(st.formula, i["task_id"]))

    # st.to("poscar", path+"{}.vasp".format(i["formula_pretty"]))
    data[i["chemsys"]] = {"a": round(st.lattice.a,3), "b": round(st.lattice.b,3),
                          "c": round(st.lattice.c,3), "gamma": round(st.lattice.gamma,3),
                          "space_gp": i["output"]["spacegroup"]["symbol"]
                          }
    print(data)
    if "Te" in i["chemsys"]:
        data[i["chemsys"]].update({
            "x-x": round(st.get_distance(0, 1),3),
            "m-x": round(st.get_distance(2, 0),3),
            "x-x-m": round(st.get_angle(1,0,2),3),
            "task_id":i["task_id"]})
    else:
        data[i["chemsys"]].update({
            "x-x": round(st.get_distance(1, 2),3),
            "m-x": round(st.get_distance(0, 1),3),
            "x-x-m": round(st.get_angle(2,1,0),3),
            "task_id":i["task_id"]
        })

# data2 = pd.DataFrame(data)
data = pd.DataFrame(data).reindex(columns=["S-W", "Mo-S", "Se-W", "Mo-Se", "Te-W", "Mo-Te"])
data.to_clipboard()
# host = pd.concat([data2, data1], axis=0).reindex(columns=["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"])
# mo_angles = host.iloc[-2, 3:]
# w_angles = host.iloc[-2, :3]
# fig, ax = plt.subplots()
# ax.plot(["S", "Se", "Te"], mo_angles, label="Mo", marker="o")
# ax.plot(["S", "Se", "Te"], w_angles, label="W", marker="x")

# host.to_clipboard()
# print(host)

#%% mx2_antisite_pc
"""
SCAN after relax
"""
from atomate.vasp.database import VaspCalcDb
from pymatgen import Structure
import pandas as pd
from matplotlib import pyplot as plt

col = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc/"
                              "db.json").collection
data = {}
data1 = {}
for i in col.find({
    "chemsys": {"$in":["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"]},
    "task_label": "SCAN_relax"
}):
    # st = Structure.from_dict(i["calcs_reversed"][0]["output"]["structure"])
    st = Structure.from_dict(i["input"]["structure"])


    data[i["chemsys"]] = {"a": round(st.lattice.a,3), "b": round(st.lattice.b,3),
                          "c": round(st.lattice.c,3), "gamma": round(st.lattice.gamma,3),
                          "space_gp": i["output"]["spacegroup"]["symbol"]
                          }
    print(data)
    if "Te" in i["chemsys"]:
        data[i["chemsys"]].update({
            "x-x": round(st.get_distance(0, 1),3),
            "m-x": round(st.get_distance(2, 0),3),
            "x-x-m": round(st.get_angle(1,0,2),3),
            "task_id":i["task_id"]})
    else:
        data[i["chemsys"]].update({
            "x-x": round(st.get_distance(1, 2),3),
            "m-x": round(st.get_distance(0, 1),3),
            "x-x-m": round(st.get_angle(2,1,0),3),
            "task_id":i["task_id"]
        })

# data2 = pd.DataFrame(data)
data = pd.DataFrame(data).reindex(columns=["S-W", "Mo-S", "Se-W", "Mo-Se", "Te-W", "Mo-Te"])
data.to_clipboard()
# host = pd.concat([data2, data1], axis=0).reindex(columns=["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"])
# mo_angles = host.iloc[-2, 3:]
# w_angles = host.iloc[-2, :3]
# fig, ax = plt.subplots()
# ax.plot(["S", "Se", "Te"], mo_angles, label="Mo", marker="o")
# ax.plot(["S", "Se", "Te"], w_angles, label="W", marker="x")

# host.to_clipboard()
# print(host)
#%% mx2_antisite_pc
"""
HSE before relax
"""
from atomate.vasp.database import VaspCalcDb
from pymatgen import Structure
import pandas as pd
from matplotlib import pyplot as plt

col = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_pc/db.json").collection
data = {}
data1 = {}
for i in col.find({
    "task_id": {"$in": [3091, 3083, 3093, 3097, 3094, 3102]}
}):
    # st = Structure.from_dict(i["calcs_reversed"][0]["output"]["structure"])
    st = Structure.from_dict(i["input"]["structure"])

    # st.to("poscar", "/Users/jeng-yuantsai/Research/project/qubit/calculations/"
    #                 "mx2_antisite_pc/structures/{}_{}.vasp".format(st.formula, i["task_id"]))

    # st.to("poscar", path+"{}.vasp".format(i["formula_pretty"]))
    data[i["chemsys"]] = {"a": round(st.lattice.a,3), "b": round(st.lattice.b,3),
                          "c": round(st.lattice.c,3), "gamma": round(st.lattice.gamma,3),
                          "space_gp": i["output"]["spacegroup"]["symbol"]
                          }
    print(data)
    if "Te" in i["chemsys"]:
        data[i["chemsys"]].update({
            "x-x": round(st.get_distance(0, 1),3),
            "m-x": round(st.get_distance(2, 0),3),
            "x-x-m": round(st.get_angle(1,0,2),3),
            "task_id":i["task_id"]})
    else:
        data[i["chemsys"]].update({
            "x-x": round(st.get_distance(1, 2),3),
            "m-x": round(st.get_distance(0, 1),3),
            "x-x-m": round(st.get_angle(2,1,0),3),
            "task_id":i["task_id"]
        })

# data2 = pd.DataFrame(data)
data = pd.DataFrame(data).reindex(columns=["S-W", "Mo-S", "Se-W", "Mo-Se", "Te-W", "Mo-Te"])
data.to_clipboard()
# host = pd.concat([data2, data1], axis=0).reindex(columns=["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"])
# mo_angles = host.iloc[-2, 3:]
# w_angles = host.iloc[-2, :3]
# fig, ax = plt.subplots()
# ax.plot(["S", "Se", "Te"], mo_angles, label="Mo", marker="o")
# ax.plot(["S", "Se", "Te"], w_angles, label="W", marker="x")

# host.to_clipboard()
# print(host)

#%% mx2_antisite_basic_aexx0.25_sigma_test
"""
HSE defect relaxed structures
"""
from atomate.vasp.database import VaspCalcDb
from pymatgen import Structure
import pandas as pd
import os
from matplotlib import pyplot as plt

# p = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_sigma_test"
p = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/search_triplet_from_defect_db"

col = VaspCalcDb.from_db_file(os.path.join(p, "db.json")).collection
data = {}
data1 = {}
for i in col.find({
    "task_id": {"$in": [4228,
                        4266,
                        4267,
                        4270,
                        4269,
                        4271]}
}):
    # st = Structure.from_dict(i["calcs_reversed"][0]["output"]["structure"])
    st = Structure.from_dict(i["input"]["structure"])


    data[i["chemsys"]] = {"a": round(st.lattice.a,3), "b": round(st.lattice.b,3),
                          "c": round(st.lattice.c,3), "gamma": round(st.lattice.gamma,3),
                          "space_gp": i["output"]["spacegroup"]["symbol"]
                          }
    print(data)
    if "Te" in i["chemsys"]:
        data[i["chemsys"]].update({
            "x-x-m": round(st.get_angle(24,74,55),3),
            "task_id":i["task_id"]})
    else:
        data[i["chemsys"]].update({
            "x-x-m": round(st.get_angle(50,25,0),3),
            "task_id":i["task_id"]
        })

# data2 = pd.DataFrame(data)
data = pd.DataFrame(data).reindex(columns=["S-W", "Mo-S", "Se-W", "Mo-Se", "Te-W", "Mo-Te"])
data.to_clipboard()
# host = pd.concat([data2, data1], axis=0).reindex(columns=["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"])
# mo_angles = host.iloc[-2, 3:]
# w_angles = host.iloc[-2, :3]
# fig, ax = plt.subplots()
# ax.plot(["S", "Se", "Te"], mo_angles, label="Mo", marker="o")
# ax.plot(["S", "Se", "Te"], w_angles, label="W", marker="x")

# host.to_clipboard()
# print(host)

#%% search_triplet_from_defect_db
"""
HSE defect relaxed structures
"""
from atomate.vasp.database import VaspCalcDb
from pymatgen import Structure
import pandas as pd
import os
from matplotlib import pyplot as plt

# p = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_sigma_test"
p = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/search_triplet_from_defect_db"

col = VaspCalcDb.from_db_file(os.path.join(p, "db.json")).collection
data = {}
data1 = {}
for i in col.find({
    "task_id": {"$in": [347,
                        435,
                        303,
                        322,
                        256,
                        446]}
}):
    # st = Structure.from_dict(i["calcs_reversed"][0]["output"]["structure"])
    st = Structure.from_dict(i["input"]["structure"])


    data[i["chemsys"]] = {"a": round(st.lattice.a,3), "b": round(st.lattice.b,3),
                          "c": round(st.lattice.c,3), "gamma": round(st.lattice.gamma,3),
                          "space_gp": i["output"]["spacegroup"]["symbol"]
                          }
    print(data)
    if "Te" in i["chemsys"]:
        data[i["chemsys"]].update({
            "x-x-m": round(st.get_angle(24,74,55),3),
            "task_id":i["task_id"]})
    else:
        data[i["chemsys"]].update({
            "x-x-m": round(st.get_angle(50,25,0),3),
            "task_id":i["task_id"]
        })

data = pd.DataFrame(data).reindex(columns=["S-W", "Mo-S", "Se-W", "Mo-Se", "Te-W", "Mo-Te"])
data.to_clipboard()
