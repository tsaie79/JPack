#%% PC
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
    st = Structure.from_dict(i["calcs_reversed"][0]["output"]["structure"])
    st.to("poscar", "/Users/jeng-yuantsai/Research/project/qubit/calculations/"
                    "mx2_antisite_pc/structures/{}_{}.vasp".format(st.formula, i["task_id"]))

    # st.to("poscar", path+"{}.vasp".format(i["formula_pretty"]))
    data[i["chemsys"]] = {"a": round(st.lattice.a,3), "b": round(st.lattice.b,3),
                          "c": round(st.lattice.c,3), "gamma": round(st.lattice.gamma,3),
                          "space_gp": i["output"]["spacegroup"]["symbol"]
                          }
    if "Te" in i["chemsys"]:
        data1[i["chemsys"]] = {
            "x-x": round(st.get_distance(0, 1),3),
            "m-x": round(st.get_distance(2, 0),3),
            "x-x-m": round(st.get_angle(1,0,2),3),
            "task_id":i["task_id"]}
    else:
        data1[i["chemsys"]] = {
            "x-x": round(st.get_distance(1, 2),3),
            "m-x": round(st.get_distance(0, 1),3),
            "x-x-m": round(st.get_angle(2,1,0),3),
            "task_id":i["task_id"]
        }

data2 = pd.DataFrame(data)
data1 = pd.DataFrame(data1)
host = pd.concat([data2, data1], axis=0).reindex(columns=["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"])
mo_angles = host.iloc[-2, 3:]
w_angles = host.iloc[-2, :3]
fig, ax = plt.subplots()
ax.plot(["S", "Se", "Te"], mo_angles, label="Mo", marker="o")
ax.plot(["S", "Se", "Te"], w_angles, label="W", marker="x")

# host.to_clipboard()
# print(host)

