#%% print ZPL
from pymatgen import Structure
from qubitPack.tool_box import get_db
import pandas as pd
from qubitPack.tool_box import get_db
from pymatgen import Structure

col = get_db("single_photon_emitter", "soc_cdft").collection

chemsys = ["S-W", "Se-W", "Te-W"]
chem, E_ab, E_bc, E_cd, E_da, ZPL =[], [], [], [], [], []
for chemsys in chemsys:
    b = col.find_one({"task_label":{"$regex": "B-HSE_scf_soc"}, "chemsys":chemsys})
    c = col.find_one({"task_label":{"$regex": "C-HSE_relax_soc"}, "chemsys":chemsys})
    d = col.find_one({"task_label":{"$regex": "D-HSE_scf_soc"}, "chemsys":chemsys})

    E_ab.append(b["output"]["energy"] - b["source"]["total_energy"])
    E_bc.append(b["output"]["energy"] - c["output"]["energy"])
    E_cd.append(c["output"]["energy"] - d["output"]["energy"])
    E_da.append(d["output"]["energy"] - d["source"]["total_energy"])
    ZPL.append(c["output"]["energy"] - b["source"]["total_energy"])
    chem.append(chemsys)

df = pd.DataFrame({"chem": chem,"E_ab":E_ab, "E_bc":E_bc, "E_cd":E_cd, "E_da": E_da, "ZPL":ZPL})
df = df.round(4)
df.to_clipboard(excel=True)

#%% print raw data task ids
from pymatgen import Structure
from qubitPack.tool_box import get_db
import pandas as pd
from qubitPack.tool_box import get_db
from pymatgen import Structure

col = get_db("single_photon_emitter", "soc_cdft").collection

es = col.aggregate(
    [
        {"$match": {"dx": {"$exists": False}, "$or": [{"task_label": {"$regex": "HSE_relax_soc"}}, {"task_label": {"$regex": "HSE_scf_soc"}}]}},
        {"$group": {"_id": {"chemsys": "$chemsys"}, "tk_id": {"$push":"$task_id"}, "task_label": {"$push":"$task_label"},
                    "e": {"$push":"$output.energy"}, "A": {"$push": "$source.total_energy"}}},
        {"$project": {"_id":0, "chemsys": "$_id.chemsys", "task_id": "$tk_id", "task_labe":"$task_label", "e": "$e", "A": "$A"}}
    ]
)

df = pd.DataFrame(es)
df.to_clipboard(excel=True)

#%% SOC cdft AEPS

from pymatgen import Structure
from qubitPack.tool_box import get_db
import pandas as pd
from qubitPack.tool_box import get_db
from pymatgen.core.units import eV_to_Ha
from pymatgen import Structure
import matplotlib.pyplot as plt
import numpy as np

col = get_db("single_photon_emitter", "soc_cdft").collection

es = col.aggregate(
    [
        {"$match": {"dx": {"$exists":1}, "chemsys":"Te-W"}},
        {"$group": {"_id": {"task_label": "$task_label"}, "tk_id": {"$push":"$task_id"},
                    "dx": {"$push": "$dx"}, "E": {"$push": "$output.energy"}, "info": {"$push": "$info"}}},
        {"$project": {"_id":0, "tk_label": "$_id.task_label", "task_id": "$tk_id", "dx": "$dx", "energy": "$E", "info":"$info"}}
    ]
)
es = list(es)
info_g = es[1]["info"][1]
info_e = es[1]["info"][0]
info_e = pd.DataFrame(info_e)
info_e.to_clipboard()

Q_g = np.array(es[1]["dx"])*info_g["Delta_Q"]
E_g = np.array(es[1]["energy"])*eV_to_Ha
qe_g = dict(zip(es[1]["dx"], E_g))
Eg_relax = qe_g[1.0]-qe_g[0.0]
fit_normalized = np.polyfit(Q_g, E_g, 2)
f_g = np.poly1d(fit_normalized)



Q_e = np.array(es[0]["dx"])*info_e["Delta_Q"]
E_e = np.array(es[0]["energy"])*eV_to_Ha
print(es[0]["dx"])
qe_e = dict(zip(es[0]["dx"], E_e))
Ee_relax = qe_e[0.0]-qe_e[1.0]
fit_normalized = np.polyfit(Q_e, E_e, 2)
f_e = np.poly1d(fit_normalized)

Q_e, Q_g, E_e, E_g = map(pd.DataFrame, zip([Q_e, Q_g, E_e, E_g]))
Q = pd.concat([Q_e, Q_g], axis=1)
E = pd.concat([E_e, E_g], axis=1)

data = pd.concat([Q,E], axis=0, ignore_index=True)
data = data.transpose()
df = data.sort_values(0)


x, y = df[0], df[1]
plt.scatter(x, y, marker='o', c='red')
plt.plot(Q_g,f_g(Q_g), c='black')
plt.plot(Q_e,f_e(Q_e), c='black')
plt.show()

#%%
from qubitPack.tool_box import get_db
import pandas as pd
es = get_db("single_photon_emitter", "soc_cdft", port=1234).collection.aggregate([
    {"$match": {"dx": {"$exists":True}, "chemsys":"Te-W", "task_label": {"$regex": "B-HSE_scf_soc"}}},
    {"$project": {"energy": "$output.energy", "dx":"$dx", "task_id":1, "info":1, "_id":0}}
])

df = pd.DataFrame(es)
df = df.sort_values("dx")
print(df["info"][0])
df.to_clipboard()

#%%
import numpy as np
from qubitPack.tool_box import get_db
from atomate.vasp.jpowerups import *
from qubitPack.tool_box import get_interpolate_sts


tk_id = 268
db_a = get_db("single_photon_emitter", "soc_standard_defect", port=1234)
entry_a = db_a.collection.find_one({"task_id":tk_id})
a_path = entry_a["dir_name"].split("/")[-1]

nonsoc_a_entry = entry_a["nonsoc_from"].split("/")
db_nonsoc_a = get_db(nonsoc_a_entry[0], nonsoc_a_entry[1], port=1234)
nonsoc_a_path = db_nonsoc_a.collection.find_one({"task_id":int(nonsoc_a_entry[-1])})["dir_name"].split("/")[-1]


soc_cdft = get_db("single_photon_emitter", "soc_cdft", port=1234)
b_st = Structure.from_dict(soc_cdft.collection.find_one(
    {"chemsys":entry_a["chemsys"], "task_label": {"$regex": "B-HSE_scf_soc"}})["input"]["structure"])
c_st = Structure.from_dict(soc_cdft.collection.find_one(
    {"chemsys":entry_a["chemsys"], "task_label": {"$regex": "C-HSE_relax_soc"}})["output"]["structure"])

dx = np.linspace(0, 10,11)
sts, info = get_interpolate_sts(
    b_st,
    c_st,
    disp_range=[1],
    output_dir=None
)
for frac, st in sts.items():
    st.to("poscar", "/Users/jeng-yuantsai/Research/project/single_photon/calculations/"
                    "soc_cdft/WX2_Ch/S-W/{}.vasp".format("_".join(str(frac).split("."))))

#%%
from qubitPack.qc_searching.analysis.main import get_defect_state
from qubitPack.tool_box import get_db

db = get_db("single_photon_emitter", "standard_defect")
tot, proj, d_df = get_defect_state(
    db,
    {"task_id": 22},
    1,-2,
    None,
    True,
    "dist",
    None,
    0.1
)