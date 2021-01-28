#%% WS2 with C3V
from qubitPack.tool_box import get_db
from pymatgen import Structure

col = get_db("single_photon_emitter","cdft").collection

b = col.find_one({"task_id": 23})
c = col.find_one({"task_id": 24})
d = col.find_one({"task_id": 25})


b_st = Structure.from_dict(b["input"]["structure"])
c_st = Structure.from_dict(c["output"]["structure"])

E_ab = b["output"]["energy"] - b["source"]["total_energy"]
E_bc = c["output"]["energy"] - b["output"]["energy"]
E_cd = d["output"]["energy"] - c["output"]["energy"]
E_da = d["output"]["energy"] - d["source"]["total_energy"]

#%% WS2_C3V_IN_PAPER
from pymatgen import Structure

col = db("mx2_antisite_basic_aexx0.25_cdft").collection

b = col.find_one({"task_id": 4179})
c = col.find_one({"task_id": 4180})
d = col.find_one({"task_id": 4181})


b_st = Structure.from_dict(b["input"]["structure"])
c_st = Structure.from_dict(c["output"]["structure"])

E_ab = b["output"]["energy"] - b["source"]["total_energy"]
E_bc = c["output"]["energy"] - b["output"]["energy"]
E_cd = d["output"]["energy"] - c["output"]["energy"]
E_da = d["output"]["energy"] - d["source"]["total_energy"]

#%%
from qubitPack.tool_box import get_db
from pymatgen import Structure

col = get_db("single_photon_emitter","soc_standard_defect").collection
a = col.find_one({"task_id": 70})

col = get_db("single_photon_emitter","soc_cdft").collection

b = col.find_one({"task_id": 73})
c = col.find_one({"task_id": 76})
d = col.find_one({"task_id": 75})

a_st = Structure.from_dict(b["input"]["structure"])
b_st = Structure.from_dict(b["input"]["structure"])
c_st = Structure.from_dict(c["output"]["structure"])

E_ab = b["output"]["energy"] - b["source"]["total_energy"]
E_bc = c["output"]["energy"] - b["output"]["energy"]
E_cd = d["output"]["energy"] - c["output"]["energy"]
E_da = d["output"]["energy"] - d["source"]["total_energy"]

#%% interpolate bewteen excited eq and ground eq
from qubitPack.tool_box import get_interpolate_sts
from pymatgen import Structure
import numpy as np

dx = np.linspace(-1,-21,5)
get_interpolate_sts(b_st,
                    c_st,
                    disp_range=dx,
                    output_dir='/Users/jeng-yuantsai/Research/project/single_photon/calculations/'
                               'cdft/WS2_c3v/interpolate_st/ground_state/extra/')


#%%
from qubitPack.tool_box import get_db
import pandas as pd
es = get_db("single_photon_emitter", "cdft", port=1234).collection.aggregate([
    {"$match": {"cdft_info": "ground_state_C3V", "poscar_idx": {"$exists":False}}},
    {"$project": {"energy": "$output.energy", "poscar_idx":"$interpolate_st_idx", "task_id":1, "_id":0}}
])

df = pd.DataFrame(es)
df = df.sort_values("poscar_idx")
df.to_clipboard()

#%%
from qubitPack.qc_searching.analysis.main import get_defect_state
from qubitPack.tool_box import get_db
from pymatgen.electronic_structure.plotter import DosPlotter

cdft = get_db("single_photon_emitter", "soc_standard_defect")

dos = cdft.get_dos(70)
ds_plotter = DosPlotter()
ds_plotter.add_dos("t", dos)
plt = ds_plotter.get_plot(xlim=[-5,5])
plt.show()
# tot, proj, d_df = get_defect_state(
#     cdft,
#     {"task_id": 72},
#     1,-3,
#     None,
#     True,
#     "dist",
#     None,
#     0.1

#%%
from qubitPack.tool_box import get_db
soc = get_db("single_photon_emitter", "soc_standard_defect")
nosoc = get_db("single_photon_emitter", "standard_defect")

chemsys = "S-W"

en_soc = soc.collection.find_one({"chemsys":chemsys, "task_label":"HSE_soc",
                                  "output.spacegroup.symbol":"P3m1"})["output"]["energy"]
en_nosoc = nosoc.collection.find_one({"chemsys":chemsys, "task_label":"HSE_scf"})["output"]["energy"]
