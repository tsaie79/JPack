#%% IMPORT
from atomate.vasp.database import VaspCalcDb

db_dict = {
    "host": "localhost",
    "port": 1234,
    "database": "single_photon_emitter",
    "collection": "standard_defect",
    "admin_user": "Jeng",
    "admin_password": "qimin",
    "readonly_user": "Jeng_ro",
    "readonly_password": "qimin",
    "aliases": {}
}
def db(collection):
    return VaspCalcDb(host="localhost", port=1234, database="single_photon_emitter",
                collection=collection, user="Jeng_ro", password="qimin", authsource="single_photon_emitter")


#%% WS2 with C3V
from pymatgen import Structure

col = db("cdft").collection

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
from pymatgen import Structure

col = db("soc_standard_defect").collection
a = col.find_one({"task_id": 70})

col = db("soc_cdft").collection

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

dx = np.linspace(0,2,11)
get_interpolate_sts(a_st, c_st, disp_range=dx, output_dir='/Users/jeng-yuantsai/Research/project/single_photon/calculations/soc_cdft/Ws_c3v/interpolate_st/excited_state')


#%%

import pandas as pd
es = db("cdft").collection.aggregate([
    {"$match": {"cdft_info": "excited_state_C3V_in_paper"}},
    {"$project": {"energy": "$output.energy", "idx": "$interpolate_st_idx", "_id":0}}
])

df = pd.DataFrame(es)
df = df.sort_values("idx")
df.to_clipboard()

#%%
from qubitPack.qc_searching.analysis.main import get_defect_state

cdft = db("standard_defect")

tot, proj, d_df = get_defect_state(
    cdft,
    {"task_id": 22},
    1,-3,
    None,
    True,
    "dist",
    None,
    0.1
)