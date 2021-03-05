from qubitPack.tool_box import get_db
from pymatgen import Structure


db = get_db("2dMaterials_from_icsd_v1", "VaspTasks_v1", password="qiminyan", user="readUser")
e = db.collection.find_one({"icsd_id": "166578"})
st = Structure.from_dict(e["output"]["structure"])
st.to("poscar","/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/t.vasp")

#%%
from pymatgen import MPRester, Structure

with MPRester() as mp:
    e = mp.query({"task_id": "mp-553942"}, ["cifs.conventional_standard"])
st = Structure.from_str(e[0]["cifs.conventional_standard"], fmt="cif").to(
    "poscar", "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/t.vasp")

#%%
import pandas as pd
import os
from pymatgen import MPRester, Structure
from monty.serialization import dumpfn

dir_path = "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/"
os.chdir(dir_path)
tci = pd.read_csv("icsd_num_layer_materials_TCI.dat", header=None)
ti = pd.read_csv("icsd_num_layer_materials_TI.dat", header=None)

er = []
collection = []
for row in ti[0]:
    icsd_id = row.split("_")[1]
    print(icsd_id)
    with MPRester() as mp:
        try:
            e = mp.query({"icsd_ids": int(icsd_id)}, ["cifs.conventional_standard", "task_id",
                                                      "band_gap", "magnetic_type", "total_magnetization"])
            collection.append((icsd_id, e[0]))
            print(collection)
        except Exception as err:
            er.append((icsd_id, row))
            print("$$"*10)
            print(err)
dumpfn(dict(er), "ti/err_icsd.json", indent=4)
dumpfn(dict(collection), "ti/orig.json", indent=4)

#%%
import pandas as pd
import os, glob
from pymatgen import MPRester, Structure
from monty.serialization import dumpfn

dir_path = "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/"
os.chdir(dir_path)

gap = []
for row in glob.glob("ti_slab_hex/*.vasp"):
    row = row.split("/")[-1]
    icsd_id = row.split("_")[0]
    icsd_id = icsd_id.split("-")[1]
    print(icsd_id)
    with MPRester() as mp:
        try:
            e = mp.query({"icsd_ids": int(icsd_id)}, ["band_gap", "magnetic_type", "total_magnetization"])
            gap.append((row, e[0]))
            print(e[0])
        except Exception as err:
            print("$$"*10)
            print(err)
dumpfn(dict(gap), "ti_slab_hex/pbe_gap.json", indent=4)
#%% check angles
from pymatgen import Structure
from glob import glob
from mpinterfaces.utils import *
from monty.serialization import loadfn


dir_path = "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/"
os.chdir(dir_path)
es = loadfn("ti/orig.json")
col = []
for icsd_id, v in es.items():
    print(icsd_id)
    st = Structure.from_str(v["cifs.conventional_standard"], fmt="cif")
    if st.lattice.is_hexagonal():
        print(st.num_sites)
        slab = get_ase_slab(st, hkl=[0,0,1], min_thick=40, min_vac=15)
        os.makedirs("ti_slab_hex", exist_ok=True)
        slab.sort()
        v["cifs.conventional_standard"] = slab.to("cif")
        col.append((icsd_id, v))
dumpfn(dict(col), "ti_slab_hex/ti_hex_slabs.json", indent=4)

#%%
