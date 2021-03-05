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
for row in ti[0]:
    icsd_id = row.split("_")[1]
    print(icsd_id)
    with MPRester() as mp:
        try:
            e = mp.query({"icsd_ids": int(icsd_id)}, ["cifs.conventional_standard", "task_id"])
            st = Structure.from_str(e[0]["cifs.conventional_standard"], fmt="cif").\
                to("poscar", "ti/icsd-{}_{}.vasp".format(icsd_id, e[0]["task_id"]))
        except Exception as err:
            er.append((icsd_id, row))
            print("$$"*10)
            print(err)
dumpfn(dict(er), "ti/err_icsd.json")
#%%
import pandas as pd
import os, glob
from pymatgen import MPRester, Structure
from monty.serialization import dumpfn

dir_path = "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/"
os.chdir(dir_path)
# tci = pd.read_csv("icsd_num_layer_materials_TCI.dat", header=None)
# ti = pd.read_csv("icsd_num_layer_materials_TI.dat", header=None)

gap = []
for row in glob.glob("ti_slab_hex/*"):
    row = row.split("/")[-1]
    icsd_id = row.split("_")[0]
    icsd_id = icsd_id.split("-")[1]
    print(icsd_id)
    with MPRester() as mp:
        try:
            e = mp.query({"icsd_ids": int(icsd_id)}, ["band_gap"])
            gap.append((row, e[0]["band_gap"]))
            print(e[0]["band_gap"])
        except Exception as err:
            print("$$"*10)
            print(err)
dumpfn(dict(gap), "ti_slab_hex/pbe_gap.json", indent=4)
#%% check angles
from pymatgen import Structure
import glob
from mpinterfaces.utils import *


dir_path = "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/"
os.chdir(dir_path)

for i in glob.glob("ti/*.vasp"):
    st = Structure.from_file(i)
    angles = [round(i, 0)for i in st.lattice.angles]
    if st.lattice.is_hexagonal():
        print(i)
        print(angles)
        print(st.num_sites)
        slab = get_ase_slab(st, hkl=[0,0,1], min_thick=40, min_vac=15)
        os.makedirs("ti_slab", exist_ok=True)
        slab.sort()
        slab.to("poscar", "ti_slab/{}".format(i.split("/")[-1]))

#%%
