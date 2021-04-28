from qubitPack.tool_box import get_db
from pymatgen import Structure


db = get_db("2dMaterials_from_icsd_v1", "VaspTasks_v1", password="qiminyan", user="readUser")
e = db.collection.find_one({"icsd_id": "027394"})
st = Structure.from_dict(e["output"]["structure"])
st.to("poscar","/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/ti_slab_nonhex/t.vasp")

#%%
from pymatgen import MPRester, Structure

with MPRester() as mp:
    e = mp.query({"task_id": "mp-13963"}, ["cifs.conventional_standard"])
st = Structure.from_str(e[0]["cifs.conventional_standard"], fmt="cif").to(
    "poscar", "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/ti_slab_nonhex/.vasp")

#%% grab both bulk and monolayer
from pymatgen import MPRester, Structure
from qubitPack.tool_box import get_db
from pymatgen import Structure
from monty.serialization import loadfn
import os

os.chdir("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/ti_slab_nonhex")
es = loadfn("ti_nonhex_thinner_slabs.json")

db = get_db("2dMaterials_from_icsd_v1", "VaspTasks_v1", password="qiminyan", user="readUser")

spg_list = []
for i in es:
    print(es[i]["mp_spg"])
    spg_list.append(es[i]["mp_spg"])
spg_list = set(spg_list)

for spg in spg_list:
    for icsd in es:
        if es[icsd]["mp_spg"] == spg and es[icsd]["monolayer_scan_gap"] != 0:
            if "/" in spg:
                dir_name = "-".join(spg.split("/"))
            else:
                dir_name = spg
            os.makedirs(dir_name, exist_ok=True)
            e = db.collection.find_one({"icsd_id": icsd})
            Structure.from_dict(e["output"]["structure"]).to("poscar","{}/monolayer_icsd_{}.vasp".format(dir_name,icsd))
            Structure.from_str(es[icsd]["cifs.conventional_standard"], fmt="cif").to(
                "poscar", "{}/bulk_mpid_{}.vasp".format(dir_name, es[icsd]["task_id"]))

#%% create slab from bulk
from pymatgen.core.surface import generate_all_slabs, SlabGenerator
from pymatgen import Structure

bulk = Structure.from_file('/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/ti_slab_nonhex/Cmcm/bulk_mpid_mp-605.vasp')
slab_gen = SlabGenerator(bulk, [0,1,0], 10, 10)
slab = slab_gen.get_slab()
slab.to("poscar", '/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/ti_slab_nonhex/Cmcm/slab.vasp')
# slabs = generate_all_slabs(bulk, 1, 10, 10)
# for i in slabs:
#     print(i.miller_index, i.shift)


#%%
import pandas as pd
import os
from pymatgen import MPRester, Structure
from monty.serialization import dumpfn
from qubitPack.tool_box import get_db

dir_path = "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/"
os.chdir(dir_path)

db = get_db("2dMaterials_from_icsd_v1", "VaspTasks_v1", user="readUser", password="qiminyan")
col = db.collection

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
            mp_e = e[0]
            en = col.find_one({"icsd_id": icsd_id, "functional_type":"SCAN", "task_type":"static"})
            try:
                m_gap = en["output"]["bandgap"]
            except Exception:
                m_gap = None
            mp_e.update({"monolayer_scan_gap": m_gap, "formula_pretty": en["formula_pretty"]})
            collection.append((icsd_id, mp_e))
            print(collection)
        except Exception as err:
            er.append((icsd_id, row))
            print("$$"*10)
            print(err)
dumpfn(dict(er), "ti/err_icsd.json", indent=4)
dumpfn(dict(collection), "ti/orig.json", indent=4)


#%% create slab for hexagonal lattice
from pymatgen import Structure
from glob import glob
from mpinterfaces.utils import *
from monty.serialization import loadfn
from pymatgen.core.surface import SlabGenerator
from mpinterfaces.utils import ensure_vacuum

dir_path = "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/"
os.chdir(dir_path)
es = loadfn("ti/orig.json")
col = []
for icsd_id, v in es.items():
    st = Structure.from_str(v["cifs.conventional_standard"], fmt="cif")
    if not st.lattice.is_hexagonal():
        print(icsd_id)
        # slab = SlabGenerator(st, [0,0,1], 30, 15)
        # slab = slab.get_slabs()[-1]
        # print(icsd_id)
        # print(slab.is_symmetric(), slab.is_polar())
        # slab = ensure_vacuum(slab, 15)
        #
        # slab.sort()
        # v["cifs.conventional_standard"] = slab.to("cif")
        col.append((icsd_id, v))


os.makedirs("ti_slab_nonhex", exist_ok=True)
dumpfn(dict(col), "ti_slab_nonhex/ti_nonhex_thinner_slabs.json", indent=4)



#%%
from qubitPack.tool_box import get_db
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen import Structure
import os
from matplotlib import pyplot as plt

os.chdir("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/ti_slab_hex")
d = get_db("mag_insulator", "slab")
e = d.collection.find_one({"task_id": 249})

# locpot = e["calcs_reversed"][0]["output"]["locpot"]["2"]
# fig = plt.figure()
# ax = fig.add_subplot(1,1,1)
# ax.plot(locpot)
# plt.show()

# st = Structure.from_dict(e["output"]["structure"])
# st.to("poscar", "poscar.vasp")

bs = d.get_band_structure(249)

plotter = BSPlotter(bs)
plt = plotter.get_plot(vbm_cbm_marker=True)
plt.show()


#%%
import ase.build
from ase.visualize import view
from ase.geometry.dimensionality import analyze_dimensionality
from pymatgen import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.geometry.dimensionality import isolate_components

st = Structure.from_file("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/"
                         "ti_slab_nonhex/Cmcm/monolayer_icsd_653212.vasp")

#%%
from mpinterfaces.calibrate import CalibrateSlab
from mpinterfaces.interface import Interface
from mpinterfaces.transformations import *
from mpinterfaces.utils import *
from mpinterfaces.old_transformations import generate_all_configs

# os.chdir("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/tci_slab_nonhex/C2-m/")
os.chdir('/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/ti_slab_hex/graphite')

separation = 3  # in angstroms
nlayers_2d = 1
nlayers_substrate = 1

substrate_bulk = Structure.from_file("t.vasp")
st = Structure.from_file("t.vasp")
# substrate_bulk = get_struct_from_mp('Ag')
sa_sub = SpacegroupAnalyzer(substrate_bulk)
# substrate_bulk = sa_sub.get_conventional_standard_structure()
substrate_slab = Interface(substrate_bulk,
                           hkl=[0, 0, 1],
                           min_thick=10,
                           min_vac=25,
                           primitive=False, from_ase=True)
# substrate_slab = slab_from_file([0,0,1], 'POSCAR_substrate')
mat2d_slab = Structure(st.lattice, st.species, [np.multiply(site.frac_coords, [1,1,-1]) for site in st.sites])

mat2d_slab = Interface(mat2d_slab,
                           hkl=[0, 0, 1],
                           min_thick=10,
                           min_vac=25,
                           primitive=False, from_ase=True)

# get the in-plane lattice aligned slabs
# substrate_slab.to(fmt='poscar', filename='POSCAR_substrate_slab.vasp')
# mat2d_slab.to(fmt='poscar', filename='POSCAR_mat2d_slab.vasp')
# selective dynamics flag
# sd_flags = CalibrateSlab.set_sd_flags(
#     interface=substrate_slab,
#     n_layers=nlayers_substrate,
#     top=True, bottom=False)
# poscar = Poscar(substrate_slab, selective_dynamics=sd_flags)
# poscar.write_file(filename='POSCAR_substrate_slab.vasp')
# get aligned lattices
substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    substrate_slab,
    mat2d_slab,
    max_area=400,
    max_mismatch=0.05,
    max_angle_diff=1,
    r1r2_tol=0.01)
substrate_slab_aligned.to(fmt='poscar',
                          filename='POSCAR_substrate_aligned.vasp')
mat2d_slab_aligned.to(fmt='poscar',
                      filename='POSCAR_mat2d_aligned.vasp')
# merge substrate and mat2d in all possible ways
hetero_interfaces = generate_all_configs(mat2d_slab_aligned,
                                         substrate_slab_aligned,
                                         nlayers_2d, nlayers_substrate,
                                         separation)

# generate all poscars
for i, iface in enumerate(hetero_interfaces):
    sd_flags = CalibrateSlab.set_sd_flags(
        interface=iface,
        n_layers=nlayers_2d + nlayers_substrate,
        top=True, bottom=False)
    poscar = Poscar(iface, selective_dynamics=sd_flags)
    poscar.write_file(
        filename='POSCAR_final_{}.vasp'.format(i))

#%%
os.chdir("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/tci_slab_nonhex/C2-m/")
st = Structure.from_file("bulk_icsd_657711.vasp")

from pymatgen.core.surface import generate_all_slabs, SlabGenerator
from mpinterfaces.utils import center_slab

slab_gen = SlabGenerator(st, [1,0,0], 20, 10)
slabs = slab_gen.get_slabs()
# slabs = generate_all_slabs(st, 1, 10, 10)
for slab in slabs:
    slab = center_slab(slab)
    if not slab.is_polar() and slab.is_symmetric():
        print(slab.miller_index, slab.shift, slab.is_polar(), slab.is_symmetric())
        slab.to("poscar", "slab_{}_{}.vasp".format(slab.miller_index, slab.shift))




#%%
from monty.serialization import loadfn
from qubitPack.tool_box import get_db
from pymatgen import Structure
import os

db = get_db("2dMaterials_from_icsd_v1", "VaspTasks_v1", password="qiminyan", user="readUser")
es = loadfn("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/tci_slab_nonhex/tci_nonhex_thinner_slabs.json")

spg_list = set([es[e]["spacegroup"] for e in es])
for spg in spg_list:
    for e in es:
        print(es)
        if es[e]["spacegroup"] == spg:
            if "/" in spg:
                dir_name = "-".join(spg.split("/"))
            else:
                dir_name = spg
            os.makedirs("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/tci_slab_nonhex/{}".format(dir_name), exist_ok=True)
            os.chdir("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/tci_slab_nonhex/{}".format(dir_name))
            # st = Structure.from_str(es[e]["five_layers"], fmt="cif")
            # st.to("poscar", "bulk_icsd_{}.vasp".format(e))

            ml = Structure.from_dict(db.collection.find_one({"icsd_id": e})["output"]["structure"])
            ml.to("poscar", "mono_icsd_{}.vasp".format(e))




#%%
from mpinterfaces.utils import ensure_vacuum
from monty.serialization import loadfn, dumpfn
ti = "ti"
db = get_db("2dMaterials_from_icsd_v1", "VaspTasks_v1", password="qiminyan", user="readUser")
es = loadfn("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/{}_slab_hex/{}_hex_thinner_slabs.json".format(ti, ti))

for e in es:
    ml = Structure.from_dict(db.collection.find_one({"icsd_id": e})["output"]["structure"])
    ml = ensure_vacuum(ml, 200)
    species, coords = ml.species, ml.cart_coords
    for l in [0,1,2,3]:
        for n, (specie, coor) in enumerate(zip(species, coords)):
            ml.append(specie, coor+[0,0,10*(l+1)], coords_are_cartesian=True)
            print(st.num_sites)
    ml = ensure_vacuum(ml, 15)
    es[e].update({"five_layers": ml.to(fmt="cif")})

dumpfn(es, "/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/{}_slab_hex/{}_hex_thinner_slabs.json".format(ti, ti), indent=4)

#%%
from monty.serialization import loadfn
from qubitPack.tool_box import get_db
from pymatgen import Structure
import os
ti = "ti"
db = get_db("2dMaterials_from_icsd_v1", "VaspTasks_v1", password="qiminyan", user="readUser")
es = loadfn("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/{}_slab_hex/{}_hex_thinner_slabs.json".format(ti, ti))

for e in es:
    os.makedirs("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/{}_slab_hex/five_layers".format(ti), exist_ok=True)
    os.chdir("/Users/jeng-yuantsai/Research/code/JPack/projects/ferromag_insulator_Lin/{}_slab_hex/five_layers".format(ti))
    st = Structure.from_str(es[e]["five_layers"], fmt="cif")
    st.to("poscar", "bulk_icsd_{}.vasp".format(e))


#%%
from qubitPack.tool_box import get_db

db = get_db("mag_insulator", "")