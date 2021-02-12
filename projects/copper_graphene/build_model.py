#%% hetero
from pymatgen import Structure, Molecule
from pymatgen.transformations.standard_transformations import RotationTransformation

import os

from mpinterfaces.interface import Ligand, Interface
from mpinterfaces.old_transformations import *
from mpinterfaces.utils import *

p = "/Users/jeng-yuantsai/Research/project/copper_graphene/model"

graphene = Structure.from_file("/Users/jeng-yuantsai/Research/project/copper_graphene/model/C_mp-568806_conventional_standard.cif")

graphene = Interface(graphene, hkl=[0,0,1], min_thick=1, min_vac=40, primitive=False, from_ase=True)

copper = Structure.from_file("/Users/jeng-yuantsai/Research/project/copper_graphene/model/Cu_mp-30_conventional_standard.cif")
copper = Interface(copper, hkl=[1,1,1], min_thick=5, min_vac=40, primitive=False, from_ase=True)

substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    graphene,
    copper,
    max_area=200
)

substrate_slab_aligned.to("poscar", os.path.join(p, "subs.vasp"))
mat2d_slab_aligned.to("poscar", os.path.join(p, "mat2d.vasp"))

hetero_interfaces = generate_all_configs(mat2d_slab_aligned, substrate_slab_aligned,
                                         nlayers_substrate=1, nlayers_2d=1, seperation=3.26)
for idx, i in enumerate(hetero_interfaces):
    i.to("poscar", os.path.join(p, "interface", "{}.vasp".format(idx)))


#%% hetero
from pymatgen import Structure, Molecule
from pymatgen.transformations.standard_transformations import RotationTransformation

import os

from mpinterfaces.interface import Ligand, Interface
from mpinterfaces.old_transformations import *
from mpinterfaces.utils import *

p = "/Users/jeng-yuantsai/Research/project/copper_graphene/model"

copper = Structure.from_file("/Users/jeng-yuantsai/Research/project/copper_graphene/model/Cu_mp-30_conventional_standard.cif")
copper = Interface(copper, hkl=[1,1,1], min_thick=5, min_vac=35, primitive=False, from_ase=True, center_slab=True)
#30-2.5-0.05
copper_gr = Structure.from_file("/Users/jeng-yuantsai/Research/project/copper_graphene/model/interface/gr_cop.vasp")
# copper_gr = Structure(copper_gr.lattice, copper_gr.species, copper_gr.frac_coords*[1,1,-1])
copper_gr = Interface(copper_gr, hkl=[0,0,1], min_thick=1, min_vac=5, primitive=False, from_ase=True, center_slab=True)

substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    copper,
    copper_gr,
    max_area=200
)

# substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
#     copper,
#     copper_gr,
#     max_area=200
# )


substrate_slab_aligned.to("poscar", os.path.join(p, "subs.vasp"))
mat2d_slab_aligned.to("poscar", os.path.join(p, "mat2d.vasp"))

hetero_interfaces = generate_all_configs(mat2d_slab_aligned, substrate_slab_aligned,
                                         nlayers_substrate=1, nlayers_2d=3, seperation=3.26)
for idx, i in enumerate(hetero_interfaces):
    i.to("poscar", os.path.join(p, "interface", "{}.vasp".format(idx)))

#%%
from mpinterfaces.interface import Ligand, Interface
from mpinterfaces.old_transformations import *
from mpinterfaces.utils import *

p = "/Users/jeng-yuantsai/Research/project/copper_graphene/model"

copper = Structure.from_file("/Users/jeng-yuantsai/Research/project/copper_graphene/model/Cu_mp-30_conventional_standard.cif")
copper = Interface(copper, hkl=[1,1,1], min_thick=5, min_vac=20, primitive=False, from_ase=True, center_slab=True)

copper.to("poscar", "/Users/jeng-yuantsai/Research/project/copper_graphene/model/3-Cu/3-Cu.vasp")

#%% random sampling
from random import sample


rnd_site = sample(range(192), 64)
for site in rnd_site:
    st.replace(site, "C")



#%%
from qubitPack.tool_box import modify_vacuum
p = "/Users/jeng-yuantsai/Research/project/copper_graphene/model"
st = Structure.from_file(os.path(p, "interface", "gr_cop.vasp"))
st = modify_vacuum(st, )

#%%
from pymatgen import Structure
from qubitPack.tool_box import get_db

db = get_db("copper_graphene", "modeling")

# for e in db.collection.find({"stacking":"head", "input.incar.TEEND":{"$exists":1}}):
#     Structure.from_dict(e["output"]["structure"]).to("poscar","/Users/jeng-yuantsai/"
#                                                               "Research/project/copper_graphene/"
#                                                               "calculations/structures/stacking_head/{}K.vasp".format(e["input"]["incar"]["TEEND"]))

for e in db.collection.find({"formula_pretty":"Cu3C2", "input.incar.TEEND":{"$exists":1}}):
    Structure.from_dict(e["output"]["structure"]).to("poscar","/Users/jeng-yuantsai/"
                                                              "Research/project/copper_graphene/"
                                                              "calculations/structures/c-cu/{}K.vasp".format(e["input"]["incar"]["TEEND"]))
#%% Cu
from pymatgen import Structure
from qubitPack.tool_box import get_db

db = get_db("copper_graphene", "modeling")

for e in db.collection.find({"chemsys":"Cu", "nsites":384, "input.incar.TEEND":{"$exists":1}}):
    Structure.from_dict(e["output"]["structure"]).to(
        "poscar",
        "/Users/jeng-yuantsai/"
        "Research/project/copper_graphene/"
        "calculations/structures/cu/{}K.vasp".format(e["input"]["incar"]["TEEND"])
    )