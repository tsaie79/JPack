#%% db
from atomate.vasp.database import VaspCalcDb

db_dict = {
    "host": "localhost",
    "port": 1234,
    "database": "2dMat_from_cmr_fysik",
    "collection": "standard_defect",
    "admin_user": "adminUser",
    "admin_password": "qiminyan",
    "readonly_user": "readUser",
    "readonly_password": "qiminyan",
    "aliases": {}
}
def db(collection):
    return VaspCalcDb(host="localhost", port=1234, database="2dMat_from_cmr_fysik",
                      collection=collection, user="readUser", password="qiminyan", authsource="2dMat_from_cmr_fysik")
#%% get pc from c2db
e = db("2dMaterial_v1").collection.find_one({"uid":"WSe2-MoS2-NM"})

from pymatgen import Structure
pc = Structure.from_dict(e["structure"])
pc.to("poscar", "/Users/jeng-yuantsai/Research/project/twisted/calculations/model/WSe2_c2db.vasp")

#%% rot pc
from pymatgen.core.operations import SymmOp
twist_op = SymmOp.from_axis_angle_and_translation(axis=[0,0,1],
                                                  angle=5,
                                                  angle_in_radians=False,
                                                  translation_vec=(
                                                      0, 0, 1))
pc.apply_operation(twist_op)
pc.to("poscar", "/Users/jeng-yuantsai/Research/project/twisted/calculations/model/pc_rot.vasp")



#%% hetero
from pymatgen import Structure, Molecule
from pymatgen.transformations.standard_transformations import RotationTransformation

import os

from mpinterfaces.interface import Ligand, Interface
from mpinterfaces.old_transformations import *
from mpinterfaces.utils import *


rot = Structure.from_file("/Users/jeng-yuantsai/Research/project/twisted/calculations/model/pc_rot.vasp")
rot = Interface(rot, hkl=[0,0,1], min_thick=1, min_vac=40, primitive=False, from_ase=True)

no_rot = slab_from_file([0,0,1], "/Users/jeng-yuantsai/Research/project/twisted/calculations/model/WSe2_c2db.vasp")
no_rot = Interface(no_rot, hkl=[0,0,1], min_thick=1, min_vac=40, primitive=False, from_ase=True)
#
# AlN = Structure.from_file('/Users/jeng-yuantsai/Research/project/MgB2/model/AlN/AlN_mp-661_computed.cif')
# AlN = Interface(AlN, hkl=[0,0,1], min_thick=1, min_vac=40, primitive=False, from_ase=True)

# graphene = Structure.from_file('/Users/jeng-yuantsai/Research/project/MgB2/model/new/hetero/graphene/C_mp-568806_computed.cif')
# graphene = Interface(graphene, hkl=[0,0,1], min_thick=1, min_vac=40, primitive=False, from_ase=True)

# mgb2_slab.to("poscar", os.path.join(p, "mgb2_slab.vasp"))
# sic_slab.to("poscar", os.path.join(p, "sic_slab.vasp"))

substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    rot,
    no_rot,
    max_area=20,
    max_angle_diff=0,
    r1r2_tol=0.00001,
)
# s = "/Users/jeng-yuantsai/Research/project/MgB2/model/new/hetero/si_terminate_top"
# s = '/Users/jeng-yuantsai/Research/project/MgB2/model/new/hetero/graphene'
s = "/Users/jeng-yuantsai/Research/project/twisted/calculations/model/test"
substrate_slab_aligned.to("poscar", os.path.join(s, "rot.vasp"))
mat2d_slab_aligned.to("poscar", os.path.join(s, "no_rot.vasp"))

hetero_interfaces = generate_all_configs(mat2d_slab_aligned, substrate_slab_aligned,
                                         nlayers_substrate=1, nlayers_2d=1, seperation=2)
for idx, i in enumerate(hetero_interfaces):
    i.to("poscar", os.path.join(s, "{}.vasp".format(idx)))

#%% WF
import numpy as np
from pymatgen.io.vasp.sets import MPRelaxSet
from atomate.vasp.fireworks import OptimizeFW, StaticFW
from atomate.vasp.powerups import *
from fireworks import Workflow, LaunchPad, FileTransferTask
import os

CATEGORY = "hetero"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/mgb2/{}/my_launchpad.yaml".format(CATEGORY)))
#%% db
from atomate.vasp.database import VaspCalcDb

db_dict = {
    "host": "localhost",
    "port": 1234,
    "database": "2dMat_from_cmr_fysik",
    "collection": "standard_defect",
    "admin_user": "adminUser",
    "admin_password": "qiminyan",
    "readonly_user": "readUser",
    "readonly_password": "qiminyan",
    "aliases": {}
}
def db(collection):
    return VaspCalcDb(host="localhost", port=1234, database="2dMat_from_cmr_fysik",
                      collection=collection, user="readUser", password="qiminyan", authsource="2dMat_from_cmr_fysik")
#%% get pc from c2db
e = db("2dMaterial_v1").collection.find_one({"uid":"WSe2-MoS2-NM"})

from pymatgen import Structure
pc = Structure.from_dict(e["structure"])
print(pc.lattice.abc)
pc.to("poscar", "/Users/jeng-yuantsai/Research/project/twisted/calculations/model/WSe2_c2db.vasp")

#%% build twisted bilayer
"""
WARNING: the bond lengths are wrong!
Correct bond lengths by changing lattice constants (here, change only lattice vector C)
"""
import supercell_core as sc
import numpy as np

# graphene = sc.read_POSCAR("/Users/jeng-yuantsai/Research/project/twisted/calculations/model/WSe2_c2db.vasp")
graphene = sc.read_POSCAR("/Users/jeng-yuantsai/Research/project/twisted/calculations/model/Bi2Se3/layer.vasp")

h = sc.heterostructure().set_substrate(graphene).add_layer(graphene)

# res = h.opt(max_el=20, thetas=[np.arange(0.1*sc.DEGREE, 30*sc.DEGREE, 0.1*sc.DEGREE)])
res = h.opt(max_el=20, thetas=[[5.1*sc.DEGREE]]) #5.1
print("{}, {}".format(res.thetas()[0]*1/sc.DEGREE, res.max_strain()))

# Save supercell to VASP POSCAR
res.superlattice().save_POSCAR("/Users/jeng-yuantsai/Research/project/twisted/calculations/model/Bi2Se3/test/POSCAR_1.vasp")



#%% setup interlayer separation and vacuum
from pymatgen import Structure
from qubitPack.tool_box import modify_vacuum
dx = 1.64
rot = Structure.from_file("/Users/jeng-yuantsai/Research/project/twisted/calculations/model/Bi2Se3/test/2.vasp")
sites = []
for idx, i in enumerate(rot.cart_coords):
    print(i)
    if i[2] < 31.12 and i[2] > 24:
        sites.append(idx)
print(sites)

rot.translate_sites(sites, [0,0,-dx], frac_coords=False)
# st = Structure(rot.lattice, rot.species, sites)
# rot = modify_vacuum(rot, 40)
rot.to("poscar", "/Users/jeng-yuantsai/Research/project/twisted/calculations/model/Bi2Se3/test/2_sh.vasp")

#%% Bi2Se3

