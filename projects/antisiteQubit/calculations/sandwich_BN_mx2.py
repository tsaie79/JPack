#%% wflow
from fireworks import LaunchPad, Workflow

from atomate.vasp.fireworks import OptimizeFW, StaticFW
from atomate.vasp.powerups import *
from atomate.vasp.fireworks.jcustom import *

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPHSEBSSet

from projects.antisiteQubit.wf_defect import Sandwich

import os
from monty.serialization import loadfn

CATEGORY = "sandwich_Ws"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/antisiteQubit/{}/my_launchpad.yaml".format(CATEGORY)))

# heters_vacu_44_antisite.vasp
defect_st = Structure.from_file("/home/tug03990/work/antisiteQubit/{}/stage_1_0_221antisite.vasp".format(CATEGORY))
host = Structure.from_file("/home/tug03990/work/antisiteQubit/{}/stage_1_0_221sc.vasp".format(CATEGORY))

wf = Sandwich(structure=defect_st).wf()
fws1 = wf.fws

wf = Sandwich(structure=host).wf()
fws2 = wf.fws

wf = Workflow(fws2+fws1)
wf.name = "BN-WS2-BN:PBEvdw:HSE_scf"
# wf = add_additional_fields_to_taskdocs(wf, {"NN": [56, 67, 62, 72], "defect_type": "dz=4.4"})

wf = set_execution_options(wf, category=CATEGORY)
wf = preserve_fworker(wf)


LPAD.add_wf(wf)



#%% antisiteQubit/sandwich_Mos
"""
build supercells
"""
from pymatgen import Structure
from pymatgen.analysis.substrate_analyzer import ZSLGenerator
import os

from mpinterfaces.interface import Interface
from mpinterfaces.transformations import *
from mpinterfaces.old_transformations import *
from mpinterfaces.utils import *

# p = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/sandwich_Mos/"
p = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/sandwich_Ws/"

# bn = slab_from_file([0,0,1],os.path.join(p, "modeling", "BN_hse_relax_pc.vasp"))
# mos2 = slab_from_file([0,0,1],os.path.join(p, "modeling", "Mo1 S2_4237.vasp"))

bn = slab_from_file([0,0,1],os.path.join(p, "modeling", "BN_hse_relax_pc.vasp"))
mos2 = slab_from_file([0,0,1],os.path.join(p, "modeling", "ws2_hse_relax_pc.vasp"))

substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    mos2,
    bn,
    max_area=200
)

substrate_slab_aligned.to("poscar", os.path.join(p,"modeling", "MoS2_exp.vasp"))
mat2d_slab_aligned.to("poscar", os.path.join(p,"modeling", "BN_exp.vasp"))

mis = get_mismatch(mat2d_slab_aligned.lattice.a, substrate_slab_aligned.lattice.a)


hetero_interfaces = generate_all_configs(mat2d_slab_aligned,
                                         substrate_slab_aligned,
                                         1, 1, 3.4)

for idx, i in enumerate(hetero_interfaces):
    i.to("poscar", os.path.join(p,"modeling", "stage_0_{}.vasp".format(idx)))

#%% above
"""
make sandwich
"""
# bn_exp = slab_from_file([0,0,1], os.path.join(p,"modeling", "BN_exp.vasp"))
# mos2_bn_exp = slab_from_file([0,0,1], os.path.join(p,"modeling", "0.vasp"))

bn = slab_from_file([0,0,1],os.path.join(p,"modeling", "BN_exp.vasp"))
mos2_bn = slab_from_file([0,0,1],os.path.join(p,"modeling", "stage_0_0.vasp"))
# bn_exp = Interface(bn, hkl=[0,0,1], from_ase=False)
# mos2_bn_exp = Interface(mos2_bn, hkl=[0,0,1], from_ase=False)

substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    bn,
    mos2_bn,
    max_area=200
)

hetero_interfaces = generate_all_configs(mat2d_slab_aligned, substrate_slab_aligned,
                                         1, 1, 3.4)


for idx, i in enumerate(hetero_interfaces):
    i.to("poscar", os.path.join(p, "modeling", "stage_1_{}.vasp".format(idx)))
# hetero_interfaces[1].to("poscar", os.path.join(p,"structures", "heter_bn_mos2_bn.vasp"))

#%% above
"""
make antisite structure
"""
import os
from pymatgen import Structure
from monty.serialization import dumpfn,loadfn

p = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/sandwich_Ws/"

host_pc = Structure.from_file(os.path.join(p, "modeling","stage_1_0.vasp"))

host_pc.make_supercell([2,2,1])
host_sc = host_pc

host_sc.to("poscar", os.path.join(p, "structures","stage_1_0_221sc.vasp"))

host_sc.replace(74, "W")
old_antisite = host_sc[74]

host_sc.sort()
host_sc.to("poscar", os.path.join(p, "structures","stage_1_0_221antisite.vasp"))

nn_include_self = host_sc.get_sites_in_sphere(old_antisite.coords, 3, include_index=True, include_image=True)

nn_idx = []
self_idx = []
for nn in nn_include_self:
    if old_antisite == nn[0]:
        self_idx += [nn[-2]]
        continue
    else:
        nn_idx.append(nn[-2])
    nn_idx += self_idx
dumpfn(nn_idx, os.path.join(p, "structures", "nn.json"))








#%% Read defect state

import os
from qubitPack.qc_searching.analysis.main import get_defect_state

proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/sandwich_Mos"
# proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/perturbed"
db_json = os.path.join(proj_path, "db.json")

# host_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_pc"
# host_path = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc"
host_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/sandwich_Mos"
db_host_json = os.path.join(host_path, "db.json")


perturb = None
tid = 540


tot, proj, d_df = get_defect_state(
    db_json,
    {"task_id": tid},
    0.2,
    -2.5,
    None,
    True,
    "dist",
    None,
    0.0
)
tot.to_clipboard()

#%% antisiteQubit/sandwich_Mos
"""
get PDOS on B & N
"""
from atomate.vasp.database import VaspCalcDb
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.inputs import Element

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/sandwich_Mos/db.json')
dos = db.get_dos(540)
dos_plt = DosPlotter(stack=False, sigma=0)
# tdos = DosPlotter(stack=True, zero_at_efermi=False)
# tdos.add_dos("tdos", dos)
# fig2 = tdos.get_plot(xlim=[-5,5])
# fig2.show()
dos_plt.add_dos("Total DOS", dos)

for element, dos in dos.get_element_dos().items():
    if element == Element("B"):
        dos_plt.add_dos("{}".format(element), dos)
    elif element == Element("N"):
        dos_plt.add_dos("{}".format(element), dos)
    elif element == Element("W"):
        dos_plt.add_dos("{}".format(element), dos)
    # elif element == Element("S"):

fig = dos_plt.get_plot(xlim=[-5,5])

fig.show()
# fig.savefig("/Users/jeng-yuantsai/Research/project/qubit/plt/sandwich_Mos.eps", format="eps")


