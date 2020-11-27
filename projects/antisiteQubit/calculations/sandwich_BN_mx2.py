#%% wflow
from fireworks import LaunchPad, Workflow

from atomate.vasp.fireworks import OptimizeFW, StaticFW
from atomate.vasp.powerups import *
from atomate.vasp.fireworks.jcustom import *

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPHSEBSSet

from projects.antisiteQubit.wf_defect import Sandwich

import os

CATEGORY = "sandwich_Mos"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/antisiteQubit/{}/my_launchpad.yaml".format(CATEGORY)))

defect_st = Structure.from_file("/home/tug03990/work/antisiteQubit/sandwich_Mos/mo_s_331.vasp")

Sandwich(structure=defect_st, lpad=LPAD, category=CATEGORY).wf()



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

p = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/sandwich_Mos/"

# bn = Structure.from_file(os.path.join(p, "structures", "BN_hse_relax_pc.vasp"))
# mos2 = Structure.from_file(os.path.join(p, "structures", "Mo1 S2_4237.vasp"))

bn = slab_from_file([0,0,1], os.path.join(p, "structures", "BN_hse_relax_pc.vasp"))
mos2 = slab_from_file([0,0,1], os.path.join(p, "structures", "Mo1 S2_4237.vasp"))

substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    mos2,
    bn,
    max_area=200
)

substrate_slab_aligned.to("poscar", os.path.join(p,"structures", "MoS2_exp.vasp"))
mat2d_slab_aligned.to("poscar", os.path.join(p,"structures", "BN_exp.vasp"))

mis = get_mismatch(mat2d_slab_aligned.lattice.a, substrate_slab_aligned.lattice.a)


hetero_interfaces = generate_all_configs(mat2d_slab_aligned,
                                         substrate_slab_aligned,
                                         1, 1, 3.4)


for idx, i in enumerate(hetero_interfaces):
    i.to("poscar", os.path.join(p,"structures", "{}.vasp".format(idx)))

#%% above
"""
make sandwich
"""
bn_exp = slab_from_file([0,0,1], os.path.join(p,"structures", "BN_exp.vasp"))
mos2_bn_exp = slab_from_file([0,0,1], os.path.join(p,"structures", "0.vasp"))

substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    bn_exp,
    mos2_bn_exp,
    max_area=200
)

hetero_interfaces = generate_all_configs(mat2d_slab_aligned, substrate_slab_aligned,
                                         1, 1, 3.4)

hetero_interfaces[0].to("poscar", os.path.join(p,"structures", "heter_bn_mos2_bn.vasp"))

#%% Read defect state

import os
from qubitPack.qc_searching.analysis.main import get_defect_state

proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/sandwich_Mos"
# proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/perturbed"
db_json = os.path.join(proj_path, "db.json")

host_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_pc"
# host_path = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc"
db_host_json = os.path.join(host_path, "db.json")



perturb = None
tid = 502

tot, proj, d_df = get_defect_state(
    db_json,
    {"task_id": tid},
    2.656,
    -0.1,
    None,
    True,
    "dist",
    db_host_json,
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
dos = db.get_dos(502)
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
    elif element == Element("Mo"):
        dos_plt.add_dos("{}".format(element), dos)
    # elif element == Element("S"):
fig = dos_plt.get_plot(xlim=[-5,5])

fig.show()
fig.savefig("/Users/jeng-yuantsai/Research/project/qubit/plt/sandwich_Mos.eps", format="eps")



