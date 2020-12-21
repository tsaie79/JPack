#%% hetero
from pymatgen import Structure, Molecule
from pymatgen.transformations.standard_transformations import RotationTransformation

import os

from mpinterfaces.interface import Ligand, Interface
from mpinterfaces.old_transformations import *
from mpinterfaces.utils import *

p = "/Users/jeng-yuantsai/Research/project/MgB2/model/new"

sic = Structure.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/SiC_6H/SiC_mp-7631_computed.cif")
sic.remove_sites([5, 11])
sic.to("poscar", "/Users/jeng-yuantsai/Research/project/MgB2/model/SiC_6H/SiC_mp-7631_computed_cut.vasp")

# si_terminate
sic = Structure(sic.lattice, sic.species, [np.multiply(site.frac_coords, [1,1,-1]) for site in sic.sites])

sic_slab = Interface(sic, hkl=[0,0,1], min_thick=10, min_vac=40, primitive=False, from_ase=True)

# mgb2_slab = slab_from_file([0,0,1], "/Users/jeng-yuantsai/Research/project/MgB2/model/MgB2/MgB2_mp-763_computed.cif")
# mgb2_slab = slab_from_file([0,0,1], '/Users/jeng-yuantsai/Research/project/MgB2/model/new/mgb2_112.vasp')
mgb2 = Structure.from_file('/Users/jeng-yuantsai/Research/project/MgB2/model/new/mgb2_112.vasp')
mgb2 = Structure.from_file('/Users/jeng-yuantsai/Research/project/MgB2/model/new/hetero/si_terminate_b/mat2d.vasp')
# mgb2 = Structure(mgb2.lattice, mgb2.species, [np.multiply(site.frac_coords, [1,1,-1]) for site in mgb2.sites])
mgb2_slab = Interface(mgb2, hkl=[0,0,1], min_thick=1, min_vac=40, primitive=False, from_ase=True)



# mgb2_slab.to("poscar", os.path.join(p, "mgb2_slab.vasp"))
# sic_slab.to("poscar", os.path.join(p, "sic_slab.vasp"))

substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    sic_slab,
    mgb2_slab,
    max_area=200
)
s = "/Users/jeng-yuantsai/Research/project/MgB2/model/new/hetero/si_terminate_b"
substrate_slab_aligned.to("poscar", os.path.join(s, "subs.vasp"))
mat2d_slab_aligned.to("poscar", os.path.join(s, "mat2d.vasp"))

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

terminate = "c_terminate"
st_p = "/home/tug03990/work/mgb2/hetero"

heter = Structure.from_file(os.path.join(st_p, terminate, "0.vasp"))
mgb2 = Structure.from_file(os.path.join(st_p, terminate, "mat2d.vasp"))
sic = Structure.from_file(os.path.join(st_p, terminate, "subs.vasp"))

def dp(structure):
    weights = [s.species.weight for s in structure]
    center_of_mass = np.average(
        structure.frac_coords, weights=weights, axis=0
    )

    return {"IDIPOL":3, "LDIPOL": True, "DIPOL": center_of_mass}


uis = {"ISIF": 2, "EDIFFG":-0.02}
uis.update(dp(heter))

opt_vdw = OptimizeFW(heter, override_default_vasp_params={"user_incar_settings": uis})
uis = {}
uis.update({"EDIFF":1E-5})
static = StaticFW(heter, parents=opt_vdw, vasp_input_set_params={"user_incar_settings": uis})

hetero_fws = [opt_vdw, static]
wf = Workflow(hetero_fws, name=terminate)

add_additional_fields_to_taskdocs(wf, {"terminate": terminate})
wf = add_modify_incar(wf, {"incar_update": dict(NUPDOWN=0, MAGMOM=MPRelaxSet(heter).incar.get("MAGMOM"))})
wf = preserve_fworker(wf)
wf = set_execution_options(wf, category=CATEGORY)

print(wf)
# LPAD.add_wf(wf)


for st in [sic]:
    fws = []
    uis = {"EDIFF":1E-5, "EDIFFG":-0.02, "ISIF":2}
    uis.update(dp(st))
    opt = OptimizeFW(st, override_default_vasp_params={"user_incar_settings":uis})

    uis = {"EDIFF":1E-6}
    uis.update(dp(st))
    static = StaticFW(st, parents=opt, vasp_input_set_params={"user_incar_settings": uis})

    fws.extend([static, opt])
    wf = Workflow(fws, name=terminate+":{}".format(st.formula))

    add_additional_fields_to_taskdocs(wf, {"terminate": terminate})
    wf = add_modify_incar(wf, {"incar_update": dict(NUPDOWN=0, MAGMOM=MPRelaxSet(st).incar.get("MAGMOM"))})
    wf = preserve_fworker(wf)
    wf = set_execution_options(wf, category=CATEGORY)

    print(wf)
    LPAD.add_wf(wf)


#%% check locpot
from atomate.vasp.database import VaspCalcDb
from matplotlib import pyplot as plt
import os

CATEGORY = "hetero"
DB = VaspCalcDb.from_db_file(
    os.path.join(os.path.expanduser("~"), "config/project/mgb2/{}/db.json".format(CATEGORY)))

#locpot
# e = DB.collection.find_one({"task_id":129})
# locpot_z = e["calcs_reversed"][0]["output"]["locpot"]["2"]
# plt.plot(locpot_z)
# plt.show()

def energy(tids):
    sic = DB.collection.find_one({"task_id":tids[0]})["output"]["energy"]
    mgb2 = DB.collection.find_one({"task_id":tids[1]})["output"]["energy"]
    sic_mgb2 = DB.collection.find_one({"task_id":tids[2]})["output"]["energy"]
    dE = sic_mgb2-(sic+mgb2)
    return sic, mgb2, sic_mgb2, dE

#si
si = energy([141, 137, 145])
#c
c = energy([140, 137, 149])



print(si,c)
