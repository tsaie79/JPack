from atomate.vasp.workflows.base.adsorption import get_wf_slab
from atomate.vasp.powerups import *

from fireworks import LaunchPad, Workflow

from pymatgen import Structure, Molecule
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.core.surface import SlabGenerator

import os


CATEGORY = "H_ads"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/CO2_hydrogenation/{}/my_launchpad.yaml".format(CATEGORY))))


Rh = Structure.from_file("/home/tug03990/config/project/CO2_hydrogenation/H_ads/structures/"
                         "Rh_mp-74_conventional_standard.cif")
slab_100 = SlabGenerator(Rh, [1, 0, 0], 10, 15, center_slab=True).get_slab()

H = Molecule(["H", "H"], [[0,0,0], [0.74, 0, 0]])

wf = get_wf_slab(slab_100, False, [H], dict(find_args={"distance": 1.6}, repeat=[3,3,1]), add_molecules_in_box=True)

for fw in wf.fws:
    # print(fw.name)
    if not "structure optimization" in fw.name:
        print(fw.name)
        wf = add_modify_incar(
            wf,
            {"incar_update":
                 {
                     "EDIFFG": -0.01,
                     "EDIFF":1E-5,
                     "ENCUT": 500,
                     "ISIF":2,
                     "LVHAR":True,
                     "LVTOT":False,
                     "ISMEAR":1,
                     "IBRION":2,
                     "POTIM":0.1,
                     "SIGMA":0.1,
                     "LREAL": False
                 }
            },
            fw_name_constraint=fw.name)
for fw in wf.fws:
    st = fw.tasks[0]["structure"]
    st.to("POSCAR", "/home/tug03990/h2_rh/poscar_{}.vasp".format(fw.name))


wf = add_modify_incar(wf, {"incar_update": {"IVDW": 1}}) #DFT-D2 vdw correction
wf = preserve_fworker(wf)
wf = set_execution_options(wf, category=CATEGORY)
wf = add_modify_incar(wf)

# LPAD.add_wf(wf)

#%% Run lobster for COHP
from fireworks import LaunchPad
from atomate.vasp.workflows.base.lobster import get_wf_lobster
from atomate.vasp.powerups import *
from qubitPack.tool_box import get_db
import numpy as np
import os
from pymatgen import Structure

CATEGORY = "H_ads"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/CO2_hydrogenation/{}/my_launchpad.yaml".format(CATEGORY))))

d = get_db("CO2_hydrogenation", "H_ads", port=12345)




for tkid, name, n in zip([23], ["top"], range(1)):
    e = d.collection.find_one({"task_id": tkid})
    st = Structure.from_dict(e["output"]["structure"])

    weights = [s.species.weight for s in st]
    center_of_mass = np.average(
        st.frac_coords, weights=weights, axis=0
    )
    dipole = {"LDIPOL": True, "IDIPOL": 3, "DIPOL": center_of_mass, "IVDW": 1}
    uis = dipole.copy()
    uis.update({"ENCUT":500})

    wf = get_wf_lobster(st, user_incar_settings=uis)
    wf.name = "{}:{}:".format(st.formula, name)+wf.name

    wf = add_additional_fields_to_taskdocs(wf, {"ads_site": name})
    wf = preserve_fworker(wf)
    wf = set_execution_options(wf, category=CATEGORY)
    wf = add_modify_incar(wf)
    LPAD.add_wf(wf)
