#%% WX2 soc gw test
from pymatgen import Structure
from pymatgen.io.vasp.sets import MVLGWSet

from atomate.vasp.workflows.jcustom.gw_workflow import gw_wf
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *

from fireworks import LaunchPad

import os
import numpy as np

cat = "soc_gw"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "atomate/config/project/single_photon_emitter/{}/my_launchpad.yaml".format(cat))
)
scf_path = "/global/project/projectdirs/m2663/tsai/single_photon_emitter/standard_defect/launcher_2021-01-13-04-09-57-563577"
st = Structure.from_file(os.path.join(scf_path, "POSCAR.gz"))
wf = gw_wf(structure=st, ncores=1, vis_static=None, vasp_input_set_params=None, vasptodb={"category":cat},
           wf_addition_name=None, prev_dir=scf_path, nbands_factor=5)

# related to directory

# wf = add_modify_incar(wf,
#                       {
#                           "incar_update":
#                               {
#                                    "ISMEAR": 0,
#                                    "SIGMA": 0.01,
#                                    "ENCUT": np.ceil(1.3*max([potcar.enmax for potcar in MVLGWSet(structure).potcar]))
#                               }
#                       })

wf = clear_to_db(wf)
wf = add_modify_incar(wf, {"incar_update": {"ICHARG":0}}, fw_name_constraint=wf.fws[0].name)
wf = add_modify_incar(wf)
wf = set_execution_options(wf, category=cat)
wf = preserve_fworker(wf)
LPAD.add_wf(wf)


