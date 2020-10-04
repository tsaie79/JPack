from pymatgen import Structure
from pymatgen.io.vasp.sets import MVLGWSet

from atomate.vasp.workflows.jcustom.gw_workflow import gw_wf
from atomate.vasp.powerups import *

from fireworks import LaunchPad

import os
import numpy as np

def gw_test(structure, ncores, cat):
    lpad = LaunchPad.from_file(
        os.path.join(os.path.expanduser("~"), "config/project/antisiteQubit/{}/my_launchpad.yaml".format(cat))
    )
    wf = gw_wf(structure, ncores, vis_static=None, vasp_input_set_params=None, vasptodb={"category":cat},
               wf_addition_name=None)
    # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="SCAN_relax")
    # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="SCAN_scf")
    # related to directory

    wf = add_modify_incar(wf,
                          {
                              "incar_update":
                                  {
                                       "ISMEAR": 0,
                                       "SIGMA": 0.01,
                                       "ENCUT": np.ceil(1.3*max([potcar.enmax for potcar in MVLGWSet(structure).potcar]))
                                  }
                          })
    wf = set_execution_options(wf, category=cat)
    wf = preserve_fworker(wf)
    lpad.add_wf(wf)

gw_test(Structure.from_file("/home/tug03990/work/antisiteQubit/Ws_gw/Ws_no_opt.vasp"), 80, "Ws_gw")

