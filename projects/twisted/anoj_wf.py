from atomate.vasp.fireworks.core import ScanOptimizeFW, StaticFW, NonSCFFW, OptimizeFW
from atomate.vasp.workflows.presets.core import wf_bandstructure
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    add_modify_kpoints,
    set_queue_options,
    set_execution_options,
    clean_up_files,
    modify_gzip_vasp
)
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet, MPScanStaticSet

from fireworks import LaunchPad, Workflow

from qubitPack.tool_box import phonopy_structure
import os
import numpy as np

path = "structures"
st = Structure.from_file(os.path.join(path,"POSCAR_4.0.vasp"))

st, site_info = phonopy_structure(st)

def dipole(structure):
    weights = [s.species.weight for s in structure]
    center_of_mass = np.average(
        structure.frac_coords, weights=weights, axis=0
    )
    return {"IDIPOL": 3, "LDIPOL": True, "DIPOL": center_of_mass}

encut = 1.3*max([potcar.enmax for potcar in MPRelaxSet(st).potcar])

uis ={
         "ENCUT": encut,
         "EDIFF": 1e-4,
         "EDIFFG": -0.05,
         "METAGGA": "R2SCAN",
         "ISIF": 2,
         "LCHARG": False,
         "LAECHG": False,
         "LELF": False,
         "NELM": 250,
     } # "ISMEAR": 0, "SIGMA": 0.05}
# uis.update(dipole(st))
scan_relax_incar = MPScanStaticSet(st, user_incar_settings=uis, force_gamma=True, bandgap=0).incar
print(scan_relax_incar)


for pop_input in ["ENAUG", "MAGMOM", "KSPACING"]:
    scan_relax_incar.pop(pop_input)

fw = OptimizeFW(structure=st, override_default_vasp_params=dict(user_incar_settings=scan_relax_incar),
                job_type="normal")
wf = Workflow([fw], name="{}:{}:metal:no_magmom".format(st.formula, st.num_sites))

wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"MAGMOM":"1*{}".format(st.num_sites)}})

cat = "anoj"
wf = set_execution_options(wf, category=cat, fworker_name="gpu_nersc")
wf = preserve_fworker(wf)

lpad = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/twisted/{}/my_launchpad.yaml".format(cat))))
# lpad.add_wf(wf)
