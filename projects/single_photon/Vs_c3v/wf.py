#%%
from projects.single_photon.wf import Defect
from fireworks import LaunchPad
from atomate.vasp.powerups import *
import os

category = "standard_defect"
wf = Defect.standard_defect(defect_type="vacancies", dtort=0)
# wf = add_modify_incar(wf, {"incar_update": {"ENCUT":320, "NBANDS": 950}})
wf = set_execution_options(wf, category=category)
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                         "my_launchpad.yaml".format(category))))
LPAD.add_wf(wf)

#%%
from projects.single_photon.wf import Defect
from fireworks import LaunchPad
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from qubitPack.tool_box import get_db
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Poscar
import os

category = "soc_standard_defect"
wf = Defect.soc_standard_defect(
    374, "/gpfs/scratch/tug03990/single_photon_emitter/standard_defect/block_2021-01-11-17-39-35-861559/", scp=True)

# db = get_db("single_photon_emitter", "standard_defect", port=12345)
# st = Structure.from_dict(db.collection.find_one({"task_id":22})["output"]["structure"])
# c = st[25].coords
# st.translate_sites(range(st.num_sites), -1*c, frac_coords=False)
# wf = write_PMGObjects(wf, dict(poscar=Poscar(st)))

wf = add_modify_incar(wf, {"incar_update": {"ISYM":3, "ENCUT":320}})
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                         "my_launchpad.yaml".format(category))))

LPAD.add_wf(wf)
#%%
from projects.single_photon.wf import Defect
from fireworks import LaunchPad
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from qubitPack.tool_box import get_db
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Poscar
import os

category = "soc_cdft"
wf = Defect.soc_cdft(
    "B-C-D",
    375,
    950,
    "/home/tug03990/scratch/single_photon_emitter/soc_standard_defect/block_2021-01-22-20-56-50-305107",
    "/home/tug03990/scratch/single_photon_emitter/standard_defect/block_2021-01-11-17-39-35-861559",
    secondary=False
)
wf = add_modify_incar(wf, {"incar_update":{"ENCUT":320}})
wf = add_additional_fields_to_taskdocs(wf, {"PS": "direct"})
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                         "my_launchpad.yaml".format(category))))
LPAD.add_wf(wf)

#%% CDFT C fails
from projects.single_photon.wf import Defect
from fireworks import LaunchPad
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from qubitPack.tool_box import get_db
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Poscar
import os, glob

category = "soc_cdft"
sts_p = '/home/tug03990/work/single_photon_emitter/soc_cdft/block_2021-01-20-04-45-04-585710/launcher_2021-02-09-22-52-38-100010/inte'
fws = []
for st in glob.glob(os.path.join(sts_p, "*.vasp")):
    wf = Defect.soc_cdft(
        "B",
        354,
        950,
        "/home/tug03990/work/single_photon_emitter/soc_standard_defect/block_2021-01-17-04-46-01-108825",
        "/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869,",
        specific_poscar=Poscar.from_file(st),
    )
    wf = add_additional_fields_to_taskdocs(wf, {"poscar_idx": st, "start":"2_3"})
    fws.extend(wf.fws)

wf = Workflow(fws, name="MoS2-cdft-B-to-C-interpolate")
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                         "my_launchpad.yaml".format(category))))
LPAD.add_wf(wf)
#%%
from qubitPack.tool_box import get_db
from qubitPack.qc_searching.analysis.main import get_defect_state

tot, proj, d_df = get_defect_state(
    get_db("single_photon_emitter", "soc_standard_defect"),
    {"task_id": 375},
    1,-2.5,
    None,
    True,
    "dist",
    None,
    0.1)



