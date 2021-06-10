#%%
from projects.single_photon.wf import Defect
from fireworks import LaunchPad
from atomate.vasp.powerups import *
import os

category = "standard_defect"
wf = Defect.standard_defect(defect_type="substitutions", dtort=0)
wf = add_modify_incar(wf, {"incar_update": {"LMAXMIX": 4, "ENCUT":320}})
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

category = "soc_excited_state"
wf = Defect.soc_standard_defect(
    425, "/gpfs/scratch/tug03990/single_photon_emitter/excited_state/block_2021-02-18-21-01-22-891897/", scp=False,
    saxis=(0,0,1))

# db = get_db("single_photon_emitter", "standard_defect", port=12345)
# st = Structure.from_dict(db.collection.find_one({"task_id":22})["output"]["structure"])
# c = st[25].coords
# st.translate_sites(range(st.num_sites), -1*c, frac_coords=False)
# wf = write_PMGObjects(wf, dict(poscar=Poscar(st)))

wf = add_modify_incar(wf, {"incar_update": {"ISYM":3, "ENCUT":320, "LMAXMIX":4}})
wf = add_additional_fields_to_taskdocs(wf, {"PS": "singlet_excited_state"})
# wf = clear_to_db(wf)
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                         "my_launchpad.yaml".format(category))))

LPAD.add_wf(wf)
#%% Triplet
from projects.single_photon.wf import Defect
from fireworks import LaunchPad
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from qubitPack.tool_box import get_db
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Poscar
import os
# now, saxis=110
category = "soc_excited_state"
for i,j, ir in zip(["656*1 1*1 1*0 1*0 1*1 1*0 1*0 138*0", "656*1 1*0 1*1 1*1 1*0 1*0 1*0 138*0",
                    "656*1 1*1 1*0 1*1 1*0 1*0 1*0 138*0", "656*1 1*0 1*1 1*0 1*1 1*0 1*0 138*0"],
               ["100100", "011000", "101000", "010100"], ["g3", "g1", "g2", "g1"]):
    wf = Defect.soc_cdft(
        "B",
        399,
        800,
        i,
        None,
        None,
        secondary=False,
        category=category,
    )
    wf.name = wf.name+":ea:{}:{}".format(j, ir)
    wf = add_modify_incar(wf, {"incar_update":{"LMAXMIX":4, "ENCUT":320}})
    wf = clear_to_db(wf)
    wf = add_additional_fields_to_taskdocs(wf, {"excite_config": j, "ks_config": "ea", "ir_c3":ir, "PS":"ZPL"})
    LPAD = LaunchPad.from_file(
        os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                             "my_launchpad.yaml".format(category))))
    LPAD.add_wf(wf)
#%% Triplet testing
from projects.single_photon.wf import Defect
from fireworks import LaunchPad
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from qubitPack.tool_box import get_db
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Poscar
import os
# now, saxis=110
category = "soc_excited_state"
for i,j, ir in zip(["656*1 1*0.5 1*0.5 1*1 1*0 1*0 1*0 138*0", "656*1 1*0.5 1*0.5 1*0 1*1 1*0 1*0 138*0"],
                   ["551000", "550100"], ["g6", "g5"]):
    wf = Defect.soc_cdft(
        "B",
        399,
        800,
        i,
        None,
        None,
        secondary=False,
        category=category,
    )
    wf.name = "ea:{}:{}:".format(j, ir)+wf.name
    wf = add_modify_incar(wf, {"incar_update":{"LMAXMIX":4, "ENCUT":320}})
    wf = clear_to_db(wf)
    wf = add_additional_fields_to_taskdocs(wf, {"excite_config": j, "ks_config": "ea", "ir_c3":ir, "PS":"ZPL"})
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

category = "excited_state"
wf = Defect.cdft(
    task="B",
    std_d_tkid=400,
    nbands=400,
    std_d_base_dir="/gpfs/scratch/tug03990/single_photon_emitter/standard_defect/block_2021-01-11-17-39-35-861559/",
    category=category,
    occ={"1":"328*1 1*0 1*0 1*1 69*0", "-1":"329*1 71*0"}
)

wf = add_modify_incar(wf, {"incar_update": {"LMAXMIX":4, "ENCUT":320, "LWAVE":True, "LCHARG": True, "ISYM":3}})
# wf = add_modify_incar(wf, {"incar_update": {"LCHARG":True}}, fw_name_constraint="CDFT-C-HSE_relax")
wf = add_additional_fields_to_taskdocs(wf, {"spin":"singlet", "PS":"nonsoc"})
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