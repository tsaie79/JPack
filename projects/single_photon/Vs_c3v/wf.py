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
    354,
    950,
    "/home/tug03990/work/single_photon_emitter/soc_standard_defect/block_2021-01-17-04-46-01-108825",
    "/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869,"
)
wf = add_additional_fields_to_taskdocs(wf, {"PS": "secondary"})
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
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os
from pymatgen import Structure
d = get_db("single_photon_emitter", "soc_cdft", port=12345)
es = d.collection.find({"poscar_idx":{"$exists":1}, "chemsys":"Mo-S", "start":"2_3"})

def force(outcar):
    from pymatgen.io.vasp.outputs import Outcar

    outcar = Outcar(outcar)

    forces = outcar.read_table_pattern(
        header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
        row_pattern=r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        footer_pattern=r"\s--+",
        postprocess=lambda x: float(x),
        last_one_only=False
    )
    forces = np.array(forces[0])
    f = np.max(forces)
    return f


en, dxs, fs, z = [], [], [], []
for e in es:
    dir_name = e["calcs_reversed"][0]["dir_name"]
    energy = e["output"]["energy"]
    dx = e["poscar_idx"]
    f = None
    try:
        f = force(os.path.join(dir_name, "OUTCAR.gz"))
    except:
        f = force(os.path.join(dir_name, "OUTCAR"))
    dx = dx.split("_")[-1].split(".")[0]
    dx = int(dx)*0.1
    en.append(energy)
    dxs.append(dx)
    fs.append(f)

    st = Structure.from_dict(e["output"]["structure"])
    z.append(st[25].coords[-1])

df = pd.DataFrame({"dx": dxs, "f": fs, "en":en, 'z': z})
df.sort_values("dx", inplace=True)

fig, ax = plt.subplots(3,1, sharex=True)
ax[0].scatter(df["dx"], df["z"])
ax[1].scatter(df["dx"], df["f"])
ax[2].scatter(df["dx"], df["en"])

plt.show()


