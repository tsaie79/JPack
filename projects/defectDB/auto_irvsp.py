import mpinterfaces.utils
from atomate.vasp.workflows.base.core import get_wf
from fireworks import LaunchPad, Workflow
import os
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    set_queue_options,
    set_execution_options,
    clean_up_files,
    add_modify_kpoints
)
from qubitPack.tool_box import get_db
import numpy as np
from pymatgen.core.structure import Structure, SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MPRelaxSet
from subprocess import call
from pytopomat.workflows.fireworks import IrvspFW


c2db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12345, user="adminUser", password="qiminyan").collection
for spg in c2db.distinct("spacegroup"):
    print(spg)
    e = c2db.find_one({"spacegroup": spg, "magstate":"NM"})
    try:
        st = e["structure"]
    except Exception:
        continue

    os.makedirs("symmetrized_st", exist_ok=True)
    os.chdir("symmetrized_st")
    st = Structure.from_dict(st)
    st.to("poscar", "POSCAR")
    call("phonopy --symmetry --tolerance 0.01 -c POSCAR".split(" "))
    st = Structure.from_file("PPOSCAR")
    st.to("poscar", "POSCAR")
    call("pos2aBR")
    st = Structure.from_file("POSCAR_std")
    os.chdir("..")

    wf = get_wf(st, "/home/tug03990/site-packages/pytopomat/pytopomat/workflows/irvsp_hse_sp.yaml")
    fws = wf.fws[:2]
    fw_to_db = IrvspFW(structure=st, parents=fws[1], additional_fields={"c2db_uid": e["uid"],
                                                                        "spg_c2db": e["spacegroup"],
                                                                        "spg": SpacegroupAnalyzer(st).get_space_group_symbol()
                                                                        })
    fws.append(fw_to_db)
    wf = Workflow(fws, name=wf.name)
    lpad = LaunchPad.from_file(os.path.expanduser(
        os.path.join("~", "config/project/testIR/irvsp_test/my_launchpad.yaml")))
    wf = clean_up_files(wf, ("WAVECAR*", "CHGCAR*"), wf.fws[-1].name, task_name_constraint=wf.fws[-1].tasks[-1].fw_name)
    magmom = MPRelaxSet(st).incar.get("MAGMOM", None)
    wf = add_modify_incar(wf, {"incar_update": {"LWAVE": True, "ISYM":3, "MAGMOM": magmom}}, "line")
    wf = add_modify_incar(wf, {"incar_update": {"ISPIN":2}})
    wf = set_execution_options(wf, category="irvsp_test")
    wf = preserve_fworker(wf)
    wf.name = wf.name + ":{}".format(spg)
    # lpad.add_wf(wf)

#%%
from atomate.vasp.workflows.base.core import get_wf
from fireworks import LaunchPad, Workflow
import os
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    set_queue_options,
    set_execution_options,
    clean_up_files,
    add_modify_kpoints
)
from qubitPack.tool_box import get_db
import numpy as np
from pymatgen.core.structure import Structure, SymmOp, PeriodicSite
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from subprocess import call
from pytopomat.workflows.fireworks import IrvspFW


c2db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12345, user="adminUser", password="qiminyan").collection

e = c2db.find_one({"uid": "WS2-MoS2-NM"})

os.makedirs("symmetrized_st", exist_ok=True)
os.chdir("symmetrized_st")
st = Structure.from_dict(e["structure"])
st.to("poscar", "POSCAR")
call("phonopy --symmetry --tolerance 0.01 -c POSCAR".split(" "))
st = Structure.from_file("PPOSCAR")
st.to("poscar", "POSCAR")
call("pos2aBR")
st = Structure.from_file("POSCAR_std")
os.chdir("..")

wf = get_wf(st, "/home/tug03990/site-packages/pytopomat/pytopomat/workflows/irvsp_hse_sp.yaml")
fws = wf.fws[:3]
fw_to_db = IrvspFW(structure=st, parents=fws[-1], additional_fields={"c2db_uid": e["uid"],
                                                                    "spg_c2db": e["spacegroup"],
                                                                    "spg": SpacegroupAnalyzer(st).get_space_group_symbol()
                                                                    })
fws.append(fw_to_db)
wf = Workflow(fws, name=wf.name)
lpad = LaunchPad.from_file(os.path.expanduser(
    os.path.join("~", "config/project/testIR/irvsp_test/my_launchpad.yaml")))
wf = clean_up_files(wf, ("WAVECAR*", "CHGCAR*"), wf.fws[-1].name, task_name_constraint=wf.fws[-1].tasks[-1].fw_name)
magmom = MPRelaxSet(st).incar.get("MAGMOM", None)
wf = add_modify_incar(wf, {"incar_update": {"LWAVE": True, "ISYM":3, "MAGMOM": magmom}}, "line")
wf = add_modify_incar(wf, {"incar_update": {"ISPIN":2}})
wf = set_execution_options(wf, category="irvsp_test")
wf = preserve_fworker(wf)
print(wf)
# wf.name = wf.name + ":{}".format(spg)
lpad.add_wf(wf)

#%% defect
from atomate.vasp.workflows.base.core import get_wf
from fireworks import LaunchPad, Workflow
import os
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    set_queue_options,
    set_execution_options,
    clean_up_files,
    add_modify_kpoints
)
from qubitPack.tool_box import get_db
import numpy as np
from pymatgen.core.structure import Structure, SymmOp, PeriodicSite
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from subprocess import call
from pytopomat.workflows.fireworks import IrvspFW
from pymatgen.analysis.defects.core import Vacancy, Substitution, Interstitial
from mpinterfaces.utils import ensure_vacuum
from qubitPack.tool_box import get_rand_vec, move_site

c2db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12345, user="adminUser", password="qiminyan").collection

e = c2db.find_one({"uid": "WS2-MoS2-NM"})
host = Structure.from_dict(e["structure"])
# defect_site = host.sites[1]
defect_site = PeriodicSite("W", host.sites[1].frac_coords, host.lattice)
vac = Substitution(host, defect_site, charge=0)
st = vac.generate_defect_structure((5,5,1))
move_site(st, 47, get_rand_vec(0.1))
st = ensure_vacuum(st, 20)


os.makedirs("symmetrized_st", exist_ok=True)
os.chdir("symmetrized_st")
# st = Structure.from_dict(e["structure"])
st.to("poscar", "POSCAR")
call("phonopy --symmetry --tolerance 0.01 -c POSCAR".split(" "))
st = Structure.from_file("PPOSCAR")
st.to("poscar", "POSCAR")
call("pos2aBR")
st = Structure.from_file("POSCAR_std")
os.chdir("..")


wf = get_wf(st, "/home/tug03990/site-packages/pytopomat/pytopomat/workflows/irvsp.yaml")
fws = wf.fws[:3]
fw_to_db = IrvspFW(structure=st, parents=fws[-1], additional_fields={"c2db_uid": e["uid"],
                                                                     "spg_c2db": e["spacegroup"],
                                                                     "spg": SpacegroupAnalyzer(st).get_space_group_symbol()
                                                                     })
fws.append(fw_to_db)
wf = Workflow(fws, name=wf.name)
lpad = LaunchPad.from_file(os.path.expanduser(
    os.path.join("~", "config/project/testIR/irvsp_test/my_launchpad.yaml")))
wf = clean_up_files(wf, ("WAVECAR*", "CHGCAR*"), wf.fws[-1].name, task_name_constraint=wf.fws[-1].tasks[-1].fw_name)
magmom = MPRelaxSet(st).incar.get("MAGMOM", None)
wf = add_modify_incar(wf, {"incar_update": {"LWAVE": True, "ISYM":3, "MAGMOM": magmom}}, "line")
wf = add_modify_incar(wf, {"incar_update": {"ISPIN":2}})
wf = set_execution_options(wf, category="irvsp_test")
wf = preserve_fworker(wf)
print(wf)
# wf.name = wf.name + ":{}".format(spg)
lpad.add_wf(wf)

#%%
from atomate.vasp.workflows.base.core import get_wf
from fireworks import LaunchPad, Workflow
import os
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    set_queue_options,
    set_execution_options,
    clean_up_files,
    add_modify_kpoints
)
from qubitPack.tool_box import get_db
import numpy as np
from pymatgen.core.structure import Structure, SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MPRelaxSet
from subprocess import call
from pytopomat.workflows.fireworks import IrvspFW
import pandas as pd

c2db = get_db("ML_data", "PBE_bulk", port=12345)
problems = pd.read_csv("/home/tug03990/scripts/JPack/projects/defectDB/ids.txt", header=None)[0]
for spg in problems:
    print(spg)
    e = c2db.collection.find_one({"task_id": spg})
    try:
        st = e["output"]["structure"]
    except Exception:
        continue

    os.makedirs("symmetrized_st", exist_ok=True)
    os.chdir("symmetrized_st")
    st = Structure.from_dict(st)
    st.to("poscar", "POSCAR")
    call("phonopy --symmetry --tolerance 0.01 -c POSCAR".split(" "))
    st = Structure.from_file("PPOSCAR")
    st.to("poscar", "POSCAR")
    call("pos2aBR")
    st = Structure.from_file("POSCAR_std")
    os.chdir("..")

    wf = get_wf(st, "/home/tug03990/site-packages/pytopomat/pytopomat/workflows/irvsp_hse_sp.yaml")
    fws = wf.fws[:3]
    fw_to_db = IrvspFW(structure=st, parents=fws[-1], additional_fields={"PBE_bulk_tkid": e["task_id"],
                                                                        # "spg_c2db": e["spacegroup"],
                                                                        "spg": SpacegroupAnalyzer(st).get_space_group_symbol()
                                                                        })
    fws.append(fw_to_db)
    wf = Workflow(fws, name=wf.name)
    lpad = LaunchPad.from_file(os.path.expanduser(
        os.path.join("~", "config/project/testIR/irvsp_test/my_launchpad.yaml")))
    wf = clean_up_files(wf, ("WAVECAR*", "CHGCAR*"), wf.fws[-1].name, task_name_constraint=wf.fws[-1].tasks[-1].fw_name)
    magmom = MPRelaxSet(st).incar.get("MAGMOM", None)
    wf = add_modify_incar(wf, {"incar_update": {"LWAVE": True, "ISYM":3, "MAGMOM": magmom}}, "line")
    wf = add_modify_incar(wf, {"incar_update": {"ISPIN":2}})
    wf = set_execution_options(wf, category="irvsp_test")
    wf = preserve_fworker(wf)
    wf.name = wf.name + ":{}".format(spg)
    print(wf)
    lpad.add_wf(wf)
    break

#%% no line test
from atomate.vasp.workflows.base.core import get_wf
from fireworks import LaunchPad, Workflow
import os
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    set_queue_options,
    set_execution_options,
    clean_up_files,
    add_modify_kpoints
)
from qubitPack.tool_box import get_db
import numpy as np
from pymatgen.core.structure import Structure, SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MPRelaxSet
from subprocess import call
from pytopomat.workflows.fireworks import IrvspFW
import pandas as pd
from mpinterfaces.utils import ensure_vacuum

c2db = get_db("ML_data", "PBE_bulk", port=12345)
problems = pd.read_csv("/home/tug03990/scripts/JPack/projects/defectDB/ids.txt", header=None)[0]
for spg in problems:
    print(spg)
    e = c2db.collection.find_one({"task_id": spg})
    try:
        st = e["output"]["structure"]
    except Exception:
        continue

    os.makedirs("symmetrized_st", exist_ok=True)
    os.chdir("symmetrized_st")
    st = Structure.from_dict(st)
    st = ensure_vacuum(st, 15)
    st.to("poscar", "POSCAR")
    call("phonopy --symmetry --tolerance 0.01 -c POSCAR".split(" "))
    st = Structure.from_file("PPOSCAR")
    st.to("poscar", "POSCAR")
    call("pos2aBR")
    st = Structure.from_file("POSCAR_std")
    os.chdir("..")

    wf = get_wf(st, "/home/tug03990/site-packages/pytopomat/pytopomat/workflows/irvsp_hse_sp.yaml")

    lpad = LaunchPad.from_file(os.path.expanduser(os.path.join("~", "config/project/testIR/irvsp_test/my_launchpad.yaml")))
    wf = clean_up_files(wf, ("WAVECAR*", "CHGCAR*"), wf.fws[-1].name, task_name_constraint=wf.fws[-1].tasks[-1].fw_name)
    # magmom = MPRelaxSet(st).incar.get("MAGMOM", None)
    # wf = add_modify_incar(wf, {"incar_update": {"MAGMOM": magmom}}, "line")
    uis_encut = MPRelaxSet(st).incar.get("ENCUT", None)*1.3
    wf = add_modify_incar(wf, {"incar_update": {"ENCUT": uis_encut}})
    wf = set_execution_options(wf, category="irvsp_test")
    wf = preserve_fworker(wf)
    wf.name = wf.name + ":{}".format(spg)
    print(wf)
    # lpad.add_wf(wf)
    break