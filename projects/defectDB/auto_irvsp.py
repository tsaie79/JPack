from atomate.vasp.workflows.base.core import get_wf
from pymatgen import Structure
from fireworks import LaunchPad
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
from phonopy.structure.symmetry import ge



st = {'@module': 'pymatgen.core.structure',
      '@class': 'Structure',
      'charge': None,
      'lattice': {'matrix': [[5.468728, 0.0, 3.348630120303753e-16],
                             [-3.348630120303753e-16, 5.468728, 3.348630120303753e-16],
                             [0.0, 0.0, 5.468728]],
                  'a': 5.468728,
                  'b': 5.468728,
                  'c': 5.468728,
                  'alpha': 90.0,
                  'beta': 90.0,
                  'gamma': 90.0,
                  'volume': 163.55317139465933},
      'sites': [{'species': [{'element': 'Si', 'occu': 1.0}],
                 'abc': [0.25, 0.75, 0.25],
                 'xyz': [1.3671819999999997, 4.101546, 1.3671820000000003],
                 'label': 'Si',
                 'properties': {}},
                {'species': [{'element': 'Si', 'occu': 1.0}],
                 'abc': [0.0, 0.0, 0.5],
                 'xyz': [0.0, 0.0, 2.734364],
                 'label': 'Si',
                 'properties': {}},
                {'species': [{'element': 'Si', 'occu': 1.0}],
                 'abc': [0.25, 0.25, 0.75],
                 'xyz': [1.367182, 1.367182, 4.101546],
                 'label': 'Si',
                 'properties': {}},
                {'species': [{'element': 'Si', 'occu': 1.0}],
                 'abc': [0.0, 0.5, 0.0],
                 'xyz': [-1.6743150601518765e-16, 2.734364, 1.6743150601518765e-16],
                 'label': 'Si',
                 'properties': {}},
                {'species': [{'element': 'Si', 'occu': 1.0}],
                 'abc': [0.75, 0.75, 0.75],
                 'xyz': [4.101546, 4.101546, 4.101546000000001],
                 'label': 'Si',
                 'properties': {}},
                {'species': [{'element': 'Si', 'occu': 1.0}],
                 'abc': [0.5, 0.0, 0.0],
                 'xyz': [2.734364, 0.0, 1.6743150601518765e-16],
                 'label': 'Si',
                 'properties': {}},
                {'species': [{'element': 'Si', 'occu': 1.0}],
                 'abc': [0.75, 0.25, 0.25],
                 'xyz': [4.101546, 1.367182, 1.3671820000000003],
                 'label': 'Si',
                 'properties': {}},
                {'species': [{'element': 'Si', 'occu': 1.0}],
                 'abc': [0.5, 0.5, 0.5],
                 'xyz': [2.734364, 2.734364, 2.7343640000000002],
                 'label': 'Si',
                 'properties': {}}]}


st = Structure.from_dict(st)

c2db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12345, user="adminUser", password="qiminyan").collection
for spg in c2db.distinct("spacegroup"):
    print(spg)
    e = c2db.find_one({"spacegroup": spg, "magstate":"NM"})
    try:
        st = e["structure"]
    except Exception:
        continue
    st = Structure.from_dict(st)


    wf = get_wf(st, "/home/tug03990/site-packages/pytopomat/pytopomat/workflows/irvsp_hse.yaml")
    lpad = LaunchPad.from_file(os.path.expanduser(
        os.path.join("~", "config/project/testIR/irvsp_test/my_launchpad.yaml")))

    wf = add_modify_incar(wf, {"incar_update": {"LWAVE": True, "ISYM":3, "ISPIN":1}}, "line")
    wf = add_modify_incar(wf, {"incar_update": {"ISPIN":1}})
    wf = add_additional_fields_to_taskdocs(wf, {"c2db_uid": e["uid"], "spg": e["spacegroup"]})
    wf = set_execution_options(wf, category="irvsp_test")
    wf = preserve_fworker(wf)
    wf.name = wf.name + ":{}".format(spg)
    lpad.add_wf(wf)

#%%
from atomate.vasp.workflows.base.core import get_wf
from pymatgen import Structure
from fireworks import LaunchPad
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
from subprocess import call

st = Structure.from_dict(st)

c2db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12345, user="adminUser", password="qiminyan").collection

e = c2db.find_one({"uid":"MoS2-MoS2-NM"})

st = Structure.from_dict(e["structure"])
st.to("poscar", "POSCAR")
call("phonopy --symmetry --tolerance 0.01 -c POSCAR".split(" "))
st = Structure.from_file("POSCAR")

wf = get_wf(st, "/home/tug03990/site-packages/pytopomat/pytopomat/workflows/irvsp_hse.yaml")
lpad = LaunchPad.from_file(os.path.expanduser(
    os.path.join("~", "config/project/testIR/irvsp_test/my_launchpad.yaml")))

wf = add_modify_incar(wf, {"incar_update": {"LWAVE": True, "ISYM":3, "ISPIN":1}}, "line")
wf = add_modify_incar(wf, {"incar_update": {"ISPIN":1}})
wf = add_additional_fields_to_taskdocs(wf, {"c2db_uid": e["uid"], "spg": e["spacegroup"]})
wf = set_execution_options(wf, category="irvsp_test")
wf = preserve_fworker(wf)
wf.name = wf.name + ":{}".format(spg)
lpad.add_wf(wf)