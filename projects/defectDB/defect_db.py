#%% IMPORT
from fireworks import LaunchPad, Workflow

from atomate.vasp.fireworks.core import ScanOptimizeFW
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    set_queue_options,
    set_execution_options,
    clean_up_files,
)

from my_atomate_jyt.vasp.powerups import *
from my_atomate_jyt.vasp.workflows.wf_full import get_wf_full_scan

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet


from qubitPack.tool_box import *

from mpinterfaces.utils import *

from monty.serialization import loadfn
import os
import pandas as pd

#%% looking for candidate hosts
from phonopy.phonon.irreps import character_table

db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/code/JPack/qubitPack/database_profile/db_dk_local.json")
col = db.collection

filter = {
    "gap_hse_nosoc":{"$gte":1},
    # "ehull":{"$lt":0.4},
    "nkinds":2,
    "magstate": "NM",
    "uid": "Ga2Cl6-TiCl3-NM"
}
show = {
    "formula":1,
    "class":1,
    "ehull":1,
    "spacegroup":1,
    "gap_hse_nosoc":1
}

data = []
for e in col.find(filter):
    st = Structure.from_dict(e["structure"])
    space_sym_analyzer = SpacegroupAnalyzer(st, symprec=1e-4)
    sites = space_sym_analyzer.get_symmetry_dataset()["std_positions"]
    site_syms = space_sym_analyzer.get_symmetry_dataset()["site_symmetry_symbols"]
    point_gp = space_sym_analyzer.get_symmetry_dataset()["pointgroup"]
    vw = space_sym_analyzer.get_symmetry_dataset()["wyckoffs"]
    site_irreps = []
    for specie, site, site_sym in zip(st.species, sites, site_syms):
        print(site_sym)
        site = list(site.round(4))
        site_sym = [x for x in site_sym.split(".") if x][0]
        if site_sym == "-4m2":
            site_sym = "-42m"
        if site_sym == "2mm" or site_sym == "m2m":
            site_sym = "mm2"
        if site_sym == "1":
            continue
        irreps = character_table[site_sym][0]["character_table"]
        for irrep, char_vec in irreps.items():
            print(char_vec)
            if char_vec[0] >= 2:
                site_irreps.append((str(specie), site, irrep, char_vec[0]))
                # site_irreps.append(str(specie))
                break




    print(e["uid"])
    print(site_irreps)
    data.append(
        {
            "formula": e["formula"],
            "species": [str(sp) for sp in st.species],
            "site_sym": site_syms,
            # "irreps": site_irreps,
            "Wf_position": vw,
            "spacegroup": e["spacegroup"],
            "pmg_point_gp": point_gp,
            "gap_hse": e["gap_hse"],
            "gap_hse_nosoc": e["gap_hse_nosoc"],
            "gap_nosoc": e["gap_nosoc"],
            "soc_intensity": e["gap_hse_nosoc"] - e["gap_hse"],
            "ehull": e["ehull"],
            "uid":e["uid"],
            # "structure": e["structure"]
        }
    )



df = pd.DataFrame(data)
df = df.round(3).sort_values(["soc_intensity", "ehull"], ascending=True)
df.set_index("uid", inplace=True)
# df.to_json("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/gap_gt1-binary-NM-full.json",
#            orient="index", indent=4, index=True)
# df.to_excel("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/gap_gt1-binary-NM-full.xlsx", index=False)

#%%
db = VaspCalcDb.from_db_file("/global/homes/t/tsaie79/code/qubitPack/database_profile/db_mongo_cori.json")
col = db.collection
data = []
for e in col.find():
    data.append(
        {
            "formula": e["formula_pretty"],
            "spacegroup": e["output"]["spacegroup"]["symbol"],
            "bandgap": e["output"]["bandgap"],
            "task_id": e["task_id"]
        }
    )

df = pd.DataFrame(data).sort_values(["bandgap"], ascending=False)
df.to_json("cori_relax_lc.json", orient="records")


#%% get structure from SCAN_relax

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'scan_relax_pc/db.json')

# task_id = 34
filter = {"nsites":4}
es = db.collection.find(filter)[:5]
for e in es:
    st = Structure.from_dict(e["output"]["structure"])
    st.to("POSCAR", '/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                    'scan_relax_pc/structures/tid_{}_ehull_{}.vasp'.format(e["task_id"], e["c2db_info"]["ehull"]))

#%% look for defect states
from qubitPack.qc_searching.analysis.main import get_defect_state
from qubitPack.tool_box import get_db

defect_db = get_db("defect_db", "binary_defect")
host_db = get_db("defect_db", "scan_bs_pc")
c2db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan")


tk_id = 14
# pc_from_id = defect_db.collection.find_one({"task_id": tk_id})["pc_from_id"]
# c2db_uid = host_db.collection.find_one({"task_id": pc_from_id})["c2db_info"]["uid"]

tot, proj, d_df = get_defect_state(
    defect_db,
    {"task_id": tk_id},
    1, -5,
    None,
    True,
    "proj",
    (host_db, 78, 0), #(get_db("antisiteQubit", "W_S_Ef"), 312, 0.),
    0.1,
    locpot_c2db=None #(c2db, c2db_uid, 0)
)

#%% nsites vs
