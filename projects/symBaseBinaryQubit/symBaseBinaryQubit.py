from atomate.vasp.database import VaspCalcDb

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen import Structure

import pandas as pd
#%%
col = VaspCalcDb.from_db_file(
    "/Users/jeng-yuantsai/Research/project/qubit/calculations/symBaseBinaryQubit/BN_vac/db.json").collection

data = []
for e in col.find({"task_label":{"$in": ["SCAN_scf", "SCAN_relax"]}, "input.incar.SIGMA":{"$in":[0.01, 0.005]}}):
    entry = {}
    # entry["kpoints"] = e["calcs_reversed"][0]["input"]["kpoints"]["kpoints"][0]
    # entry["dir"] = e["calcs_reversed"][0]["dir_name"]
    # entry["complete"] = e["calcs_reversed"][0]["completed_at"].split(" ")[0]
    # entry["evac"] = max(e["calcs_reversed"][0]["output"]["locpot"]["2"])
    # entry["vbm"] = e["calcs_reversed"][0]["output"]["vbm"]
    # entry["cbm"] = e["calcs_reversed"][0]["output"]["cbm"]
    # entry["Band Gap (eV)"] = "{:.3f}".format(e["calcs_reversed"][0]["output"]["bandgap"])
    entry["task_id"] = e["task_id"]
    entry["task_label"] = e["task_label"]
    entry["chemsys"] = e["chemsys"]
    entry["defect"] = e["defect_entry"]["name"]
    entry["formula"] = e["formula_pretty"]
    entry["energy"] = e["output"]["energy"]
    # entry["LVHAR"] = bool(e["input"]["parameters"]["LVHAR"])
    # entry["ENCUT"] = e["input"]["incar"]["ENCUT"]
    # entry["nsites"] = e["nsites"]
    # entry["Lz"] = e["input"]["structure"]["lattice"]["c"]
    entry["charge_state"] = MPRelaxSet(Structure.from_dict(e["input"]["structure"])).nelect - e["input"]["incar"]["NELECT"]
    entry["magmom"] = e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]
    entry["point_group"] = e["output"]["spacegroup"]["point_group"]
    data.append(entry)

df = pd.DataFrame(data)
df = df.sort_values(["formula", "charge_state"])
df.to_clipboard()

#%%
col = VaspCalcDb.from_db_file(
    "/Users/jeng-yuantsai/Research/project/qubit/calculations/symBaseBinaryQubit/TMDCs_complex/db.json").collection

data = []
for e in col.find({"task_label":{"$in":["SCAN_scf", "SCAN_relax"]}, "input.incar.SIGMA":{"$in":[0.01, 0.005]}}):
    entry = {}
    # entry["kpoints"] = e["calcs_reversed"][0]["input"]["kpoints"]["kpoints"][0]
    # entry["dir"] = e["calcs_reversed"][0]["dir_name"]
    # entry["complete"] = e["calcs_reversed"][0]["completed_at"].split(" ")[0]
    # entry["evac"] = max(e["calcs_reversed"][0]["output"]["locpot"]["2"])
    # entry["vbm"] = e["calcs_reversed"][0]["output"]["vbm"]
    # entry["cbm"] = e["calcs_reversed"][0]["output"]["cbm"]
    # entry["Band Gap (eV)"] = "{:.3f}".format(e["calcs_reversed"][0]["output"]["bandgap"])
    entry["task_id"] = e["task_id"]
    entry["task_label"] = e["task_label"]
    entry["chemsys"] = e["chemsys"]
    entry["defect"] = e["defect_entry"]["name"]
    entry["formula"] = e["formula_pretty"]
    entry["energy"] = e["output"]["energy"]
    # entry["LVHAR"] = bool(e["input"]["parameters"]["LVHAR"])
    # entry["ENCUT"] = e["input"]["incar"]["ENCUT"]
    # entry["nsites"] = e["nsites"]
    # entry["Lz"] = e["input"]["structure"]["lattice"]["c"]
    entry["charge_state"] = MPRelaxSet(Structure.from_dict(e["input"]["structure"])).nelect - e["input"]["incar"]["NELECT"]
    entry["magmom"] = e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]
    entry["point_group"] = e["output"]["spacegroup"]["point_group"]
    data.append(entry)

df = pd.DataFrame(data)
df = df.sort_values(["formula", "charge_state"])
