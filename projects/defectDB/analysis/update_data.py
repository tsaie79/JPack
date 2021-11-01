#%% 1
from qubitPack.tool_box import get_db, get_good_ir_sites
from monty.json import jsanitize
from pymatgen.io.vasp.inputs import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

db = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=1236)

for e in db.collection.find():
    st = Structure.from_dict(e["output"]["structure"])
    sym = SpacegroupAnalyzer(st, symprec=1e-1).get_symmetry_dataset()
    d = jsanitize(sym)
    print(d)
    d.update({"symprec":1e-1})
    species = [str(sp) for sp in st.species]
    combine = dict(zip(d["wyckoffs"], zip(species, d["site_symmetry_symbols"])))
    species, syms = get_good_ir_sites(dict(combine.values()).keys(), dict(combine.values()).values())
    combine.update({"species": species, "syms": syms})

    d.update({"high_dim_ir_info": combine})
    print(d)
    db.collection.update_one({"task_id": e["task_id"]}, {"$set":{"sym_data":d}})
#%% 2
from qubitPack.tool_box import get_db, get_good_ir_sites, get_unique_sites_from_wy
from monty.json import jsanitize
from pymatgen.io.vasp.inputs import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

db = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=1236)

for e in db.collection.find({"task_label": "SCAN_scf"}):
    symprec = 1e-1
    st = Structure.from_dict(e["output"]["structure"])
    species, site_syms, wyckoffs, pg, site_idx, spg, spg_number = get_unique_sites_from_wy(st, symprec=symprec).values()
    good_ir = get_good_ir_sites(st, symprec=symprec)

    print(e["c2db_info"]["uid"])
    d = {}
    wy_fields = {
        "species": species,
        "site_sym": site_syms,
        "wyckoffs": wyckoffs,
        "site_idx": site_idx,
    }

    if good_ir["site_syms"]:
        ir = {
            "species": good_ir["species"],
            "site_sym": good_ir["site_syms"],
            "wyckoffs": good_ir["wys"],
            "site_idx": good_ir["site_idx"],
        }
    else:
        ir = {
            "species": [],
            "site_sym": [],
            "wyckoffs": [],
            "site_idx": [],
        }

    d["unique_wyckoff"] = wy_fields
    d["good_ir_info"] = ir

    d = jsanitize(d)
    d.update({"symprec":symprec, "pmg_pg": pg, "pmg_spg": spg, "pmg_spg_number": spg_number})

    db.collection.update_one({"task_id": e["task_id"]}, {"$set":{"sym_data":d}})
#%% 3
from qubitPack.tool_box import get_db, get_band_edges_characters
from monty.json import jsanitize

db = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=1236)
col = db.collection
filter = {"task_label": "SCAN_nscf line"}

for e in list(col.find(filter)):
    print(e["task_id"], e["c2db_info"]["uid"])
    bs = db.get_band_structure(e["task_id"])
    if bs.get_band_gap()["energy"] == 0:
        print("skip")
        continue
    data = get_band_edges_characters(bs)
    d = jsanitize(data)
    db.collection.update_one({"task_id": e["task_id"]}, {"$set":{"band_edges":d}})
#%% 4
from qubitPack.tool_box import get_db, get_band_edges_characters
from monty.json import jsanitize

src = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=1236)
tgt = get_db("Scan2dDefect", "calc_data", user="Jeng", password="qimin", port=1236)

col = tgt.collection
# filter = {"task_label": "SCAN_scf",}
# filter = {"task_id": {"$in":[1309, 1311, 1315, 1349, 1351]}}
filter = {"task_label": {"$regex": "uniform"}}
for e in list(col.find(filter)):
    try:
        print(e["task_id"], e["host_info"]["c2db_info"]["uid"])
        src_entry = src.collection.find_one({"task_id": e["pc_from_id"]})
        locpot = src_entry["calcs_reversed"][0]["output"]["locpot"]["2"]
        vacuum = max(locpot)
        sym_data = src_entry["sym_data"]
        scan_bs_output = src_entry["output"]
        for remove in ["structure", "density", "energy", "energy_per_atom", "forces", "stress", "spacegroup"]:
            scan_bs_output.pop(remove)
        band_edges = src_entry["band_edges"]
        scan_bs_output.update({"band_edges": band_edges, "vacuum_level": vacuum})
        tgt.collection.update_one({"task_id": e["task_id"]}, {"$set":{
            "host_info.sym_data":sym_data,
            "host_info.scan_bs": scan_bs_output
        }})
    except Exception as er:
        print(er)

#%% 5
from qubitPack.tool_box import get_db, get_band_edges_characters
from monty.json import jsanitize

tgt = get_db("Scan2dDefect", "calc_data", user="Jeng", password="qimin", port=1236)

col = tgt.collection
filter = {"group_id": 4}

for e in list(col.find(filter)):
    try:
        print(e["task_id"], e["host_info"]["c2db_info"]["uid"])
        tgt.collection.update_one({"task_id": e["task_id"]}, {"$set":{
            "site_symmetry_uniform": True,
        }})
    except Exception as er:
        print(er)


#%%5
import pandas as pd
import os, json
from qubitPack.tool_box import remove_entry_in_db

tgt = get_db("Scan2dDefect", "ir_data", user="Jeng", password="qimin", port=1236)

gp_id = [0,
         3,
         5,
         7,
         9,
         10,
         12,
         13,
         14,
         15,
         16,
         17,
         18,
         19,
         20]
for  i in tgt.collection.find({"group_id": {"$in": gp_id}}):
        tk_id = i["task_id"]
        remove_entry_in_db(tk_id, tgt, pmg_file=False)


#%%6
from qubitPack.tool_box import get_db, get_band_edges_characters
from monty.json import jsanitize

src = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="readUser", password="qiminyan", port=12345)
tgt = get_db("Scan2dDefect", "calc_data", user="Jeng", password="qimin", port=12347)

col = tgt.collection
filter = {}

for e in col.find({"host_info.c2db_info.prototype": {"$exists": 0}}):
    try:
        print(e["task_id"], e["host_info"]["c2db_info"]["uid"])
        src_entry = src.collection.find_one({"uid": e["host_info"]["c2db_info"]["uid"]})
        tgt.collection.update_one({"task_id": e["task_id"]}, {"$set":{
            "host_info.c2db_info.prototype": src_entry["prototype"],
        }})
    except Exception as er:
        print(er)

