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
    st = Structure.from_dict(e["output"]["structure"])
    species, site_syms, wyckoffs, pg, site_idx = get_unique_sites_from_wy(st, symprec=1e-2).values()
    good_ir = get_good_ir_sites(st, symprec=1e-2)

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
    d.update({"symprec":1e-2, "pmg_point_gp": pg})

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
