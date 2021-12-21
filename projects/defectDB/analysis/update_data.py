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

db = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=12347)
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

src = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=12347)
tgt = get_db("Scan2dDefect", "calc_data", user="Jeng", password="qimin", port=12347)

col = tgt.collection
# filter = {"task_label": "SCAN_scf",}
# filter = {"task_id": {"$in":[1309, 1311, 1315, 1349, 1351]}}
filter = {"task_label": {"$regex": "SCAN_scf"}}#, "host_info.scan_bs.band_edges.vbm_is_up_dn_band_idx_equal": {
    #"$exists": 0}}
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
            "host_info.sym_data": sym_data,
            "host_info.scan_bs": scan_bs_output
        }})
    except Exception as er:
        print(er)

#%% 5
from qubitPack.tool_box import get_db, get_band_edges_characters
from monty.json import jsanitize

tgt = get_db("Scan2dDefect", "calc_data", user="Jeng", password="qimin", port=12347)

col = tgt.collection
filter = {"task_label": "SCAN_scf"}

for e in list(col.find(filter)):
    print(e["task_id"], e["host_info"]["c2db_info"]["uid"])
    site_syms = e["host_info"]["sym_data"]["unique_wyckoff"]["site_sym"]
    if site_syms[0] == site_syms[1]:
        site_symmetry_uniform = True
    else:
        site_symmetry_uniform = False

    tgt.collection.update_one({"task_id": e["task_id"]}, {"$set":{"site_symmetry_uniform": site_symmetry_uniform}})




#%% 6
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

#%% 7 update SCAN_scf with eigenstate info
import json, os
from monty.json import MontyEncoder, MontyDecoder
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.drones import VaspDrone

local_scf_path = "/home/qimin/sdc_tsai/Research/projects/Scan2dDefect/calc_data/scf"

# run this code in db1
mmdb = VaspCalcDb.from_db_file('/mnt/sdb/tsai/scripts/update_eigen/db.json', admin=True)
tgt = list(mmdb.collection.find({"task_label": "SCAN_scf"}))

for e in tgt:
    t_id = e["task_id"]
    print("++"*20 + str(t_id))
    if e["calcs_reversed"][0].get("eigenvalues_fs_id", None):
        continue
    #    if t_id != 3257:
    #        continue

    path = e["dir_name"].split("/")[-1]
    drone = VaspDrone(parse_eigenvalues=True, parse_dos=True)
    task_doc = drone.assimilate(os.path.join(local_scf_path, path))

    eigenvals = {}
    dos = None
    if "calcs_reversed" in task_doc:
        for eigenvalue in ("eigenvalues", "projected_eigenvalues"):
            if eigenvalue in task_doc["calcs_reversed"][0]["output"]:  # only store idx=0 data
                eigenvals[eigenvalue] = json.dumps(task_doc["calcs_reversed"][0]["output"][eigenvalue],
                                                   cls=MontyEncoder)
                del task_doc["calcs_reversed"][0]["output"][eigenvalue]

            if "dos" in task_doc["calcs_reversed"][0]:  # only store idx=0 (last step)
                dos = json.dumps(task_doc["calcs_reversed"][0]["dos"], cls=MontyEncoder)
                del task_doc["calcs_reversed"][0]["dos"]

    if eigenvals:
        for name, data in eigenvals.items():
            data_gfs_id, compression_type = mmdb.insert_gridfs(data, "{}_fs".format(name), task_id=t_id)
            mmdb.collection.update_one(
                {"task_id": t_id}, {"$set": {"calcs_reversed.0.{}_compression".format(name): compression_type}})
            mmdb.collection.update_one({"task_id": t_id}, {"$set": {"calcs_reversed.0.{}_fs_id".format(name): data_gfs_id}})
    if dos:
        dos_gfs_id, compression_type = mmdb.insert_gridfs(
            dos, "dos_fs", task_id=t_id
        )
        mmdb.collection.update_one(
            {"task_id": t_id},
            {"$set": {"calcs_reversed.0.dos_compression": compression_type}},
        )
        mmdb.collection.update_one(
            {"task_id": t_id}, {"$set": {"calcs_reversed.0.dos_fs_id": dos_gfs_id}}
        )

#%% 7*
from monty.json import jsanitize
from qubitPack.tool_box import get_db
from pymatgen import Structure

db = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=12347)
col = db.collection
filter = {"task_label": "SCAN_nscf line"}

for e in list(col.find(filter))[:]:
    print(e["task_id"], e["c2db_info"]["uid"])
    structure = Structure.from_dict(e["output"]["structure"]).copy()
    structure.add_oxidation_state_by_guess()
    site_oxi_state = []
    for site in structure.sites:
        specie = site.specie.as_dict()
        site_oxi_state.append((specie["element"], specie["oxidation_state"]))
    data = {"site_oxi_state": site_oxi_state}
    d = jsanitize(data)
    db.collection.update_one({"task_id": e["task_id"]}, {"$set": data})

#%% 8
from monty.json import jsanitize
from qubitPack.tool_box import get_db
from pymatgen import Structure

host = get_db("Scan2dMat", "calc_data", user="Jeng", password="qimin", port=12347)
db = get_db("Scan2dDefect", "calc_data", user="Jeng", password="qimin", port=12347)

col = db.collection
filter = {"task_label": "SCAN_scf"}

for e in list(col.find(filter))[:]:
    print(e["task_id"], e["host_info"]["c2db_info"]["uid"])
    host_entry = host.collection.find_one({"c2db_info.uid": e["host_info"]["c2db_info"]["uid"], "task_label":
        "SCAN_nscf line"})
    site_oxi_state = host_entry["site_oxi_state"]
    db.collection.update_one({"task_id": e["task_id"]}, {"$set": {"host_info.scan_bs.site_oxi_state": site_oxi_state}})
