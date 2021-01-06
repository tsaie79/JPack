#%% search_triplet_from_defect_db
from atomate.vasp.database import VaspCalcDb

from pymatgen import Element

import pandas as pd

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'search_triplet_from_defect_db/db.json')

pc_db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc/db.json")

d = []

filter = {
    "$or":[{"calcs_reversed.output.outcar.total_magnetization":{"$gte":1.5, "$lte":2.1}},
           {"calcs_reversed.output.outcar.total_magnetization":{"$lte":-1.5, "$gte":-2.1}}],
    # "nsites":75,
    "task_label":"SCAN_scf"
}

for e in db.collection.find(filter):
    entry = {}
    entry["task_id"] = e["task_id"]
    entry["formula_pretty"] = e["formula_pretty"]
    entry["mag"] = abs(e["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"])
    entry["sym_data"] = pc_db.collection.find_one({"task_id":int(e["pc_from"].split("/")[-1])})["sym_data"]["site_symmetry_symbols"]
    entry["species"] = pc_db.collection.find_one({"task_id":int(e["pc_from"].split("/")[-1])})["sym_data"]["species"]
    species = list(dict.fromkeys(entry["species"]))
    entry["elec_configs"] = (Element(species[0]).electronic_structure, Element(species[1]).electronic_structure)
    entry["group"] = (Element(species[0]).group, Element(species[1]).group)
    entry["NN"] = e["NN"]
    entry["bandgap"] = pc_db.collection.find_one({"task_id":int(e["pc_from"].split("/")[-1])})["output"]["bandgap"]
    entry["pc_from"] = e["pc_from"]
    entry["a"] = e["input"]["structure"]["lattice"]["a"]
    entry["b"] = e["input"]["structure"]["lattice"]["b"]



    d.append(entry)

df = pd.DataFrame(d).sort_values(["sym_data"])# "group", "species"]) #"mag", "pc_from", "bandgap"])
df = df.set_index("task_id")
# df.to_clipboard()

#%% search for triplet and create sheets
from atomate.vasp.database import VaspCalcDb
from pymatgen import Element
import pandas as pd
import os
from collections import Counter
import numpy as np
from qubitPack.qc_searching.analysis.main import get_defect_state

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'search_triplet_from_defect_db/db.json')

proj_path = '/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/search_triplet_from_defect_db'
db_json = os.path.join(proj_path, "db.json")

host_path = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc"
db_host_json = os.path.join(host_path, "db.json")

es = db.collection.aggregate(
    [
        {"$match": {"task_label":"SCAN_scf"}},
        {"$addFields":{"mag": {"$abs": {"$arrayElemAt": ["$calcs_reversed.output.outcar.total_magnetization",0]}}}},
        {"$addFields":{"defect_name":"$defect_entry.name", "point_group":"$output.spacegroup.point_group"}},
        {"$project": {"_id":0, "formula_pretty":1, "nsites":1, "elements":1,
                      "mag":1, "charge_state":1, "chemsys":1, "defect_name":1,
                      "point_group":1, "task_id":1
                      }},
        {"$match":{"mag":{"$gte":1.9, "$lte":2.1}}},
        {"$sort":{"mag":-1}},
        # {"$group": {"_id":"$chemsys", "q": {"$push":"$charge_state"},"m": {"$push":"$formula_pretty"}, "nsites":{"$push":"$nsites"}}}
    ]
)

ess = []
for e in list(es):
    els = [Element(e["elements"][0]).group, Element(e["elements"][1]).group]
    els.sort()
    try:
        tot, proj, d_df = get_defect_state(
            db_json,
            {"task_id": e["task_id"]},
            1.5,-0.581,
            None,
            False,
            "dist",
            db_host_json,
            0.4 # this can affect the value of from_vbm
        )
        if d_df.loc["triplet_from"][0] == "up" and d_df.loc["up_from_vbm"][0][-1] >= 0.4 and \
                d_df.loc["up_tran_en"][0] < 0.5 and d_df.loc["up_tran_en"][0] >= 0.2:
            e.update(d_df.to_dict()[0])
            e.update({"group": "{}_{}".format(els[0], els[1])})
        elif d_df.loc["triplet_from"][0] == "dn" and d_df.loc["dn_from_vbm"][0][-1] >= 0.4 and \
                d_df.loc["dn_tran_en"][0] < 0.5 and d_df.loc["dn_tran_en"][0] >= 0.2:
            e.update(d_df.to_dict()[0])
            e.update({"group": "{}_{}".format(els[0], els[1])})
        else:
            continue
        ess.append(e)

    except Exception as er:
        # print(er)
        # e.update(
        #     {
        #         "up_from_vbm": None,
        #         "up_occ": None,
        #         "dn_from_vbm": None,
        #         "dn_occ": None
        #     }
        # )
        continue
    # e.pop("elements")

df = pd.DataFrame(ess).sort_values(["group", "point_group", "charge_state"], inplace=False)
df.to_clipboard("/t")





#%% find structure
from atomate.vasp.database import VaspCalcDb
import os
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np

from matplotlib import pyplot as plt

tid = 54
host_path = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc"
db_host_json = os.path.join(host_path, "db.json")

host = VaspCalcDb.from_db_file(db_host_json)
e = host.collection.find_one({"task_id":tid})
st = Structure.from_dict(e["input"]["structure"])
formula = e["formula_pretty"]
sym = e["c2db_info"]["spacegroup"]
sym_data = SpacegroupAnalyzer(st).get_symmetry_dataset()

st.to("poscar", os.path.join(host_path, "structures", "tid_{}_{}.vasp".format(tid, formula)))

#%% print information for defect triplet in db
from atomate.vasp.database import VaspCalcDb
from pymatgen import Element
import pandas as pd
import os
from collections import Counter
import numpy as np
from qubitPack.qc_searching.analysis.read_eigen import DetermineDefectState

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'search_triplet_from_defect_db/db.json')


es = db.collection.aggregate(
    [
        {"$match": {"task_label":"SCAN_scf"}},
        {"$group": {"_id":"$pc_from",
                    "tid": {"$push": "$task_id"},
                    "defect": {"$push": "$defect_entry.name"},
                    "point_group": {"$push": "$output.spacegroup.point_group"},
                    "mag": {"$push": "$calcs_reversed.output.outcar.total_magnetization"}
                    }
         },
    ]
)

df = pd.DataFrame(es)
df = df.transpose()



#%% plot dos
from pymatgen import Orbital
from pymatgen.electronic_structure.plotter import DosPlotter
import numpy as no

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/'
                             'search_triplet_from_defect_db/db.json')

comp_dos = db.get_dos(303)

projection_d = dict(zip(["s", "xy", "yz", "z2", "xz", "x2-y2"],
                      [comp_dos.get_site_orbital_dos(
                          comp_dos.structure.sites[25], Orbital(j)) for j in [0]+ list(range(4, 9))]))

projection_sp = dict(zip(["s", "y", "z", "x"],
                      [comp_dos.get_site_orbital_dos(
                          comp_dos.structure.sites[28], Orbital(j)) for j in range(0, 5)]))
print(comp_dos.structure.sites[27])
dos_plotter = DosPlotter(stack=False)
dos_plotter.add_dos_dict(projection_d)
plt = dos_plotter.get_plot([-8, 8])
plt.show()

#%% look defect states in detail
from qubitPack.qc_searching.analysis.main import get_defect_state
import os

proj_path = '/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/search_triplet_from_defect_db'
db_json = os.path.join(proj_path, "db.json")

host_path = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc"
db_host_json = os.path.join(host_path, "db.json")


tot, proj, d_df = get_defect_state(
    db_json,
    {"task_id": 2183},
    2,-2,
    None,
    True,
    "dist",
    db_host_json,
    0.2
)

#%% plot dos
from atomate.vasp.database import VaspCalcDb
from pymatgen.electronic_structure.plotter import DosPlotter

proj_path = '/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/search_triplet_from_defect_db'
db = VaspCalcDb.from_db_file(os.path.join(proj_path, "db.json"))

dos = db.get_dos(2052)
dsplot = DosPlotter(zero_at_efermi=True)
dsplot.add_dos('tdos', dos)
plt = dsplot.get_plot([-10,10])
plt.show()
