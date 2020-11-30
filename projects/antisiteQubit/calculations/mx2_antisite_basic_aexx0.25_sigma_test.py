#%% Defect states
import os
from qubitPack.qc_searching.analysis.main import get_defect_state

proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_sigma_test"
# proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/perturbed"
db_json = os.path.join(proj_path, "db.json")

host_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_pc"
db_host_json = os.path.join(host_path, "db.json")


perturb = None
tid = 4271

tot, proj, d_df = get_defect_state(
    db_json,
    {"task_id": tid},
    1,
    -1,
    None,
    True,
    "dist",
    db_host_json,
    0.0
)
tot.to_clipboard()

#%% mx2_antisite_aexx0.25_sigma_test vs perturbed
"""
compare energy
"""
"""
Q: MoS2 has little splitting on x2-y2 and xy, and z2 is occupied. I don't think this is a ground state based on 
large energetic difference between C3v and Ch. Running new MoS2 hse_scf!
A: Mag of MoS2 is 4 not 2
"""
import os
import numpy as np
import pandas as pd
from atomate.vasp.database import VaspCalcDb
from pymatgen import Structure

p1 = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_sigma_test"
db1 = VaspCalcDb.from_db_file(os.path.join(p1, "db.json"))

p2 = "/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/perturbed"
db2 = VaspCalcDb.from_db_file(os.path.join(p2, "db.json"))

p3 = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_final"
db3 = VaspCalcDb.from_db_file(os.path.join(p3, "db.json"))

e1 = db1.collection.aggregate(
    [
        {"$match": {"task_label": "HSE_scf", "input.incar.EDIFF": 1E-5, "task_id": {"$ne":4268}}},
        {"$project": {
            "_id": 0,
            "energy": "$output.energy",
            "chem": "$chemsys",
            "tid": "$task_id",
            "a": "$output.structure.lattice.a",
            "unpert_st": "$output.structure",
            "unpert_nn": "$NN"
# "ediff": "$input.incar.EDIFF",
            # "mag": "$calcs_reversed.output.outcar.total_magnetization"
        },
        }
    ]
)
e2 = db2.collection.aggregate(
    [
        {"$match": {"task_label":"HSE_scf", "input.incar.EDIFF":1E-5, "task_id":{"$nin":[491, 488]}}},
        {"$project":{
            "_id":0,
            "energy_perturb": "$output.energy",
            "chem": "$chemsys",
            # "mag_pert": "$calcs_reversed.output.outcar.total_magnetization",
            "tid": "$task_id",
            "a": "$output.structure.lattice.a",
            "pert_st": "$output.structure",
            "pert_nn": "$NN"
        }
        },
    ]
)
# e3 = db3.collection.aggregate(
#     [
#         {"$match": {"task_label":"HSE_scf", "nupdown_set":2}},
#         {"$project":{
#             "_id":0,
#             "energy_no_sigma": "$output.energy",
#             "chem": "$chemsys",
#         }},
#     ]
# )


# e1 = pd.DataFrame(list(e1)).set_index("chem")
# e2 = pd.DataFrame(list(e2)).set_index("chem")
e1 = list(e1)
e2 = list(e2)
# e3 = pd.DataFrame(list(e3)).set_index("chem")
# e = pd.concat([e1,e2], axis=1).round(5)
# e["dE"] = e["energy"] - e["energy_perturb"]
# e.to_clipboard()

dd = []
for i in e1:
    data = {}
    for j in e2:
        if i["chem"] == j["chem"]:
            unpert_st = Structure.from_dict(i["unpert_st"])
            unpert_nn = i["unpert_nn"]
            pert_st = Structure.from_dict(j["pert_st"])
            pert_nn = j["pert_nn"]
            y = unpert_st.sites[unpert_nn[-1]]
            a1 = y.coords[:-1]
            x = pert_st.sites[pert_nn[-1]]
            a2 = x.coords[:-1]
            a = a1-a2
            d1 = np.linalg.norm(a1-a2)
            data.update({"chem":i["chem"],"d":d1, "dr":a})
            dd.append(data)
            print(i["chem"], j["chem"],a1,a2,a, d1)

pd.DataFrame(dd).to_clipboard()





#%% mx2_antisite_aexx0.25_sigma_test
"""
Q1 See energetic diff between EDIFF=1E-7 and 1E-5. A1 almost the same
Q2 Compare energy of MoS2 between triplet and mag=4
"""
from atomate.vasp.database import VaspCalcDb
import pandas as pd

p1 = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_sigma_test"
db1 = VaspCalcDb.from_db_file(os.path.join(p1, "db.json"))

e1 = db1.collection.aggregate(
    [
        {"$match": {"task_label": "HSE_scf", "input.incar.EDIFF":1E-5, "task_id":{"$ne":4268}}},
        {"$project": {
            "_id": 0,
            "energy": "$output.energy",
            "sigma": "$input.incar.SIGMA",
            "chem": "$chemsys",
            "tid": "$task_id",
            "ediff": "$input.incar.EDIFF",
            "NN":"$NN",
            "mag": "$calcs_reversed.output.outcar.total_magnetization"
        },
        }
    ]
)
e1 = list(e1)
print(pd.DataFrame(e1).sort_values("chem"))

#%% WF
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.powerups import *
from fireworks import LaunchPad
from pycdt.core.defectsmaker import ChargedDefectsStructures
from qubitPack.tool_box import find_cation_anion, GenDefect
from atomate.vasp.workflows.jcustom.wf_full import get_wf_full_hse
import os


CATEGORY = "mx2_antisite_basic_aexx0.25_sigma_test"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/category/{}/my_launchpad.yaml".format(CATEGORY)))


def wf(defect_type="substitutions"):

    col = VaspCalcDb.from_db_file("/home/tug03990/config/category/mx2_antisite_pc/db.json").collection
    # 4229: S-W, 4239: Se-W, 4236: Te-W, 4237:Mo-S, 4238: Mo-Se, 4235:Mo-Te
    mx2s = col.find({"task_id":{"$in":[4237]}})

    geo_spec = {5* 5 * 3: [20]}
    for mx2 in mx2s:
        pc = Structure.from_dict(mx2["output"]["structure"])

        defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
        cation, anion = find_cation_anion(pc)

        for sub in range(len(defect[defect_type])):
            print(cation, anion)
            # cation vacancy
            if "{}_on_{}".format(cation, anion) not in defect[defect_type][sub]["name"]:
                continue
            for na, thicks in geo_spec.items():
                for thick in thicks:
                    for dtort in [0]:
                        se_antisite = GenDefect(
                            orig_st=pc,
                            defect_type=(defect_type, sub),
                            natom=na,
                            vacuum_thickness=thick,
                            distort=dtort,
                        )

                        wf = get_wf_full_hse(
                            structure=se_antisite.defect_st,
                            charge_states=[0],
                            gamma_only=False,
                            gamma_mesh=True,
                            dos=True,
                            nupdowns=[0],
                            encut=320,
                            task="hse_relax-hse_scf",
                            vasptodb={"category": CATEGORY, "NN": se_antisite.NN,
                                      "defect_entry": se_antisite.defect_entry},
                            wf_addition_name="{}:{}".format(na, thick),
                            category=CATEGORY,
                        )

                        wf = add_additional_fields_to_taskdocs(
                            wf,
                            {
                                "lattice_constant": "HSE",
                                "perturbed": se_antisite.distort,
                                "sigma": 0.002,
                                "pc_from": "mx2_antisite_pc/HSE_scf/{}".format(mx2["task_id"])
                            }
                        )

                        wf = add_modify_incar(
                            wf,
                            {
                                "incar_update": {
                                    "EDIFFG":-0.02,
                                    "NSW":150,
                                    "SIGMA":0.002
                            }
                            },
                            "HSE_relax"
                        )
                        wf = add_modify_incar(wf, {"incar_update": {"LWAVE": True,
                                                                    "SIGMA":0.002}}, "HSE_scf")

                        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="PBE_relax")
                        wf = set_queue_options(wf, "48:00:00", fw_name_constraint="HSE_relax")
                        wf = set_queue_options(wf, "32:00:00", fw_name_constraint="HSE_scf")
                        wf = add_modify_incar(wf)
                        # related to directory
                        wf = set_execution_options(wf, category=CATEGORY)
                        wf = preserve_fworker(wf)
                        wf.name = wf.name+":dx[{}]:singlet".format(se_antisite.distort)
                        LPAD.add_wf(wf)

wf()