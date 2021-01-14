#%% IMPORT
from fireworks import LaunchPad, Workflow

from atomate.vasp.database import VaspCalcDb
from atomate.vasp.workflows.jcustom.wf_full import get_wf_full_scan
from atomate.vasp.fireworks.core import ScanOptimizeFW
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    set_queue_options,
    set_execution_options,
    clean_up_files,
    add_modify_kpoints
)
from atomate.vasp.jpowerups import *

from pymatgen import Structure
from pymatgen.io.vasp.sets import MPScanRelaxSet, MPRelaxSet

from pycdt.core.defectsmaker import ChargedDefectsStructures

from qubitPack.tool_box import *

from monty.serialization import loadfn
import os
import pandas as pd


#%%
def relax_pc():
    lpad = LaunchPad.from_file(os.path.expanduser(
        os.path.join("~", "config/project/symBaseBinaryQubit/scan_relax_pc/my_launchpad.yaml")))

    mx2s = loadfn(os.path.expanduser(os.path.join("~", "config/project/symBaseBinaryQubit/"
                  "scan_relax_pc/gap_gt1-binary-NM.json")))

    for mx2 in mx2s:
        # if mx2["irreps"] and mx2["formula"] == "Rh2Br6" and mx2["spacegroup"] == "P3":
        if mx2["formula"] == "BN":
            pc = mx2["structure"]
            pc = modify_vacuum(pc, 30)
            scan_opt = ScanOptimizeFW(structure=pc, name="SCAN_relax")
            wf = Workflow([scan_opt], name="{}:SCAN_opt".format(mx2["formula"]))
            wf = add_modify_incar(wf)
            wf = add_modify_incar(
                wf,
                {
                    "incar_update": {
                        "LCHARG": False,
                        "LWAVE": False
                    }
                }
            )
            mx2.pop("structure")
            wf = add_additional_fields_to_taskdocs(
                wf,
                {"c2db_info": mx2}
            )
            wf = add_modify_incar(wf)
            wf = set_execution_options(wf, category="scan_relax_pc")
            wf = preserve_fworker(wf)
            lpad.add_wf(wf)
relax_pc()

#%%
def binary_scan_defect(defect_choice="substitutions", impurity_on_nn=None): #BN_vac

    col = VaspCalcDb.from_db_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/symBaseBinaryQubit/scan_relax_pc/db.json")).collection

    mx2s = col.find(
        {
            "nsites": 8
        }
    )
    geo_spec = None
    aexx = 0.25
    test = []
    for mx2 in mx2s[14:34]:
        pc = Structure.from_dict(mx2["output"]["structure"])
        cat = None

        if defect_choice == "substitutions":
            cat = "binary_defect_antisite_"
        elif defect_choice == "vacancies":
            cat = "bincary_defect_vac_"

        if mx2["nsites"] == 2:
            geo_spec = {25*2: [20]}
            cat += "AB"
        if mx2["nsites"] == 3:
            geo_spec = {25*3: [20]}
            cat += "AB2"
        if mx2["nsites"] == 4:
            geo_spec = {25*4: [20]}
            cat += "AB"
        if mx2["nsites"] == 8:
            geo_spec = {25*8: [20]}
            cat += "AB3"

        lpad = LaunchPad.from_file(
            os.path.join(
                os.path.expanduser("~"),
                "config/project/defect_db/{}/my_launchpad.yaml".format(cat)))
        print(cat)

        defect = ChargedDefectsStructures(pc, antisites_flag=True).defects

        cation, anion = find_cation_anion(pc)

        combo = mx2["sym_data"]["combo"]

        good_sym_site_symbols = get_good_ir_sites(combo["species"], combo["syms"])

        # Exclude those element in irrep
        # species = list(dict.fromkeys(mx2["sym_data"]["species"]))
        # for tgt in good_sym_site_symbols:
        #     species.remove(tgt)
        # good_sym_site_symbols = species

        print(good_sym_site_symbols)
        for good_sym_site_symbol in good_sym_site_symbols[0]:
            # substitutions or vacancies
            defect_type = (defect_choice, "{}".format(good_sym_site_symbol))

            for de_idx in range(len(defect[defect_type[0]])):
                print(cation, anion)

                if defect_type[1] == defect[defect_type[0]][de_idx]["name"].split("_")[-1]:
                    for na, thicks in geo_spec.items():
                        for thick in thicks:
                            for dtort in [0]:
                                if not impurity_on_nn:
                                    impurity_on_nn = []
                                gen_defect = GenDefect(pc, [defect_type[0], de_idx], na, thick, distort=dtort,
                                                       sub_on_side=list(impurity_on_nn))
                                # gen_defect.vacancies(dtort, list(impurity_on_nn))

                                defect_data = gen_defect.defect_entry

                                # add charge state regarding nelect
                                charge, nupdn = None, None
                                if MPRelaxSet(gen_defect.defect_st).nelect % 2 == 1:
                                    charge = [1,-1, 0]
                                    nupdn = [-1,-1,-1]
                                    print(charge)
                                    # if odd nelect, it used to be charge=[0]
                                else:
                                    # continue
                                    charge = [0]
                                    nupdn = [-1]

                                wf = get_wf_full_scan(
                                    structure=gen_defect.defect_st,
                                    charge_states=charge,
                                    gamma_only=False,
                                    gamma_mesh=True,
                                    dos=True,
                                    nupdowns=nupdn,
                                    vasptodb={
                                        "category": cat,
                                        "NN": gen_defect.NN,
                                        "NN_dist": gen_defect.nn_dist,
                                        "defect_entry": defect_data,
                                    },
                                    wf_addition_name="{}:{}".format(gen_defect.defect_st.num_sites, thick),
                                    task="scan_relax-scan_scf",
                                    category=cat,
                                )

                                def kpoints(kpts):
                                    kpoints_kwarg = {
                                        'comment': "Jcustom",
                                        "style": "G",
                                        "num_kpts": 0,
                                        'kpts': [kpts],
                                        'kpts_weights': None,
                                        'kpts_shift': (0, 0, 0),
                                        'coord_type': None,
                                        'labels': None,
                                        'tet_number': 0,
                                        'tet_weight': 0,
                                        'tet_connections': None
                                    }
                                    return kpoints_kwarg

                                # wf = add_modify_kpoints(wf, {"kpoints_update": kpoints([1,1,1])}, "PBE_relax")
                                # wf = add_modify_kpoints(wf, {"kpoints_update": kpoints([1,1,1])}, "HSE_relax")
                                # wf = add_modify_kpoints(wf, {"kpoints_update": kpoints([1,1,1])}, "HSE_scf")
                                wf = add_additional_fields_to_taskdocs(
                                    wf,
                                    {"pc_from": "symBaseBinaryQubit/scan_relax_pc/SCAN_relax/{}".format(mx2["task_id"]),
                                     "pc_from_id": mx2["task_id"],
                                     "combo": mx2["sym_data"]["combo"],
                                     "class": mx2["c2db_info"]["class"],
                                     "perturbed": gen_defect.distort}
                                )

                                wf = add_modify_incar(wf, {"incar_update": {"NSW":150, "LCHARG":False,
                                                                            "LWAVE":False}}, "SCAN_relax")
                                wf = add_modify_incar(wf, {"incar_update": {"LWAVE": True, "LCHARG": False}}, "SCAN_scf")
                                wf = set_queue_options(wf, "24:00:00")
                                # related to directory
                                wf = set_execution_options(wf, category=cat)
                                wf = preserve_fworker(wf)
                                wf.name = wf.name+":dx[{}]".format(gen_defect.distort)
                                # task_name_constraint=x meaning after x do powersup
                                wf = scp_files(
                                    wf,
                                    "/home/jengyuantsai/Research/projects/defect_db",
                                    fw_name_constraint="SCAN_scf",
                                )
                                wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb",
                                                    fw_name_constraint="SCAN_scf")

                                print(wf)
                                lpad.add_wf(wf)

binary_scan_defect()

#%% recalculate wavecar for successfull scan_scf's
def binary_scan_defect_gen_wavecar(): #BN_vac
    col = VaspCalcDb.from_db_file(
            os.path.expanduser(os.path.join("~",
            "config/project/defect_db/binary_defect/db.json"))).collection

    mx2s = col.find(
        {
            "task_label": "SCAN_scf",
            # "category": "binary_defect_antisite_AB3"
            "task_id": {"$in": [1190,1848, 2262, 2263, 2266, 2268, 2269, 2273, 2274, 2275]}}
    )
    geo_spec = None
    aexx = 0.25
    test = []

    for mx2 in list(mx2s)[:]: #max149
        lpad = LaunchPad.from_file(
            os.path.join(
                os.path.expanduser("~"),
                "config/project/defect_db/{}/my_launchpad.yaml".format(mx2["category"])))
        print(mx2["task_id"])
        wf = get_wf_full_scan(
            structure=Structure.from_dict(mx2["input"]["structure"]),
            charge_states=[0],
            gamma_only=False,
            gamma_mesh=True,
            dos=True,
            nupdowns=[-1],
            vasptodb={},
            wf_addition_name=mx2["task_id"],
            task="scan_scf",
            category=mx2["category"],
        )
        wf = write_inputs_from_db(
            wf,
            db_file=os.path.expanduser(
                os.path.join("~", "config/project/defect_db/binary_defect/db.json")),
            task_id=mx2["task_id"],
            modify_incar={"LWAVE":True, "LCHARG":False, "ICHARG":1}
        )
        wf = scp_files(wf, "/home/jengyuantsai/Research/projects/", "defect_db", mx2["dir_name"].split("/")[-1])
        # related to directory
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="SCAN_scf")
        wf = set_execution_options(wf, category=mx2["category"])
        wf = preserve_fworker(wf)
        # task_name_constraint=x meaning after x do powersup
        wf = clean_up_files(wf, files=["CHGCAR*"],
                            task_name_constraint="RunVasp")
        wf = clean_up_files(wf, files=["*"],
                            task_name_constraint="PassCalcLocs")
        print(wf)
        lpad.add_wf(wf)

binary_scan_defect_gen_wavecar()
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


#%% scan_relax_pc update
db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc/db.json")

for e in db.collection.find({}):
    st = Structure.from_dict(e["output"]["structure"])
    sym = SpacegroupAnalyzer(st, symprec=1e-1).get_symmetry_dataset()
    d = {}
    for k, v in sym.items():
        if type(v).__name__ == "ndarray":
            d[k] = v.tolist()
        else:
            d[k] = v
    d.update({"symprec":1e-1})
    species = [str(sp) for sp in st.species]
    combine = dict(zip(d["wyckoffs"], zip(species, d["site_symmetry_symbols"])))
    species, syms = get_good_ir_sites(dict(combine.values()).keys(), dict(combine.values()).values())
    combine.update({"species": species, "syms": syms})

    d.update({"combo": combine})
    print(d)
    db.collection.update_one({"task_id": e["task_id"]}, {"$set":{"sym_data":d}})


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
