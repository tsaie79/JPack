import json

from fireworks import LaunchPad, Workflow

from atomate.vasp.fireworks.core import ScanOptimizeFW
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    set_queue_options,
    set_execution_options,
    clean_up_files,
    add_modify_kpoints,
    add_tags
)
from atomate.vasp.workflows.base.core import get_wf

from my_atomate_jyt.vasp.powerups import *
from my_atomate_jyt.vasp.workflows.wf_full import get_wf_full_scan

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet


from qubitPack.tool_box import *

from mpinterfaces.utils import *

from monty.serialization import loadfn
from monty.json import jsanitize

import os
import pandas as pd

def relax_pc():
    lpad = LaunchPad.from_file(os.path.expanduser(
        os.path.join("~", "config/project/symBaseBinaryQubit/scan_relax_pc/my_launchpad.yaml")))

    mx2s = loadfn(os.path.expanduser(os.path.join("~", "config/project/symBaseBinaryQubit/"
                                                       "scan_relax_pc/gap_gt1-binary-NM.json")))

    for mx2 in mx2s:
        # if mx2["irreps"] and mx2["formula"] == "Rh2Br6" and mx2["spacegroup"] == "P3":
        if mx2["formula"] == "BN":
            pc = mx2["structure"]
            pc = ensure_vacuum(pc, 30)
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

def pc():
    cat = "calc_data"
    lpad = LaunchPad.from_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/Scan2dMat/{}/my_launchpad.yaml".format(cat)))

    calculated_list = loadfn(os.path.join(
        os.path.dirname(os.path.abspath("__file__")),  "projects/defectDB/wf/Scan2dMat/calculated.json"))
    wait_list = loadfn(os.path.join(
        os.path.dirname(os.path.abspath("__file__")),  "projects/defectDB/wf/Scan2dMat/not_yet_calc.json"))
    col = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12345, user="readUser", password="qiminyan").collection

    mx2s = col.find(
        {
            "gap_hse_nosoc":{"$gt": 0},
            "nkinds": 2,
            "magstate": "NM",
            "uid": {"$in": wait_list},
        }
    )

    # col = get_db("symBaseBinaryQubit", "scan_relax_pc", port=12345, user="Jeng", password="qimin").collection
    # mx2s = col.find({"c2db_info.uid": {"$in": wait_list}})

    for idx, mx2 in enumerate(mx2s[100:]): #start from 100:!myq

        pc = Structure.from_dict(mx2["structure"])
        print(pc.num_sites)
        pc = ensure_vacuum(pc, 30)
        pc, site_info = phonopy_structure(pc)
        print(pc.num_sites)
        wf = get_wf(pc, os.path.join(os.path.dirname(os.path.abspath("__file__")), "projects/defectDB/wf/Scan2dMat/scan_pc.yaml"))

        wf = add_additional_fields_to_taskdocs(
            wf,
            {
                "c2db_info": {"uid": mx2["uid"]},
                "site_info": site_info
            },
            task_name_constraint="ToDb"
        )

        kpt = Kpoints.automatic_density_by_vol(pc, 144)
        wf = add_modify_kpoints(
            wf,
            modify_kpoints_params={"kpoints_update": {"kpts": kpt.kpts}},
            fw_name_constraint=wf.fws[1].name
        )
        wf = add_modify_2d_nscf_kpoints(
            wf,
            modify_kpoints_params={"kpoints_line_density": 20, "reciprocal_density": 144},
            fw_name_constraint=wf.fws[2].name
        )
        wf = add_modify_2d_nscf_kpoints(
            wf,
            modify_kpoints_params={"kpoints_line_density": 20, "reciprocal_density": 144, "mode": "uniform"},
            fw_name_constraint=wf.fws[3].name
        )

        if idx % 2 == 1:
            wf = set_execution_options(wf, category=cat, fworker_name="efrc")
        else:
            wf = set_execution_options(wf, category=cat, fworker_name="owls")
        wf = set_queue_options(wf, "06:00:00", fw_name_constraint=wf.fws[0].name)
        wf = set_queue_options(wf, "02:00:00", fw_name_constraint=wf.fws[1].name)
        wf = set_queue_options(wf, "01:00:00", fw_name_constraint=wf.fws[2].name)
        wf = set_queue_options(wf, "01:00:00", fw_name_constraint=wf.fws[3].name)
        wf = set_queue_options(wf, "01:00:00", fw_name_constraint=wf.fws[4].name)

        wf = clean_up_files(wf, files=["WAVECAR"], task_name_constraint="IRVSPToDb", fw_name_constraint="irvsp")
        wf = clean_up_files(wf, files=["CHGCAR*"], fw_name_constraint=wf.fws[2].name)
        wf = clean_up_files(wf, files=["CHGCAR*"], fw_name_constraint=wf.fws[3].name)

        wf = preserve_fworker(wf)
        wf = add_modify_incar(wf)
        wf.name = "{}:scan".format(mx2["uid"])
        lpad.add_wf(wf)

def binary_scan_defect(defect_choice="substitutions", impurity_on_nn=None):
    os.chdir("/home/tug03990/scripts/JPack/projects/defectDB")
    col = get_db("Scan2dMat", "calc_data", port=12347).collection

    # group_json_file = "wf/Scan2dDefect/bg_gte_1_and_inhomo_07272021"
    group_json_file = "wf/Scan2dDefect/bg_gte_1_and_homo_07272021"

    groups = loadfn(group_json_file)

    wfs = []
    geo_spec = None

    test_idx = 0
    for group in groups[1:2]:
        print(group)
        # if group["not_calculated_yet"] == None:
        #     continue
        tids = json.loads(group["tid"])

        group_id = int(group["group_id"])
        for tid in tids[:1]:
            print(group_id, tid)


            mx2 = col.find_one({"task_id": tid})
            pc = Structure.from_dict(mx2["output"]["structure"])
            print(pc.sites)
            good_ir_info = mx2["sym_data"]["good_ir_info"]
            cat = None

            if defect_choice == "substitutions":
                cat = "binary_sub_"
            elif defect_choice == "vacancies":
                cat = "binary_vac_"

            scaling = find_scaling_for_2d_defect(pc, 15)[0]
            area = scaling[0]*scaling[1]
            if mx2["nsites"] == 2:
                geo_spec = {area*2: [20]}
                cat += "AB"
            if mx2["nsites"] == 3:
                geo_spec = {area*3: [20]}
                cat += "AB2"
            if mx2["nsites"] == 4:
                geo_spec = {area*4: [20]}
                cat += "AB"
            if mx2["nsites"] == 6:
                geo_spec = {area*6: [20]}
                cat += "AB2"
            if mx2["nsites"] == 8:
                geo_spec = {area*8: [20]}
                cat += "AB3"


            print(cat)

            defects = ChargedDefectsStructures(pc, antisites_flag=True).defects

            # substitutions or vacancies!
            site_idxs = {}
            defect_type = defect_choice
            for de_idx, defect in enumerate(defects[defect_type]):
                for specie, site_sym, wy, site_idx in zip(*good_ir_info.values()):
                    if defect["unique_site"] == pc[site_idx]:
                        site_idxs[de_idx] = [pc[site_idx], specie, site_sym, wy, site_idx, defect["name"]]

            print(site_idxs)

            for de_idx, defect_site in site_idxs.items():
                for na, thicks in geo_spec.items():
                    for thick in thicks:
                        for dtort in [0]:
                            if not impurity_on_nn:
                                impurity_on_nn = []
                            gen_defect = GenDefect(pc, [defect_type, de_idx], na, thick, distort=dtort,
                                                   sub_on_side=list(impurity_on_nn))
                            # gen_defect.vacancies(dtort, list(impurity_on_nn))

                            def data_to_db():
                                d = {}
                                # defect data
                                defect_data = gen_defect.defect_entry.copy()
                                defect_data.update({
                                    "from_host_sym_data":{
                                        "site": defect_site[0],
                                        "specie": defect_site[1],
                                        "site_sym": defect_site[2],
                                        "wyckoff": defect_site[3],
                                        "site_idx": defect_site[4],
                                    }
                                })
                                d.update({
                                    "NN": gen_defect.NN,
                                    "NN_dist": gen_defect.nn_dist,
                                    "site_info": gen_defect.site_info,
                                    "perturbed": gen_defect.distort,
                                    "defect_entry": defect_data,
                                    "defect_name": defect_data["name"],
                                })
                                # host data
                                host_data = {
                                    "pc_from": "Scan2dMat/calc_data/{}".format(mx2["task_id"]),
                                    "pc_from_id": mx2["task_id"],
                                    "host_info":{"c2db_info": mx2["c2db_info"], "sym_data": mx2["sym_data"]},
                                }
                                locpot = mx2["calcs_reversed"][0]["output"]["locpot"]["2"]
                                vacuum = max(locpot)

                                scan_bs_output = mx2["output"].copy()
                                for remove in ["structure", "density", "energy", "energy_per_atom",
                                               "forces", "stress", "spacegroup"]:
                                    scan_bs_output.pop(remove)
                                band_edges = mx2["band_edges"]
                                scan_bs_output.update({"band_edges": band_edges, "vacuum_level": vacuum})
                                host_data["host_info"].update({"scan_bs": scan_bs_output})

                                d.update(host_data)
                                return d

                            # add charge state regarding nelect
                            charge, nupdn, n_of_e = None, None, None
                            if MPRelaxSet(gen_defect.defect_st).nelect % 2 == 1:
                                charge = [1]
                                nupdn = [-1]
                                n_of_e = "odd"
                                # if odd nelect, it used to be charge=[0]
                            else:
                                # continue
                                charge = [0]
                                nupdn = [-1]
                                n_of_e = "even"

                            wf = get_wf_full_scan(
                                structure=gen_defect.defect_st,
                                charge_states=charge,
                                gamma_only=True,
                                dos=False,
                                nupdowns=nupdn,
                                vasptodb={
                                },
                                wf_addition_name="{}:{}".format(gen_defect.defect_st.num_sites, thick),
                                wf_yaml="/home/tug03990/scripts/JPack/projects/defectDB/wf/Scan2dDefect/scan_defect"
                                        ".yaml"
                            )
                            wf.name += ":{}".format(gen_defect.defect_entry["name"])
                            # wf = add_modify_2d_nscf_kpoints(
                            #     wf,
                            #     modify_kpoints_params={"mode": "line"},
                            #     fw_name_constraint=wf.fws[2].name
                            # )
                            wf = add_modify_2d_nscf_kpoints(
                                wf,
                                modify_kpoints_params={"mode": "uniform"},
                                fw_name_constraint=wf.fws[2].name
                            )
                            wf = add_modify_incar(wf, {"incar_update": {"METAGGA":"R2SCAN"}})
                            wf = add_modify_incar(wf)
                            d = data_to_db()
                            d.update({
                                "category": cat,
                                "group_id": group_id,
                                "site_symmetry_uniform": False if "inhomo" in group_json_file else True
                            })
                            wf = add_additional_fields_to_taskdocs(
                                wf,
                                d,
                                task_name_constraint="ToDb"
                            )
                            wf = add_tags(wf, [{"pc_from_id": d["pc_from_id"],
                                                "group_id": group_id,
                                                "defect_name": d["defect_name"],
                                                "n_of_e": n_of_e,
                                                "charge_state": charge
                                                }])


                            wf = set_queue_options(wf, "12:00:00")
                            wf = set_queue_options(wf, "01:30:00", fw_name_constraint="irvsp")
                            wf = preserve_fworker(wf)
                            wf.name = wf.name+":dx[{}]".format(gen_defect.distort)
                            wf = bash_scp_files(
                                wf,
                                dest="/mnt/sdb/tsai/Research/projects/Scan2dDefect/calc_data/scf",
                                port=12348,
                                fw_name_constraint="SCAN_scf",
                                task_name_constraint="RunVasp"
                            )
                            # wf = clean_up_files(wf, files=["WAVECAR*"], task_name_constraint="SCP",
                            #                     fw_name_constraint=wf.fws[1].name)
                            wf = clean_up_files(wf, files=["CHGCAR*"], task_name_constraint="ToDb",
                                                fw_name_constraint=wf.fws[2].name)
                            wf = clean_up_files(wf, files=["WAVECAR"], task_name_constraint="IRVSPToDb",
                                                fw_name_constraint=wf.fws[3].name)

                            wfs.append(wf)
                            print(wf.name)

    return wfs

if __name__ == '__main__':
    lpad = LaunchPad.from_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/Scan2dDefect/calc_data/my_launchpad.yaml"))
    wfs = binary_scan_defect()
    for idx, wf in enumerate(wfs):
        if idx % 2 == 0:
            wf = set_execution_options(wf,  category="calc_data", fworker_name="owls")
        else:
            wf = set_execution_options(wf,  category="calc_data", fworker_name="efrc")
        # lpad.add_wf(wf)
        print(idx, wf.name,
              wf.fws[0].tasks[-1]["additional_fields"]["group_id"],
              wf.fws[0].tasks[-1]["additional_fields"]["pc_from_id"])

