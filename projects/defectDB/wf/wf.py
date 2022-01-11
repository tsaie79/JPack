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
    add_tags,
    add_modify_potcar,
    remove_custodian
)
from atomate.vasp.workflows.base.core import get_wf

from my_atomate_jyt.vasp.powerups import *
from my_atomate_jyt.vasp.workflows.wf_full import get_wf_full_scan, get_wf_full_hse
from my_atomate_jyt.vasp.fireworks.pytopomat import IrvspFW
from my_atomate_jyt.vasp.fireworks.pyzfs import PyzfsFW

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet

from qubitPack.tool_box import *

from mpinterfaces.utils import *

from monty.serialization import loadfn
from monty.json import jsanitize

import os
import pandas as pd

import sys
sys.path.append("/home/tug03990/scripts/JPack_independent")

from projects.antisiteQubit.wf_defect import ZPLWF

INPUT_PATH = "analysis/input"
OUTPUT_PATH = "analysis/output/xlsx"



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

def hse_relax_pc():
    cat = "hse_pc"
    qubit_df = IOTools(excel_file="qubit_2021-11-18", cwd="analysis/output/xlsx").read_excel()
    col = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12345, user="readUser", password="qiminyan").collection
    mx2s = col.find({"uid":{"$in":list(qubit_df["uid"].unique())}})
    # remember to write OPTCELL 110 in folders
    for mx2 in mx2s:
        if mx2["uid"] in [
            "In2Se2-GaSe-NM",
            "ZnCl2-GeS2-NM"
            "CaCl2-GeS2-NM",
            "ZnBr2-GeS2-NM",
            "SrBr2-GeS2-NM"
            "CdCl2-GeS2-NM",
            "CdBr2-GeS2-NM",
            "CdI2-GeS2-NM",
            "CdI2-CdI2-NM",
            "OsBr2-CdI2-NM",
            "PbI2-CdI2-NM",
            "SnBr2-MoS2-NM",
            "Sc2Br6-BiI3-NM",
            "Y2Cl6-AgBr3-NM",
            "MoTe2-MoS2-NM",
            "SnS2-CdI2-NM",
            "MoS2-MoS2-NM",
            "MoSe2-MoS2-NM",
            "WS2-MoS2-NM",
            "WSe2-MoS2-NM",
        ]:
            continue
        print(mx2["uid"])
        pc = Structure.from_dict(mx2["structure"])
        pc = set_vacuum(pc, 20)
        pc, site_info = phonopy_structure(pc)
        wf = get_wf_full_hse(
            structure=pc,
            charge_states=[0],
            gamma_only=False,
            gamma_mesh=True,
            nupdowns=[-1],
            task="hse_relax",
        )
        print(wf)

        wf = add_tags(wf,
                      [
                          {
                              "category": cat,
                              "uid": mx2["uid"],
                          }
                      ]
                      )

        wf = add_modify_incar(
            wf,
            {
                "incar_update":
                    {"EDIFF": 1E-6, "EDIFFG": -0.001, "ISIF": 3, "NSW": 250, "LCHARG":False,
                     "LWAVE":False, "LREAL": False}
            }
        )
        wf = add_additional_fields_to_taskdocs(wf, {"cat":cat, "uid": mx2["uid"]})
        wf = set_execution_options(wf, category=cat)
        wf = preserve_fworker(wf)
        lpad = LaunchPad.from_file(
            os.path.join(
                os.path.expanduser("~"),
                "config/project/HSE_triplets_from_Scan2dDefect/hse_pc/my_launchpad.yaml"))
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

def binary_triplet_HSE_wf(distort=0.0, category="calc_data", pyzfs_fw=True, irvsp_fw=True): #1e-4
    wfs = []
    # db_name, col_name = "HSE_triplets_from_Scan2dDefect", "hse_pc"
    db_name, col_name = "2dMat_from_cmr_fysik", "2dMaterial_v1"
    hse_pc_db = get_db(db_name, col_name,  port=12345, user="readUser", password="qiminyan")
    pc_db = hse_pc_db
    # defect_df = IOTools(cwd=INPUT_PATH, excel_file="defect_2021-11-18").read_excel()
    # triplet_df = defect_df.loc[(defect_df["mag"] == 2), ["uid", "charge", "defect_name", "defect_type", "task_id"]]
    # triplet_df = triplet_df.iloc[24:34, :]
    # triplet_df = triplet_df.iloc[12:, :]
    qubit_candidate_df = IOTools(cwd=INPUT_PATH, excel_file="qubit_2021-11-18").read_excel()
    triplet_df = qubit_candidate_df.copy()
    triplet_df = triplet_df.loc[triplet_df["uid"].isin(["MoTe2-MoS2-NM", "WTe2-MoS2-NM"])]
    #"WTe2-MoS2-NM", "MoS2-MoS2-NM", "MoSe2-MoS2-NM",
    # "MoTe2-MoS2-NM"])]
    for uid, charge, defect_name, defect_type, task_id in zip(triplet_df["uid"], triplet_df["charge"],
                                                     triplet_df["defect_name"], triplet_df["defect_type"],
                                                              triplet_df["task_id"]
                                                     ):
        mx2 = pc_db.collection.find_one({"uid": uid})
        pc = Structure.from_dict(mx2["structure"])
        scaling = find_scaling_for_2d_defect(pc, 15)[0]
        area = scaling[0]*scaling[1]
        geo_spec = {area*pc.num_sites: [20]}

        defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
        cation, anion = find_cation_anion(pc)

        if defect_type == "vacancy":
            defect_type = "vacancies"
        else:
            defect_type = "substitutions"
        for sub in range(len(defect[defect_type])):
            print(cation, anion)
            print(defect[defect_type][sub]["name"])
            if defect_name not in defect[defect_type][sub]["name"]:
                continue
            for na, thicks in geo_spec.items():
                for thick in thicks:
                    print(sub, na, defect_type, thick, distort)
                    gen_defect = GenDefect(
                        orig_st=pc,
                        defect_type=(defect_type, sub),
                        natom=na,
                        vacuum_thickness=thick,
                        distort=distort,
                    )
                    magmom = MPRelaxSet(gen_defect.defect_st).incar.get("MAGMOM")
                    print(magmom)
                    magmom[55] = -5
                    # special_st = Structure.from_file("/home/tug03990/wse0_456.vasp")
                    # move_site(special_st, [25], dz)
                    wf = get_wf_full_hse(
                        structure=gen_defect.defect_st,
                        task_arg=dict(lcharg=True),
                        charge_states=[charge, charge],
                        gamma_only=False,
                        gamma_mesh=True,
                        nupdowns=[2, 0],
                        task="hse_relax-hse_scf",
                        vasptodb={
                            "category": category, "NN": gen_defect.NN,
                            "defect_entry": gen_defect.defect_entry,
                            "lattice_constant": "HSE",
                            "taskid_Scan2dDefect": task_id,
                            "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["uid"]),
                        },
                        wf_addition_name="{}:{}".format(na, thick),
                    )

                    if irvsp_fw:
                        HSE_fws = wf.fws.copy()
                        for idx, fw in enumerate(HSE_fws):
                            if "HSE_scf" in fw.name:
                                irvsp_fw = IrvspFW(structure=gen_defect.defect_st,
                                                   kpt_mode="single_kpt", symprec=0.001,
                                                   irvsptodb_kwargs=dict(collection_name="ir_data"),
                                                   name="irvsp",
                                                   parents=fw)
                                HSE_fws.append(irvsp_fw)
                        wf = Workflow(HSE_fws, name=wf.name)

                    if pyzfs_fw:
                        HSE_fws = wf.fws.copy()
                        for idx, fw in enumerate(HSE_fws):
                            if "HSE_scf" in fw.name and fw.tasks[5]["incar_update"]["NUPDOWN"] == 2:
                                pyzfs_fw = PyzfsFW(structure=gen_defect.defect_st, parents=fw, pyzfstodb_kwargs=dict(
                                    collection_name="zfs_data"))
                                HSE_fws.append(pyzfs_fw)
                        wf = Workflow(HSE_fws, name=wf.name)


                    wf = add_tags(wf,
                                  [
                                      {
                                          "category": category,
                                          "uid": mx2["uid"],
                                      }
                                  ]
                                  )

                    wf = bash_scp_files(
                        wf,
                        dest="/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/calc_data-pbe_pc/",
                        port=12348,
                        task_name_constraint="RunVasp"
                    )
                    wf = set_execution_options(wf, category=category, fw_name_constraint="HSE_relax")
                    wf = set_execution_options(wf, category=category, fw_name_constraint="HSE_scf")

                    if irvsp_fw:
                        wf = set_execution_options(wf, category="ir_data", fw_name_constraint="irvsp")
                        wf = bash_scp_files(
                            wf,
                            dest="/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/ir_data-pbe_pc/",
                            port=12348,
                            task_name_constraint="ToDb",
                            fw_name_constraint="irvsp"
                        )
                        wf = clean_up_files(wf, files=["WAVECAR*", "vasprun.xml*"], fw_name_constraint="irvsp",
                                            task_name_constraint="ToDb")


                    if pyzfs_fw:
                        wf = set_execution_options(wf, category="zfs_data", fw_name_constraint="pyzfs")
                        wf = bash_scp_files(
                            wf,
                            dest="/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/zfs_data-pbe_pc/",
                            port=12348,
                            task_name_constraint="ToDb",
                            fw_name_constraint="pyzfs"
                        )
                        wf = clean_up_files(wf, files=["WAVECAR*", "vasprun.xml*"], fw_name_constraint="pyzfs",
                                            task_name_constraint="ToDb")
                    wf = add_modify_incar(wf, {"incar_update": {"EDIFF": 1E-5, "MAGMOM": magmom}})
                    wf = add_modify_incar(wf)
                    wf.name = "{}:{}:{}".format(mx2["uid"], defect_name, ":".join(wf.name.split(":")[3:5]))
                    print(wf.name)
                    wfs.append(wf)
    return wfs

def cdft(task, std_taskid, nbands=None, occ=None, std_base_dir=None, category="cdft", specific_poscar=None,
             prevent_JT=True):

    std_db = get_db("HSE_triplets_from_Scan2dDefect", "calc_data-pbe_pc", port=12347)
    std_entry = std_db.collection.find_one({"task_id":std_taskid})

    std_path = std_entry["calcs_reversed"][0]["dir_name"]
    std_launchdir = std_path.split("/")[-1]
    std_path = os.path.join("/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/calc_data-pbe_pc/", std_launchdir)
    if not std_base_dir:
        std_base_dir = std_path.split("launch")[0]

    cdft = ZPLWF(std_entry, std_path, None)

    if not nbands:
        nbands = std_entry["input"]["parameters"]["NBANDS"]

    up_occ, dn_occ = None, None
    up_occ_no_excite, dn_occ_no_excite = get_lowest_unocc_band_idx_v2(std_taskid, std_db, nbands)
    wf = None
    if occ["1"] and occ["-1"] == "None":
        up_occ = occ["1"]
        dn_occ = dn_occ_no_excite
        print("up_occ: {}, dn_occ:{}".format(up_occ, dn_occ))
        wf = cdft.wf(
            task, 0, up_occupation=up_occ,
            down_occupation=dn_occ, nbands=nbands, gamma_only=True, selective_dyn=None,
            specific_structure=specific_poscar,
            nonsoc_prev_dir=None
        )
        wf = add_additional_fields_to_taskdocs(wf, {"cdft_occ": {"up": up_occ, "dn": dn_occ}})

    elif occ["-1"] and occ["1"] == "None":
        dn_occ = occ["-1"]
        up_occ = up_occ_no_excite
        print("up_occ: {}, dn_occ:{}".format(up_occ, dn_occ))
        wf = cdft.wf(
            task, 0, up_occupation=up_occ,
            down_occupation=dn_occ, nbands=nbands, gamma_only=True, selective_dyn=None,
            specific_structure=specific_poscar,
            nonsoc_prev_dir=None
        )
        wf = add_additional_fields_to_taskdocs(wf, {"cdft_occ": {"up": up_occ, "dn": dn_occ}})

    elif occ["1"] and occ["-1"]:
        up_occ = occ["1"]
        dn_occ = dn_occ_no_excite
        print("up_occ: {}, dn_occ:{}".format(up_occ, dn_occ))
        wf = cdft.wf(
            task, 0, up_occupation=up_occ,
            down_occupation=dn_occ, nbands=nbands, gamma_only=True, selective_dyn=None,
            specific_structure=specific_poscar,
            nonsoc_prev_dir=None
        )
        wf = add_additional_fields_to_taskdocs(wf, {"cdft_occ": {"up": up_occ, "dn": dn_occ}})
        fws_up = wf.fws.copy()


        dn_occ = occ["-1"]
        up_occ = up_occ_no_excite
        print("up_occ: {}, dn_occ:{}".format(up_occ, dn_occ))
        wf = cdft.wf(
            task, 0, up_occupation=up_occ,
            down_occupation=dn_occ, nbands=nbands, gamma_only=True, selective_dyn=None,
            specific_structure=specific_poscar,
            nonsoc_prev_dir=None
        )
        wf = add_additional_fields_to_taskdocs(wf, {"cdft_occ": {"up": up_occ, "dn": dn_occ}})
        fws_dn = wf.fws.copy()

        fws_up.extend(fws_dn)
        wf = Workflow(fws_up, name=wf.name)



    wf = clean_up_files(wf, files=["WAVECAR*", "vasprun.xml*"], fw_name_constraint="CDFT-B-HSE_scf",
                        task_name_constraint="ToDb")
    wf = clean_up_files(wf, files=["WAVECAR*", "vasprun.xml*"], fw_name_constraint="CDFT-D-HSE_scf",
                        task_name_constraint="ToDb")

    wf = bash_scp_files(
        wf,
        dest="/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/cdft-pbe_pc/",
        port=12348,
        task_name_constraint="ToDb",
    )


    wf = add_additional_fields_to_taskdocs(
        wf,
        {
            "source_entry":"{}/{}/{}".format(std_entry["db"],
                                             std_entry["collection"],
                                             std_entry["task_id"]),
            "prev_fw_taskid": std_entry["task_id"],
            "prev_fw_db": std_entry["db"],
            "prev_fw_collection": std_entry["collection"]
        })
    wf = add_modify_incar(wf)
    wf = set_execution_options(wf, category=category)
    wf = preserve_fworker(wf)
    return wf, category

def main():
    def run_binary_scan_defect():
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

    def run_triplet_HSE_wf():
        lpad = LaunchPad.from_file(
            os.path.join(
                os.path.expanduser("~"),
                "config/project/HSE_triplets_from_Scan2dDefect/calc_data/my_launchpad.yaml"))
        wfs = binary_triplet_HSE_wf(irvsp_fw=True)
        for idx, wf in enumerate(wfs[:]):
            wf.name += "Ch"
            wf = set_execution_options(wf, fworker_name="owls")
            wf = preserve_fworker(wf)
            lpad.add_wf(wf)
            print(idx, wf.name)

    def run_cdft():
        # qubit_df = IOTools(excel_file="hse_screened_qubit_2021-11-18", cwd=INPUT_PATH).read_excel()
        qubit_df = IOTools(excel_file="hse_screened_qubit_2021-11-18_nersc_modified-nbands",
                           cwd=INPUT_PATH).read_excel()
        fworker = "nersc"
        taskids = qubit_df.loc[qubit_df["fworker"]==fworker, "task_id"]
        # taskids = qubit_df.loc[qubit_df["task_id"]==868, "task_id"]
        if fworker == "nersc":
            # hse qubit on nersc uses gpu-vasp6, which has problem of running cdft!
            fworker = "owls"
        # fworker = "efrc"
        # taskids = qubit_df.loc[qubit_df["fworker"]==fworker, "task_id"][:]
        for taskid in taskids[1:]:
            print("{}=".format(taskid)*20)
            nbands = qubit_df.loc[qubit_df["task_id"] == taskid, "nbands"].iloc[0]
            up_occ = qubit_df.loc[qubit_df["task_id"] == taskid, "up_tran_cdft_occ"].iloc[0]
            dn_occ = qubit_df.loc[qubit_df["task_id"] == taskid, "dn_tran_cdft_occ"].iloc[0]
            wf, cat = cdft(
                "B-C-D",
                taskid,
                occ={"1": up_occ, "-1": dn_occ},
                nbands=nbands,
                category="cdft",
            )
            wf = remove_custodian(wf)
            wf = add_additional_fields_to_taskdocs(wf, {"cat":cat, "taskid": taskid})
            wf = set_execution_options(wf, fworker_name=fworker)
            lpad = LaunchPad.from_file(
                os.path.join(
                    os.path.expanduser("~"),
                    "config/project/HSE_triplets_from_Scan2dDefect/cdft/my_launchpad.yaml"))
            lpad.add_wf(wf)

    run_triplet_HSE_wf()
if __name__ == '__main__':
    main()