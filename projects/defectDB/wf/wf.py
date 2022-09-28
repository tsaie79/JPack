import json
import shutil
import subprocess

import monty.serialization

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
    remove_custodian,
    use_custodian
)
from atomate.vasp.workflows.base.core import get_wf

from my_atomate_jyt.vasp.powerups import *
from my_atomate_jyt.vasp.workflows.wf_full import get_wf_full_scan, get_wf_full_hse
from my_atomate_jyt.vasp.fireworks.pytopomat import IrvspFW
from my_atomate_jyt.vasp.fireworks.pyzfs import PyzfsFW

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.core.units import eV_to_Ha

from ase.units import _hplanck, _eps0, _c
from VaspBandUnfolding.vasp_constant import *
DEBYETOSI = 3.335640952e-30

from qubitPack.tool_box import *

from mpinterfaces.utils import *

from monty.serialization import loadfn
from monty.json import jsanitize

import os, sys
from pathlib import Path
import pandas as pd

from JPack_independent.projects.antisiteQubit.wf_defect import ZPLWF

# from analysis.data_analysis import Defect

MODULE_DIR = Path("__file__").resolve().parent
print(MODULE_DIR)

INPUT_PATH = os.path.join(MODULE_DIR, "analysis/input")
OUTPUT_PATH = os.path.join(MODULE_DIR, "analysis/output/xlsx")

HSEQubitDefect = get_db("HSE_triplets_from_Scan2dDefect", "calc_data-pbe_pc", port=12347, user="Jeng")

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


def grand_hse_defect(defect_choice="vacancies", impurity_on_nn=None, pyzfs_fw=True, irvsp_fw=True): #1e-4
    fworker = "efrc" #"gpu_nersc"
    defect_db_name = "C2DB_IR_vacancy_HSE"

    col = get_db("C2DB_IR", "calc_data", port=12345).collection

    wfs = []
    geo_spec = None
    pc_entries = [e for e in col.find(
        {
            "task_label": "hse line",
            "nelements": {"$lte": 2},
            "output.bandgap": {"$gte":1.0},
            "sym_data.good_ir_info.species": {"$ne": []},
        }
    )]
    for pc_entry in pc_entries[316:]: #345 pcs, index: 0-d344
        print(f"-------{pc_entry['task_id']}"*10)
        mx2 = pc_entry.copy()
        pc = Structure.from_dict(mx2["output"]["structure"])
        good_ir_info = mx2["sym_data"]["good_ir_info"]

        scaling = find_scaling_for_2d_defect(pc, 15)[0]
        area = scaling[0] * scaling[1]
        if mx2["nsites"] == 2:
            geo_spec = {area * 2: [20]}
        if mx2["nsites"] == 3:
            geo_spec = {area * 3: [20]}
        if mx2["nsites"] == 4:
            geo_spec = {area * 4: [20]}
        if mx2["nsites"] == 6:
            geo_spec = {area * 6: [20]}
        if mx2["nsites"] == 8:
            geo_spec = {area * 8: [20]}

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
                        try:
                            gen_defect = GenDefect(
                                pc, [defect_type, de_idx], na, thick, distort=dtort,
                                sub_on_side=list(impurity_on_nn)
                                )
                        except Exception as e:
                            print("error", type(e.__str__()))
                            from monty.serialization import loadfn, dumpfn
                            path = "/home/tug03990/scripts/JPack_independent/projects/defectDB"
                            taskid_error = loadfn(f"{path}/wf/C2DB_IR_vacancy_HSE/error.json")
                            print(taskid_error)
                            taskid_error.update({
                                str(mx2["task_id"]): {"pycdt_index": de_idx, "message": e.__str__()}})
                            dumpfn(taskid_error, f"{path}/wf/C2DB_IR_vacancy_HSE/error.json", indent=4)
                            continue


                        # add charge state depending on nelect
                        charge, nupdn, n_of_e = None, None, None
                        if MPRelaxSet(gen_defect.defect_st).nelect % 2 == 1:
                            charge = [1, 0, -1]
                            nupdn = [-1, -1, -1]
                            n_of_e = "odd"
                            # if odd nelect, it used to be charge=[0]
                        else:
                            # continue
                            charge = [0]
                            nupdn = [-1]
                            n_of_e = "even"


                        def data_to_db():
                            d = {}
                            # defect data
                            defect_data = gen_defect.defect_entry.copy()
                            defect_data.update(
                                {
                                    "from_host_sym_data": {
                                        "site": defect_site[0],
                                        "specie": defect_site[1],
                                        "site_sym": defect_site[2],
                                        "wyckoff": defect_site[3],
                                        "site_idx": defect_site[4],
                                    }
                                }
                            )
                            d.update(
                                {
                                    "NN": gen_defect.NN,
                                    "NN_dist": gen_defect.nn_dist,
                                    "site_info": gen_defect.site_info,
                                    "perturbed": gen_defect.distort,
                                    "defect_entry": defect_data,
                                    "defect_name": defect_data["name"],
                                }
                            )
                            # host data
                            host_data = {
                                "pc_from": "C2DB_IR/calc_data/{}".format(mx2["task_id"]),
                                "pc_from_id": mx2["task_id"],
                                "host_info": {"sym_data": mx2["sym_data"]},
                            }

                            hse_bs_output = mx2["output"].copy()
                            for remove in ["structure", "density", "energy", "energy_per_atom",
                                           "forces", "stress", "spacegroup"]:
                                hse_bs_output.pop(remove)
                            band_edges = mx2.get("band_edges", None)
                            hse_bs_output.update({"band_edges": band_edges})
                            host_data["host_info"].update({"c2db_ir_hse_line": hse_bs_output})

                            d.update(host_data)
                            return d

                        wf = get_wf_full_hse(
                            structure=gen_defect.defect_st,
                            task_arg=dict(lcharg=True),
                            charge_states=charge,
                            gamma_only=False,
                            gamma_mesh=True,
                            nupdowns=nupdn,
                            task="hse_relax-hse_scf",
                            vasptodb={},
                            wf_addition_name="{}:{}".format(na, thick),
                        )

                        d = data_to_db()
                        wf = add_additional_fields_to_taskdocs(wf, d, task_name_constraint="ToDb")
                        wf = add_tags(
                            wf, [{
                                     "pc_from_id": d["pc_from_id"],
                                     "defect_name": d["defect_name"],
                                     "n_of_e": n_of_e,
                                     "charge_state": charge
                                 }]
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
                            wf = set_execution_options(wf, category="ir_data", fworker_name=fworker,
                                                       fw_name_constraint="irvsp")

                        if pyzfs_fw:
                            HSE_fws = wf.fws.copy()
                            for idx, fw in enumerate(HSE_fws):
                                if "HSE_scf" in fw.name:
                                    pyzfs_fw = PyzfsFW(structure=gen_defect.defect_st,
                                                       parents=fw,
                                                       pyzfstodb_kwargs=dict(collection_name="zfs_data"))
                                    HSE_fws.append(pyzfs_fw)
                            wf = Workflow(HSE_fws, name=wf.name)
                            wf = set_execution_options(wf, category="zfs_data", fworker_name=fworker,
                                                       fw_name_constraint="pyzfs")

                        wf = bash_scp_files(
                            wf,
                            dest=f"/mnt/sdc/tsai/Research/projects/{defect_db_name}/calc_data/",
                            port=12348,
                            task_name_constraint="RunVasp"
                        )
                        wf = set_execution_options(wf, category="calc_data", fw_name_constraint="HSE_relax",
                                                   fworker_name=fworker)
                        wf = set_execution_options(wf, category="calc_data", fw_name_constraint="HSE_scf",
                                                   fworker_name=fworker)


                        if irvsp_fw:
                            wf = bash_scp_files(
                                wf,
                                dest=f"/mnt/sdc/tsai/Research/projects/{defect_db_name}/ir_data/",
                                port=12348,
                                task_name_constraint="ToDb",
                                fw_name_constraint="irvsp"
                            )
                            wf = clean_up_files(wf, files=["WAVECAR*", "vasprun.xml*"], fw_name_constraint="irvsp",
                                                task_name_constraint="ToDb")


                        if pyzfs_fw:
                            wf = bash_scp_files(
                                wf,
                                dest=f"/mnt/sdc/tsai/Research/projects/{defect_db_name}/zfs_data/",
                                port=12348,
                                task_name_constraint="ToDb",
                                fw_name_constraint="pyzfs"
                            )
                            wf = clean_up_files(wf, files=["WAVECAR*", "vasprun.xml*"], fw_name_constraint="pyzfs",
                                                task_name_constraint="ToDb")

                        wf = add_modify_incar(wf, {"incar_update": {"SIGMA": 0.01, "EDIFF": 1E-4}},
                                              fw_name_constraint="HSE_relax")
                        wf = add_modify_incar(wf, {"incar_update": {"SIGMA": 0.01}}, fw_name_constraint="HSE_scf")
                        wf = add_modify_incar(wf)

                        wf.name = f"{pc_entry['task_id']}:{gen_defect.defect_entry['name']}:" \
                                  f"{':'.join(wf.name.split(':')[3:5])}"
                        print(wf.name)
                        wfs.append(wf)
    return wfs


def binary_scan_defect(defect_choice="substitutions", impurity_on_nn=None):
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


    def antisite_TMD_pc():
        """
        1. HSE-relax pc: db_name, col_name = "owls", "mx2_antisite_pc" (in our Nature comm.)
        * taskids: 3083: Se-W, 3091: S-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te

        2. C2DB PBE pc:  db_name, col_name = "2dMat_from_cmr_fysik", "2dMaterial_v1"
        * uid: "MoS2-MoS2-NM", "MoSe2-MoS2-NM", "MoTe2-MoS2-NM", "WS2-MoTe2-NM", "WSe2-MoTe2-NM", "WTe2-MoS2-NM",
        :return:
        """
        db_name, col_name = "2dMat_from_cmr_fysik", "2dMaterial_v1"

        col = get_db(db_name, col_name, port=12345, user="readUser", password="qiminyan").collection
        taskids = ["WTe2-MoS2-NM"]
        mx2s = col.find({"uid":{"$in": taskids}})
        pcs = [Structure.from_dict(mx2["structure"]) for mx2 in mx2s]
        defect_inputs = {"pc": [], "defect_name": [], "charge": [], "defect_type": [], "scan_taskids":[], "db_name":
            [], "col_name": [], "taskid": [], "spin": []}
        for pc, taskid in zip(pcs, taskids):
            cation, anion = find_cation_anion(pc)
            defect_inputs["pc"].append(pc)
            defect_inputs["defect_name"].append("as_1_{}_on_{}".format(cation, anion))
            defect_inputs["defect_type"].append("antisite")
            defect_inputs["scan_taskids"].append(None)
            defect_inputs["db_name"].append(db_name)
            defect_inputs["col_name"].append(col_name)
            defect_inputs["taskid"].append(taskid)
            defect_inputs["charge"].append(0)
            defect_inputs["spin"].append(0)
            df = pd.DataFrame(defect_inputs)
        return df


    def group_36_pc():
        db_name, col_name = "2dMat_from_cmr_fysik", "2dMaterial_v1"
        col = get_db(db_name, col_name,  port=12345, user="readUser", password="qiminyan")
        taskids = ["In2Te2-GaSe-NM"]
        mx2s = col.find({"uid": {"$in": taskids}})
        pcs = [Structure.from_dict(mx2["structure"]) for mx2 in mx2s]
        defect_inputs = {"pc": [], "defect_name": [], "charge": [], "defect_type": [], "scan_taskids":[], "db_name": [],
                         "col_name": [], "taskid": []}
        for pc, taskid in zip(pcs, taskids):
            cation, anion = find_cation_anion(pc)
            defect_inputs["pc"].append(pc)
            defect_inputs["defect_name"].append("as_1_{}_on_{}".format(cation, anion))
            defect_inputs["charge"].append(-1)
            defect_inputs["defect_type"].append("antisite")
            defect_inputs["scan_taskids"].append(None)
            defect_inputs["db_name"].append(db_name)
            defect_inputs["col_name"].append(col_name)
            defect_inputs["taskid"].append(taskid)

            defect_inputs["pc"].append(pc)
            defect_inputs["defect_name"].append("as_1_{}_on_{}".format(anion, cation))
            defect_inputs["charge"].append(1)
            defect_inputs["defect_type"].append("antisite")
            defect_inputs["scan_taskids"].append(None)
            defect_inputs["db_name"].append(db_name)
            defect_inputs["col_name"].append(col_name)
            defect_inputs["taskid"].append(taskid)

            df = pd.DataFrame(defect_inputs)
        return df, db_name, col_name



    # defect_df = IOTools(cwd=INPUT_PATH, excel_file="defect_2021-11-18").read_excel()
    # triplet_df = defect_df.loc[(defect_df["mag"] == 2), ["uid", "charge", "defect_name", "defect_type", "task_id"]]
    # triplet_df = triplet_df.iloc[24:34, :]
    # triplet_df = triplet_df.iloc[12:, :]
    # qubit_candidate_df = IOTools(cwd=INPUT_PATH, excel_file="triplet_2021-11-18").read_excel()
    # triplet_df = qubit_candidate_df.copy()
    # triplet_df = triplet_df.loc[triplet_df["task_id"].isin([920, 1080, 866])]

    # for uid, charge, defect_name, defect_type, task_id in zip(triplet_df["uid"], triplet_df["charge"],
    #                                                  triplet_df["defect_name"], triplet_df["defect_type"],
    #                                                           triplet_df["task_id"]
    #                                                  ):
    pc_df = antisite_TMD_pc()
    wfs = []
    for pc, charge, spin, defect_name, defect_type, scan2dDefect_taskid, db_name, col_name, pc_taskid in zip(
            pc_df["pc"], pc_df["charge"], pc_df["spin"], pc_df["defect_name"], pc_df["defect_type"],
            pc_df["scan_taskids"], pc_df["db_name"], pc_df["col_name"], pc_df["taskid"]):
        print(f"@@@@@@@@@@@@{defect_name} {charge} {spin}")
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

                    wf = get_wf_full_hse(
                        structure=gen_defect.defect_st,
                        task_arg=dict(lcharg=True),
                        charge_states=[charge],
                        gamma_only=False,
                        gamma_mesh=True,
                        nupdowns=[spin],
                        task="hse_relax-hse_scf",
                        vasptodb={
                            "category": category, "NN": gen_defect.NN,
                            "defect_entry": gen_defect.defect_entry,
                            "lattice_constant": "HSE",
                            "taskid_Scan2dDefect": scan2dDefect_taskid,
                            "pc_from": "{}/{}/{}".format(db_name, col_name, pc_taskid),
                        },
                        wf_addition_name="{}:{}".format(na, thick),
                    )

                    if irvsp_fw:
                        HSE_fws = wf.fws.copy()
                        for idx, fw in enumerate(HSE_fws):
                            if "HSE_scf" in fw.name:
                                irvsp_fw = IrvspFW(structure=gen_defect.defect_st,
                                                   kpt_mode="single_kpt", symprec=0.001,
                                                   irvsptodb_kwargs=dict(collection_name="ir_data-pbe_pc"),
                                                   name="irvsp",
                                                   parents=fw)
                                HSE_fws.append(irvsp_fw)
                        wf = Workflow(HSE_fws, name=wf.name)

                    if pyzfs_fw:
                        HSE_fws = wf.fws.copy()
                        for idx, fw in enumerate(HSE_fws):
                            if "HSE_scf" in fw.name and fw.tasks[5]["incar_update"]["NUPDOWN"] == 2:
                                pyzfs_fw = PyzfsFW(structure=gen_defect.defect_st, parents=fw, pyzfstodb_kwargs=dict(
                                    collection_name="zfs_data-pbe_pc"))
                                HSE_fws.append(pyzfs_fw)
                        wf = Workflow(HSE_fws, name=wf.name)


                    wf = add_tags(wf,
                                  [
                                      {
                                          "category": category,
                                          "uid": pc_taskid,
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
                    wf = add_modify_incar(wf)
                    wf.name = "{}:{}:{}".format(pc_taskid, defect_name, ":".join(wf.name.split(":")[3:5]))
                    print(wf.name)
                    wfs.append(wf)
    return wfs

def cdft(task, std_taskid, nbands=None, occ=None, std_base_dir=None, category="cdft", specific_poscar=None,
             prevent_JT=True):
    std_cat_name = "calc_data-pbe_pc"
    std_db = get_db("HSE_triplets_from_Scan2dDefect", std_cat_name, port=12347)
    std_entry = std_db.collection.find_one({"task_id":std_taskid})

    std_path = std_entry["calcs_reversed"][0]["dir_name"]
    std_launchdir = std_path.split("/")[-1]
    std_path = os.path.join("/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/", std_cat_name,
                            std_launchdir)
    if not std_base_dir:
        std_base_dir = std_path.split("launch")[0]

    cdft = ZPLWF(std_entry, std_path, None)

    if not nbands:
        nbands = std_entry["input"]["parameters"]["NBANDS"]

    up_occ, dn_occ = None, None
    up_occ_no_excite, dn_occ_no_excite = get_lowest_unocc_band_idx_v2(std_taskid, std_db, nbands)
    wf = None
    if occ["1"] and pd.isna(occ["-1"]):
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

    elif occ["-1"] and pd.isna(occ["1"]):
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
        dest="/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/cdft-pbe_pc",
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

def cohp(structure, scf_path, lobsterin_key_dict=None):
    from atomate.vasp.fireworks.lobster import LobsterFW
    from qubitPack.tool_box import get_db
    from fireworks import LaunchPad, Workflow
    category = "cohp"
    # validator_group="no_validator" for running lobster with too few bands
    fw = LobsterFW(structure=structure, prev_calc_dir=scf_path, lobstertodb_kwargs={},
                   lobsterin_key_dict=lobsterin_key_dict, calculation_type="onlycoop", validator_group="no_validator")
    wf = Workflow([fw], name=fw.name)
    wf = add_modify_incar(wf)
    wf = set_execution_options(wf, category=category, fworker_name="db1")
    wf = preserve_fworker(wf)
    print(wf)
    return wf

def transition_dipole_moment(input_cdft_df):
    import ast
    df = input_cdft_df.copy()
    df = df.fillna("None")
    src_dft = df.loc[df["task_label"] == "CDFT-B-HSE_scf", ["gs_taskid", "up_TDM_input", "dn_TDM_input", "task_id"]]
    db = get_db("HSE_triplets_from_Scan2dDefect", "calc_data-pbe_pc",  user="Jeng_ro", password="qimin", port=12347)
    tgt_db = get_db("HSE_triplets_from_Scan2dDefect", "cdft-pbe_pc",  user="Jeng", password="qimin", port=12347)
    for idx, row in src_dft.iterrows():
        taskid = row["task_id"]
        print(taskid)
        gs_taskid = row["gs_taskid"]
        up_TDM_input = ast.literal_eval(row["up_TDM_input"]) if type(row["up_TDM_input"]) == str else row["up_TDM_input"]
        dn_TDM_input = ast.literal_eval(row["dn_TDM_input"]) if type(row["dn_TDM_input"]) == str else row["dn_TDM_input"]
        up_TDM_input = up_TDM_input if up_TDM_input != ((),) else None
        dn_TDM_input = dn_TDM_input if dn_TDM_input != ((),) else None
        print(up_TDM_input, dn_TDM_input)
        up_TDM_sq = []
        dn_TDM_sq = []
        up_atomic_rate, dn_atomic_rate = [], []

        entry = db.collection.find_one({"task_id": gs_taskid})
        base_dir = "/home/qimin/sdc_tsai/Research/projects/HSE_triplets_from_Scan2dDefect/calc_data-pbe_pc"
        dir_name = os.path.join(base_dir, entry["dir_name"].split("/")[-1])

        def execute_vaspwfc():
            from vaspwfc import vaspwfc
            print(dir_name)
            os.chdir(dir_name)
            os.makedirs("TDM", exist_ok=True)
            shutil.copyfile("WAVECAR.gz", "TDM/WAVECAR.gz")
            try:
                x = subprocess.check_output(["gunzip", "TDM/WAVECAR.gz"], input=b"Y")
            except Exception:
                pass
            wv = vaspwfc("TDM/WAVECAR")
            shutil.rmtree("TDM")
            return wv

        wv = execute_vaspwfc()
        refractive_index = 1
        if up_TDM_input:
            for up_tdm_input in up_TDM_input:
                tdm = wv.TransitionDipoleMoment((1, 1, up_tdm_input[0]), (1, 1, up_tdm_input[1]))
                tdm_vec = tdm[-1]
                tdm_vec_norm_squared = [np.linalg.norm(component)**2*DEBYETOSI**2 for component in tdm_vec]
                dE = tdm[2]*EVTOJ
                up_TDM_sq.append(tdm_vec_norm_squared)
            up_TDM_sq = np.sum(up_TDM_sq, axis=0)
            up_atomic_rate = refractive_index * dE**3 * up_TDM_sq / 3 / PI / _eps0 / _c**3 / (_hplanck/(2*PI))**4
            up_atomic_rate *= 1e-6 #MHz
            up_atomic_rate = np.array(up_atomic_rate).round(3)

        if dn_TDM_input:
            for dn_tdm_input in dn_TDM_input:
                print(dn_tdm_input)
                tdm = wv.TransitionDipoleMoment((2, 1, dn_tdm_input[0]), (2, 1, dn_tdm_input[1]), )
                tdm_vec = tdm[-1]
                tdm_vec_norm_squared_in_debye = [np.linalg.norm(component)**2 for component in tdm_vec]
                print(tdm_vec_norm_squared_in_debye)
                tdm_vec_norm_squared = [np.linalg.norm(component)**2*DEBYETOSI**2 for component in tdm_vec]
                dE = tdm[2]*EVTOJ
                dn_TDM_sq.append(tdm_vec_norm_squared)
            dn_TDM_sq = np.sum(dn_TDM_sq, axis=0)
            dn_atomic_rate =  refractive_index * dE**3 * dn_TDM_sq / 3 / PI / _eps0 / _c**3 / (_hplanck/(2*PI))**4
            dn_atomic_rate *= 1e-6
            dn_atomic_rate = np.array(dn_atomic_rate).round(3)
        print("up_atomic_rate: {}, dn_atomic_rate: {}".format(up_atomic_rate, dn_atomic_rate))
        tgt_db.collection.update_one({"task_id": taskid},
                                     {"$set":
                                          {
                                              "TDM_transition": {
                                                  "up_transition_band": up_TDM_input,
                                                  "dn_transition_band": dn_TDM_input,
                                                  "up_TDM_squared": list(up_TDM_sq),
                                                  "dn_TDM_squared": list(dn_TDM_sq),
                                                  "KS_dE": dE,
                                                  "formula": "Einstein's A coeff.",
                                                  "unit": "MHz",
                                                  "up_KS_TDM_rate": list(up_atomic_rate),
                                                  "dn_KS_TDM_rate": list(dn_atomic_rate)
                                              }}})

def finite_size_effect_HSE_wf(distort=0.0, category="finite_size_effect", pyzfs_fw=True, irvsp_fw=True): #1e-4
    wfs = []
    # db_name, col_name = "HSE_triplets_from_Scan2dDefect", "hse_pc"
    db_name, col_name = "2dMat_from_cmr_fysik", "2dMaterial_v1"
    hse_pc_db = get_db(db_name, col_name,  port=12345, user="readUser", password="qiminyan")
    pc_db = hse_pc_db

    qubit_candidate_df = IOTools(cwd=os.path.join("~/scripts/JPack_independent/projects/defectDB", INPUT_PATH),
                                 excel_file="qubit_2021-11-18").read_excel()
    triplet_df = qubit_candidate_df.copy()
    triplet_df = triplet_df.loc[triplet_df["uid"].isin(["BN-BN-NM"])]

    # defect_model_inputs = zip(triplet_df["uid"],
    #     triplet_df["charge"],
    #     triplet_df["defect_name"],
    #     triplet_df["defect_type"],
    #     triplet_df["task_id"]
    #     )
    defect_model_inputs = zip(
        ["BN-BN-NM"],
        [-1],
        ["vac_1_B"],
        ["vacancy"],
        [None]
    )

    for uid, charge, defect_name, defect_type, task_id in defect_model_inputs:

        mx2 = pc_db.collection.find_one({"uid": uid})
        pc = Structure.from_dict(mx2["structure"])
        scaling = find_scaling_for_2d_defect(pc, 15)[0] # determine sc size
        area = scaling[0]*scaling[1]
        geo_spec = {area*pc.num_sites: [20]} # vacuum

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
                        sub_on_side=None
                    )

                    wf = get_wf_full_hse(
                        structure=gen_defect.defect_st,
                        task_arg=dict(lcharg=True),
                        charge_states=[charge],
                        gamma_only=False,
                        gamma_mesh=True,
                        nupdowns=[2],
                        task="hse_scf",
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
                                                   irvsptodb_kwargs=dict(collection_name="finite_size_effect"),
                                                   name="irvsp",
                                                   parents=fw)
                                HSE_fws.append(irvsp_fw)
                        wf = Workflow(HSE_fws, name=wf.name)

                    if pyzfs_fw:
                        HSE_fws = wf.fws.copy()
                        for idx, fw in enumerate(HSE_fws):
                            if "HSE_scf" in fw.name:
                                modify_incar_task_idx = None
                                for task_idx, task in enumerate(fw.tasks):
                                    if "ModifyIncar" in task.fw_name:
                                        modify_incar_task_idx = task_idx
                                if fw.tasks[modify_incar_task_idx]["incar_update"]["NUPDOWN"] == 2:
                                    pyzfs_fw = PyzfsFW(structure=gen_defect.defect_st, parents=fw, pyzfstodb_kwargs=dict(
                                    collection_name="finite_size_effect"))
                                HSE_fws.append(pyzfs_fw)
                        wf = Workflow(HSE_fws, name=wf.name)

                    wf = add_tags(wf, [{"category": category, "uid": mx2["uid"],}])
                    wf = bash_scp_files(
                        wf,
                        dest="/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/finite_size_effect/",
                        port=12348,
                        task_name_constraint="RunVasp"
                    )
                    wf = set_execution_options(wf, category=category, fw_name_constraint="HSE_relax")
                    wf = set_execution_options(wf, category=category, fw_name_constraint="HSE_scf")

                    if irvsp_fw:
                        wf = set_queue_options(wf, walltime="06:00:00", fw_name_constraint="irvsp")
                        wf = set_execution_options(wf, category=category, fw_name_constraint="irvsp")
                        wf = bash_scp_files(
                            wf,
                            dest="/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/finite_size_effect/",
                            port=12348,
                            task_name_constraint="ToDb",
                            fw_name_constraint="irvsp"
                        )
                        wf = clean_up_files(wf, files=["WAVECAR*", "vasprun.xml*"], fw_name_constraint="irvsp",
                                            task_name_constraint="ToDb")


                    if pyzfs_fw:
                        wf = set_queue_options(wf, walltime="06:00:00", fw_name_constraint="pyzfs")
                        wf = set_execution_options(wf, category=category, fw_name_constraint="pyzfs")
                        wf = bash_scp_files(
                            wf,
                            dest="/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/finite_size_effect/",
                            port=12348,
                            task_name_constraint="ToDb",
                            fw_name_constraint="pyzfs"
                        )
                        wf = clean_up_files(wf, files=["WAVECAR*", "vasprun.xml*"], fw_name_constraint="pyzfs",
                                            task_name_constraint="ToDb")
                    wf = add_modify_incar(wf)
                    wf.name = "FSE:dist{}:{}:{}:{}".format(distort, mx2["uid"],
                                                           defect_name, ":".join(wf.name.split(":")[3:5]))
                    print(wf.name)
                    wfs.append(wf)
    return wfs


def scan_pc_in_scf_way():
    cat = "scan_pc_in_scf_way"
    lpad = LaunchPad.from_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/Scan2dMat/{}/my_launchpad.yaml".format(cat)))


    col = get_db("Scan2dMat", "calc_data", port=12347, user="Jeng_ro", password="qimin").collection

    mx2s = list(col.find({"task_label": "SCAN_nscf line"}))
    for idx, mx2 in enumerate(mx2s[:1]): #start from 100:!myq
        pc = Structure.from_dict(mx2["output"]["structure"])
        from pathlib import Path
        MODULE_DIR = Path("__file__").resolve().parent
        wf = get_wf(pc, os.path.join(MODULE_DIR, "wf/Scan2dMat/scan_pc_in_scf_way.yaml"))

        wf = add_additional_fields_to_taskdocs(
            wf,
            {
                "c2db_info": mx2["c2db_info"],
                "site_info": mx2["site_info"],
                "sym_data": mx2["sym_data"],
                "band_edges": mx2["band_edges"],
                "site_oxi_state": mx2["site_oxi_state"],
                "taskid-calc_data-Scan2dMat": mx2["task_id"],
                "bandgap-calc_data-Scan2dMat": mx2["output"]["bandgap"]
            },
            task_name_constraint="ToDb"
        )

        kpt = Kpoints.automatic_density_by_vol(pc, 144)
        wf = add_modify_kpoints(
            wf,
            modify_kpoints_params={"kpoints_update": {"kpts": kpt.kpts}},
            fw_name_constraint=wf.fws[0].name
        )
        wf = add_modify_2d_nscf_kpoints(
            wf,
            is_hse=True,
            modify_kpoints_params={"kpoints_line_density": 35, "reciprocal_density": 144},
            fw_name_constraint=wf.fws[1].name
        )
        wf = add_modify_2d_nscf_kpoints(
            wf,
            is_hse=True,
            modify_kpoints_params={"kpoints_line_density": 35, "reciprocal_density": 144, "mode": "uniform"},
            fw_name_constraint=wf.fws[2].name
        )

        if idx % 2 == 1:
            wf = set_execution_options(wf, category=cat, fworker_name="efrc")
        else:
            wf = set_execution_options(wf, category=cat, fworker_name="owls")
        wf = set_queue_options(wf, "04:00:00", fw_name_constraint=wf.fws[0].name)
        wf = set_queue_options(wf, "02:00:00", fw_name_constraint=wf.fws[1].name)
        wf = set_queue_options(wf, "01:00:00", fw_name_constraint=wf.fws[2].name)
        wf = set_queue_options(wf, "01:00:00", fw_name_constraint=wf.fws[3].name)

        wf = clean_up_files(wf, files=["WAVECAR"], task_name_constraint="IRVSPToDb", fw_name_constraint="irvsp")
        wf = clean_up_files(wf, files=["CHGCAR*"], fw_name_constraint=wf.fws[1].name)
        wf = clean_up_files(wf, files=["CHGCAR*", "WAVECAR*"], fw_name_constraint=wf.fws[2].name)

        wf = preserve_fworker(wf)
        wf = add_modify_incar(wf)
        wf.name = "{}:r2scan".format(mx2["formula_pretty"])
        lpad.add_wf(wf)


def pyzfs_wf():
    wfs = []
    for scf_taskid in [1088]:
        fws = []
        e = HSEQubitDefect.collection.find_one({"task_id": scf_taskid})
        scf_dir = e["dir_name"].split(":")[-1]
        st = Structure.from_dict(e["output"]["structure"])
        pyzfs_fw = PyzfsFW(structure=st, prev_calc_dir=scf_dir, pyzfstodb_kwargs=dict(
            collection_name="zfs_data-pbe_pc"))
        fws.append(pyzfs_fw)

        wf = Workflow(fws, name="taskid{}:{}:{}".format(scf_taskid, e["formula_pretty"], e["pc_from"].split("/")[-1]))
        wf = add_additional_fields_to_taskdocs(wf, {"prev_fw_taskid": scf_taskid, "prev_fw_db":
            HSEQubitDefect.db_name, "prev_fw_collection": HSEQubitDefect.collection.name}, task_name_constraint="ToDb")
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
        wf = add_modify_incar(wf)
        wfs.append(wf)
    return wfs


def transition_level_wf(category="charge_state"):
    from unfold import find_K_from_k
    # for analysis and plotting, use script in qubitPak/defect_formation_energy_correction

    aexx = 0.25

    db_name, col_name = "2dMat_from_cmr_fysik", "2dMaterial_v1"
    c2db = get_db(db_name, col_name,  port=12345, user="readUser", password="qiminyan")


    host_uids = [
        "In2S2-GaSe-NM", "In2Se2-GaSe-NM", "In2Te2-GaSe-NM",
        "Ga2S2-GaSe-NM", "Ga2Se2-GaSe-NM", "Ga2Te2-GaSe-NM",
        "In2S2-GaS-NM", "In2Se2-GaS-NM", "In2Te2-GaS-NM",
        "Ga2S2-GaS-NM", "Ga2Se2-GaS-NM", "Ga2Te2-GaS-NM"
    ]
    mx2s = c2db.collection.find({"uid": {"$in": host_uids}})

    geo_spec = {5*5*4: [25], 4*4*4: [20, 25]}
    wfs = []
    for mx2 in mx2s:
        pc = Structure.from_dict(mx2["structure"])
        defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
        cation, anion = find_cation_anion(pc)
        defect_type = "substitutions"
        for idx, antisite in enumerate(range(len(defect[defect_type]))):
            print(cation, anion)
            if "{}_on_{}".format(anion, cation) not in defect["substitutions"][antisite]["name"]:
                continue
            for na, thicks in geo_spec.items():
                sc_size = np.eye(3, dtype=int)*math.sqrt(na/3)

                for thick in thicks:
                    # ++++++++++++++++++++++++++++++ defect formation energy+++++++++++++++++++++++++++++++
                    def defect_wf():
                        antisite_st = GenDefect(
                            orig_st=pc,
                            defect_type=(defect_type, antisite),
                            natom=na,
                            vacuum_thickness=thick,
                            sub_on_side=None,
                        )

                        wf_antisite = get_wf_full_hse(
                            structure=antisite_st.defect_st,
                            charge_states=[-2, -3],
                            gamma_only=False,#[list(find_K_from_k(sc_size, [0, 0, 0])[0])],
                            gamma_mesh=True,
                            nupdowns=[1, 0],
                            task="hse_relax-hse_scf",
                            vasptodb={
                                "category": category,
                                "NN": antisite_st.NN,
                                "defect_entry": antisite_st.defect_entry,
                                "defect_type": "defect",
                                "pmg_obj": antisite_st.pmg_obj.as_dict(),
                                "c2db_uid": mx2["uid"]
                            },
                            wf_addition_name="{}:{}".format(na, thick),
                            task_arg={"parse_dos":True}
                        )

                        wf_antisite = bash_scp_files(
                            wf_antisite,
                            dest="/mnt/sdc/tsai/Research/projects/defect_qubit_in_36_group/charge_state/",
                            port=12348,
                            task_name_constraint="ToDb",
                            fw_name_constraint="HSE_scf"
                        )

                        wf_antisite = add_modify_incar(
                            wf_antisite,
                            {"incar_update":{"AEXX": aexx, "LWAVE": False}},
                        )
                        wf_antisite = add_modify_incar(wf_antisite)
                        wf_antisite = set_queue_options(wf_antisite, "06:00:00", fw_name_constraint="HSE_scf")
                        wf_antisite = set_queue_options(wf_antisite, "06:00:00", fw_name_constraint="HSE_relax")
                        # related to directory
                        wf_antisite = set_execution_options(wf_antisite, category=category)
                        wf_antisite = preserve_fworker(wf_antisite)
                        wfs.append(wf_antisite)

                    # +++++++++++++++++++++bulk formation energy part+++++++++++++++++++++++++++++++++++
                    def bulk_wf():
                        antisite_st = GenDefect(
                            orig_st=pc,
                            defect_type=("bulk", antisite),
                            natom=na,
                            vacuum_thickness=thick,
                            sub_on_side=None,
                        )

                        wf_bulk = get_wf_full_hse(
                            structure=antisite_st.bulk_st,
                            charge_states=[0],
                            gamma_only=False, #[list(find_K_from_k(sc_size, [0, 0, 0])[0])],
                            gamma_mesh=True,
                            nupdowns=[-1],
                            task="hse_relax-hse_scf",
                            vasptodb={"category": category, "defect_type": "host", "c2db_uid": mx2["uid"]},
                            wf_addition_name="{}:{}".format(na, thick),
                            task_arg={"parse_dos":False}
                        )

                        wf_bulk = add_modify_incar(
                            wf_bulk,
                            {"incar_update":{"AEXX": aexx, "LWAVE": False}},
                        )
                        wf_bulk = add_modify_incar(wf_bulk)
                        wf_bulk = set_queue_options(wf_bulk, "24:00:00", fw_name_constraint="HSE_scf")
                        wf_bulk = set_queue_options(wf_bulk, "24:00:00", fw_name_constraint="HSE_relax")
                        wf_bulk = preserve_fworker(wf_bulk)
                        # wfs.append(wf_bulk)

                    #+++++++++++++++++++++ bulk VBM++++++++++++++++++++++++++++++++++++++++++++++++++++
                    def vbm_wf():
                        antisite_st = GenDefect(
                            orig_st=pc,
                            defect_type=("bulk", antisite),
                            natom=na,
                            vacuum_thickness=thick,
                            sub_on_side=None,
                        )

                        wf_bulk_vbm = get_wf_full_hse(
                            structure=antisite_st.bulk_st,
                            charge_states=[0],
                            gamma_only=[list(find_K_from_k(sc_size, [0, 0, 0])[0])],
                            gamma_mesh=False,
                            nupdowns=[-1],
                            task="hse_relax-hse_scf",
                            vasptodb={"category": category, "defect_type": "vbm", "c2db_uid": mx2["uid"]},
                            wf_addition_name="{}:{}".format(na, thick),
                            task_arg={"parse_dos":False}
                        )

                        wf_bulk_vbm = add_modify_incar(
                            wf_bulk_vbm,
                            {"incar_update":{"AEXX": aexx,"LWAVE": False}},
                        )
                        wf_bulk_vbm = add_modify_incar(wf_bulk_vbm)
                        wf_bulk_vbm = set_queue_options(wf_bulk_vbm, "24:00:00", fw_name_constraint="HSE_scf")
                        wf_bulk_vbm = set_queue_options(wf_bulk_vbm, "24:00:00", fw_name_constraint="HSE_relax")
                        wf_bulk_vbm = preserve_fworker(wf_bulk_vbm)
                        wfs.append(wf_bulk_vbm)

                    defect_wf()
                    # bulk_wf()
                    # vbm_wf()
    return wfs




def main():
    def run_grand_hse_defect():
        lpad = LaunchPad.from_file(
            os.path.join(
                os.path.expanduser("~"),
                "config/project/C2DB_IR_vacancy_HSE/calc_data/my_launchpad.yaml"
            )
        )
        wfs = grand_hse_defect()
        for idx, wf in enumerate(wfs):
            lpad.add_wf(wf)
            print(idx, wf.name, wf.fws[0].tasks[-1]["additional_fields"]["pc_from_id"])

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
        for idx, wf in enumerate(wfs[:1]):
            wf = set_execution_options(wf, fworker_name="gpu_nersc")
            wf = preserve_fworker(wf)

            # wte2_magmom = [0.006, 0.0, 0.001, -0.0, -0.001, 0.0, -0.002, -0.0, -0.001,
            #                0.001, -0.002, 0.0, -0.0, -0.0, 0.0, -0.0, -0.001, 0.0, -0.002,
            #                0.006, 0.001, -0.0, -0.002, 0.001, 0.056, -0.0, -0.0, -0.001,
            #                0.005, 0.005, 0.006, -0.002, 0.0, -0.003, -0.001, -0.002,  -0.0, 0.0, 0.0, -0.0, -0.001,
            #                -0.001, -0.0, -0.002, -0.0, 0.001, -0.0, -0.002, 0.006,  -0.135, -0.004, -0.004, -0.012,
            #                0.192, -0.213, -0.142, 0.046, 0.004, 0.047, 0.045, 0.207,  0.018, -0.002, -0.004, 0.003,
            #                -0.012, -0.001, -0.002, -0.002, 0.044, -0.003, 0.014, -0.001, 0.018, 1.224]
            # mote2_magmom = [-0.001, 0.003, 0.0, 0.004, 0.002, 0.008, -0.003, 0.001, -0.005, 0.0,  -0.001, 0.001,
            #                 0.001, 0.001, 0.002, -0.0, 0.001, 0.0, -0.003, -0.001, -0.007,  -0.0, -0.0, 0.009, 0.075,
            #                 0.007, 0.002, 0.0, -0.003, -0.003, 0.007, -0.003, 0.001, -0.003,  0.0, -0.001, -0.0,
            #                 0.001, 0.001, 0.002, -0.001, 0.001, -0.0, -0.003, 0.006, -0.007,  -0.0, -0.0, 0.006,
            #                 -0.427, 0.427, -0.032, -0.005, 0.093, -0.218, -0.446, 0.024,  -0.012, 0.031, 0.026,
            #                 0.135, 0.067, -0.002, -0.001, -0.012, -0.008, -0.021, -0.004,  -0.002, 0.019, -0.028,
            #                 0.009, -0.017, 0.07, 1.843]



            # wf = add_modify_incar(wf, {"incar_update": {"MAGMOM": wte2_magmom}})
            lpad.add_wf(wf)
            print(idx, wf.name)

    def run_cdft():
        # qubit_df = IOTools(excel_file="hse_screened_qubit_2021-11-18", cwd=INPUT_PATH).read_excel()
        qubit_df = IOTools(excel_file="pure_c2db_TMD_2022-07-08", cwd=INPUT_PATH).read_excel()
        qubit_df = qubit_df[qubit_df["task_id"].isin([2739])]
        print(qubit_df.head())

        fworker = "gpu_nersc"
        taskids = qubit_df.loc[:, "task_id"]
        print(f"taskids: {taskids}")
        for taskid in taskids[:]:
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

    def run_cohp():
        a = Defect(defect_xlsx="defect_2021-11-18")
        df = a.defect_df.copy()
        pair_taskids =(172, 874)
        df = df.loc[df["task_id"].isin(pair_taskids)]

        db = get_db("Scan2dDefect", "calc_data",  user="Jeng_ro", password="qimin", port=12347)
        # db = get_db("HSE_triplets_from_Scan2dDefect", "calc_data-pbe_pc",  user="Jeng_ro", password="qimin", port=12347)
        for taskid in df["task_id"]:
            level_cat = df.loc[df["task_id"] == taskid, "level_cat"].iloc[0]
            e = db.collection.find_one({"task_id": taskid})
            dir_name = e["dir_name"].split("/")[-1]

            # dir_path = os.path.join("/home/qimin/sdc_tsai/Research/projects/HSE_triplets_from_Scan2dDefect/"
            #                         "calc_data-pbe_pc", dir_name)
            dir_path = os.path.join("/home/qimin/sdc_tsai/Research/projects/Scan2dDefect/calc_data/scf", dir_name)

            st = Structure.from_dict(e["output"]["structure"])
            element = st.sites[e["NN"][0]].specie
            wf = cohp(Structure.from_dict(e["output"]["structure"]), scf_path=dir_path, lobsterin_key_dict={
                "COHPstartEnergy": -10, "COHPendEnergy": 10,
                "cohpGenerator": "from 2 to 6 orbitalwise"}
                # "cohpGenerator": "from 0.1 to 5.0 type {} type {} orbitalwise".format(element, element)}
                      )

            # wf = add_additional_fields_to_taskdocs(
            #     wf,
            #     {
            #         "prev_fw_db": "HSE_triplets_from_Scan2dDefect", "prev_fw_collection": "calc_data-pbe_pc",
            #         "prev_fw_taskid": taskid,
            #         "NN": e["NN"], "pc_from": e["pc_from"], "charge_state": e["charge_state"], "defect_name": e[
            #         "defect_entry"]["name"]
            #     }, task_name_constraint="LobsterRunToDb")
            wf = add_additional_fields_to_taskdocs(
                wf,
                {
                    "prev_fw_db": "Scan2dDefect", "prev_fw_collection": "calc_data",
                    "prev_fw_taskid": taskid,
                    "NN": e["NN"], "pc_from": e["pc_from"], "charge_state": e["charge_state"], "defect_name": e[
                    "defect_entry"]["name"], "level_cat": level_cat, "pair": pair_taskids
                }, task_name_constraint="LobsterRunToDb")

            wf = clean_up_files(wf, files=["WAVECAR*", "vasprun.xml*", "CHG*"], task_name_constraint="ToDb")
            # lpad = LaunchPad.from_file(
            #     os.path.join("/home/qimin/sdb_tsai/config/project/HSE_triplets_from_Scan2dDefect/cohp",
            #                  "my_launchpad.yaml")
            # )

            wf.name += "-{}".format(level_cat)
            lpad = LaunchPad.from_file(
                os.path.join("/home/qimin/sdb_tsai/config/project/Scan2dDefect/cohp",
                             "my_launchpad.yaml")
            )
            lpad.add_wf(wf)

    def run_TDM():
        p_path = "/home/qimin/sdb_tsai/site-packages/JPack_independent/projects/defectDB"
        tgt_df = IOTools(excel_file="wte2_mote2", cwd=os.path.join(p_path, INPUT_PATH)).read_excel()
        transition_dipole_moment(tgt_df)


    def run_finite_size_effect_wf():
        def std_wf():
            lpad = LaunchPad.from_file(
                os.path.join(
                    os.path.expanduser("~"),
                    "config/project/HSE_triplets_from_Scan2dDefect/finite_size_effect/my_launchpad.yaml"))
            wfs = finite_size_effect_HSE_wf(distort=0.1, irvsp_fw=False)
            for idx, wf in enumerate(wfs[:]):
                wf = set_execution_options(wf, fworker_name="efrc")
                wf = preserve_fworker(wf)
                wf = add_modify_incar(wf, {"incar_update": {"ENCUT": 910, "AEXX": 0.32}})
                wf = add_modify_potcar(wf, {"potcar_symbols": {"B": "B_h", "N":"N_h"}})
                wf = add_modify_incar(wf)
                lpad.add_wf(wf)
                print(wf.name)

        def cdft_wf():
            taskid = 2518
            fworker = "owls"
            print("{}=".format(taskid)*20)
            nbands = 360
            up_occ = "None"
            dn_occ = "253*1 1*0 1*1 1*0 104*0"
            wf, cat = cdft(
                "B-C-D",
                taskid,
                occ={"1": up_occ, "-1": dn_occ},
                nbands=nbands,
                category="finite_size_effect",
            )
            wf = add_modify_incar(wf, {"incar_update": {"AEXX": 0.32}})
            wf = use_custodian(wf, custodian_params={"vasp_cmd": ">>vasp5_std<<", "handler_group":[]})
            wf = add_additional_fields_to_taskdocs(wf, {"cat":cat, "taskid": taskid})
            wf = set_execution_options(wf, fworker_name=fworker)
            lpad = LaunchPad.from_file(
                os.path.join(
                    os.path.expanduser("~"),
                    "config/project/HSE_triplets_from_Scan2dDefect/finite_size_effect/my_launchpad.yaml"))
            lpad.add_wf(wf)
        std_wf()

    def run_scan_pc_in_scf_way():
        scan_pc_in_scf_way()

    def run_pyzfs_wf():
        lpad = LaunchPad.from_file(
            os.path.join(
                os.path.expanduser("~"),
                "config/project/HSE_triplets_from_Scan2dDefect/calc_data/my_launchpad.yaml"))
        wfs = pyzfs_wf()
        for idx, wf in enumerate(wfs[:]):
            wf = set_execution_options(wf, fworker_name="owls")
            wf = preserve_fworker(wf)
            lpad.add_wf(wf)
            print(idx, wf.name)

    def run_transition_levels():
        cat = "charge_state"
        lpad = LaunchPad.from_file("/home/tug03990/config/project/defect_qubit_in_36_group/{}/"
                                   "my_launchpad.yaml".format(cat))
        wfs = transition_level_wf()
        print(*wfs)
        for wf in wfs:
            wf = set_execution_options(wf, category=cat, fworker_name="gpu_nersc")
            # add wf tags
            lpad.add_wf(wf)
        print(len(wfs))

    run_grand_hse_defect()


if __name__ == '__main__':
    main()




