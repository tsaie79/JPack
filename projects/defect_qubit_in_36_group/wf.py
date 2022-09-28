import json
import shutil
import subprocess

import monty.serialization
from monty.subprocess import Popen

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
import pandas as pd

from JPack_independent.projects.antisiteQubit.wf_defect import ZPLWF
from analysis.data_analysis import Defect

INPUT_PATH = "analysis/input"
OUTPUT_PATH = "analysis/output/xlsx"



def hse_wf(distort=0.0, category="hse_calc_data", pyzfs_fw=True, irvsp_fw=True): #1e-4
    wfs = []
    db_name, col_name = "2dMat_from_cmr_fysik", "2dMaterial_v1"
    hse_pc_db = get_db(db_name, col_name,  port=12345, user="readUser", password="qiminyan")
    pc_db = hse_pc_db

    uid = []
    charge = []
    defect_name = []
    defect_type = []

    for uid, charge, defect_name, defect_type, task_id in zip(uid, charge, defect_name, defect_type):
        mx2 = pc_db.collection.find_one({"uid": uid})
        pc = Structure.from_dict(mx2["structure"])
        scaling = find_scaling_for_2d_defect(pc, 20)[0]
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
                        dest="/mnt/sdc/tsai/Research/projects/defect_qubit_in_36_group/hse_calc_data/",
                        port=12348,
                        task_name_constraint="RunVasp"
                    )
                    wf = set_execution_options(wf, category=category, fw_name_constraint="HSE_relax")
                    wf = set_execution_options(wf, category=category, fw_name_constraint="HSE_scf")

                    if irvsp_fw:
                        wf = set_execution_options(wf, category="ir_data", fw_name_constraint="irvsp")
                        wf = bash_scp_files(
                            wf,
                            dest="/mnt/sdc/tsai/Research/projects/defect_qubit_in_36_group/hse_calc_data/",
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
                            dest="/mnt/sdc/tsai/Research/projects/defect_qubit_in_36_group/hse_calc_data/",
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
    std_cat_name = "finite_size_effect"
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
        dest="/mnt/sdc/tsai/Research/projects/HSE_triplets_from_Scan2dDefect/{}".format(std_cat_name),
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


def transition_level_wf(cls, category="charge_state"):
    # for analysis and plotting, use script in qubitPak/defect_formation_energy_correction
    lpad = LaunchPad.from_file("/home/tug03990/config/project/defect_qubit_in_36_group/{}/"
                               "my_launchpad.yaml".format(category))
    aexx = 0.25

    db_name, col_name = "2dMat_from_cmr_fysik", "2dMaterial_v1"
    c2db = get_db(db_name, col_name,  port=12345, user="readUser", password="qiminyan")

    mx2s = c2db.collection.find({"uid": {"$in": ["Ga2S2-GaS-NM"]}})
    geo_spec = {5*5*4: [20]}
    for mx2 in mx2s:
        # pc = Structure.from_dict(mx2["structure"])
        pc = Structure.from_dict(mx2["output"]["structure"])
        defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
        cation, anion = find_cation_anion(pc)
        defect_type = "substitutions"
        for idx, antisite in enumerate(range(len(defect[defect_type]))):
            print(cation, anion)
            # if "{}_on_{}".format(cation, anion) not in defect["substitutions"][antisite]["name"]:
            #     continue
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
                            charge_states=[0],
                            gamma_only=False,#[list(find_K_from_k(sc_size, [0, 0, 0])[0])],
                            gamma_mesh=True,
                            nupdowns=[-1],
                            task="hse_relax-hse_scf",
                            vasptodb={
                                "category": category,
                                "NN": antisite_st.NN,
                                "defect_entry": antisite_st.defect_entry,
                                "defect_type": "defect",
                                "pmg_obj": antisite_st.pmg_obj.as_dict()
                            },
                            wf_addition_name="{}:{}".format(na, thick),
                            task_arg={"parse_dos":False}
                        )

                        wf_antisite = add_modify_incar(
                            wf_antisite,
                            {"incar_update":{"AEXX": aexx, "LWAVE": False}},
                        )
                        wf_antisite = add_modify_incar(wf_antisite)
                        wf_antisite = set_queue_options(wf_antisite, "24:00:00", fw_name_constraint="HSE_scf")
                        wf_antisite = set_queue_options(wf_antisite, "24:00:00", fw_name_constraint="HSE_relax")
                        # related to directory
                        wf_antisite = set_execution_options(wf_antisite, category=category)
                        wf_antisite = preserve_fworker(wf_antisite)
                        return wf_antisite
                    # lpad.add_wf(defect_wf())

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
                            vasptodb={"category": category, "defect_type": "host"},
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
                        wf_bulk = set_execution_options(wf_bulk, category=category, fworker_name="efrc")
                        wf_bulk = preserve_fworker(wf_bulk)
                        return wf_bulk
                    lpad.add_wf(bulk_wf())

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
                            gamma_only=[list(find_K_from_k(sc_size, [0, 1/3, 0])[0])],
                            gamma_mesh=False,
                            nupdowns=[-1],
                            task="hse_relax-hse_scf",
                            vasptodb={"category": category, "defect_type": "vbm"},
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
                        wf_bulk_vbm = set_execution_options(wf_bulk_vbm, category=category)
                        wf_bulk_vbm = preserve_fworker(wf_bulk_vbm)
                        return wf_bulk_vbm
                    # lpad.add_wf(vbm_wf())


