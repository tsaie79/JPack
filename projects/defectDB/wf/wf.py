from fireworks import LaunchPad, Workflow

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
from atomate.vasp.workflows.base.core import get_wf

from my_atomate_jyt.vasp.powerups import *
from my_atomate_jyt.vasp.workflows.wf_full import get_wf_full_scan

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet


from qubitPack.tool_box import *

from mpinterfaces.utils import *

from monty.serialization import loadfn
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
    cat = "test"
    lpad = LaunchPad.from_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/test/{}/my_launchpad.yaml".format(cat)))

    col = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12345, user="readUser", password="qiminyan").collection

    # mx2s = col.find(
    #     {
    #         "gap_hse_nosoc":{"$gt": 0},
    #         "nkinds": 2,
    #         "magstate": "NM",
    #         "nsites": 12,
    #     }
    # )

    mx2s = col.find(
        {
            "uid": "C2-C-NM"
        }
    )



    for idx, mx2 in enumerate(mx2s):
        pc = Structure.from_dict(mx2["structure"])
        print(pc.num_sites)
        pc = ensure_vacuum(pc, 30)
        pc, site_info = phonopy_structure(pc)
        print(pc.num_sites)
        wf = get_wf(pc, os.path.join(os.path.dirname(os.path.abspath("__file__")), "projects/defectDB/wf/scan_pc.yaml"))


        wf = add_additional_fields_to_taskdocs(
            wf,
            {
                "c2db_uid": mx2["uid"],
                "c2db_nkinds": mx2["nkinds"],
                "site_info": site_info
            },
            task_name_constraint="ToDb"
        )

        kpt = Kpoints.automatic_density_by_vol(pc, 144)
        wf = add_modify_kpoints(
            wf,
            modify_kpoints_params={"kpoints_update": {"kpts": kpt.kpts}},
            fw_name_constraint=wf.fws[2].name
        )
        wf = add_modify_2d_nscf_kpoints(
            wf,
            modify_kpoints_params={"kpoints_line_density": 20, "reciprocal_density": 144},
            fw_name_constraint=wf.fws[3].name
        )
        wf = add_modify_2d_nscf_kpoints(
            wf,
            modify_kpoints_params={"kpoints_line_density": 20, "reciprocal_density": 144, "mode": "uniform"},
            fw_name_constraint=wf.fws[4].name
        )

        if idx % 2 == 1:
            wf = set_execution_options(wf, category=cat, fworker_name="owls")
        else:
            wf = set_execution_options(wf, category=cat, fworker_name="nersc")
        wf = set_queue_options(wf, "00:30:00", fw_name_constraint=wf.fws[0].name)
        wf = set_queue_options(wf, "03:00:00", fw_name_constraint=wf.fws[1].name)
        wf = set_queue_options(wf, "03:00:00", fw_name_constraint=wf.fws[2].name)
        wf = set_queue_options(wf, "02:00:00", fw_name_constraint=wf.fws[3].name)
        wf = set_queue_options(wf, "02:00:00", fw_name_constraint=wf.fws[4].name)
        wf = set_queue_options(wf, "00:30:00", fw_name_constraint=wf.fws[5].name)

        wf = clean_up_files(wf, files=["WAVECAR"], task_name_constraint="IRVSPToDb", fw_name_constraint="irvsp")
        wf = clean_up_files(wf, files=["CHGCAR*"], fw_name_constraint=wf.fws[3].name)
        wf = clean_up_files(wf, files=["CHGCAR*"], fw_name_constraint=wf.fws[4].name)

        wf = add_modify_incar(wf, {"incar_update": {"METAGGA": "R2SCAN"}})
        wf = preserve_fworker(wf)
        wf = add_modify_incar(wf)
        wf = bash_scp_files(wf, dest="/home/tsai/Research/projects/Scan2dMat/calc_data", port=12348)
        wf.name = "{}:SCAN_full".format(mx2["uid"])
        lpad.add_wf(wf)

def binary_scan_defect(defect_choice="substitutions", impurity_on_nn=None): #BN_vac

    col = get_db("symBaseBinaryQubit", "scan_relax_pc", port=12345).collection

    mx2s = col.find(
        {
            "task_id": 161,
        }
    )
    geo_spec = None
    aexx = 0.25
    test = []
    for idx, mx2 in enumerate(mx2s[:1]): #mx2s[14:34]
        pc = Structure.from_dict(mx2["output"]["structure"])
        cat = None

        if defect_choice == "substitutions":
            cat = "binary_defect_sub_"
        elif defect_choice == "vacancies":
            cat = "bincary_defect_vac_"

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
        if mx2["nsites"] == 8:
            geo_spec = {area*8: [20]}
            cat += "AB3"

        lpad = LaunchPad.from_file(
            os.path.join(
                os.path.expanduser("~"),
                "config/project/Scan2dDefect/calc_data/{}/my_launchpad.yaml".format(cat)))
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
            # substitutions or vacancies!
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
                                    dos=False,
                                    nupdowns=nupdn,
                                    vasptodb={
                                        "category": cat,
                                        "NN": gen_defect.NN,
                                        "NN_dist": gen_defect.nn_dist,
                                        "defect_entry": defect_data,
                                    },
                                    wf_addition_name="{}:{}".format(gen_defect.defect_st.num_sites, thick),
                                    task="scan_relax-scan_scf",
                                    wf_yaml=os.path.join(os.path.dirname(os.path.abspath("__file__")), "projects/defectDB/wf/scan_defect.yaml")
                                )

                                wf = add_modify_2d_nscf_kpoints(
                                    wf,
                                    modify_kpoints_params={"mode": "line"},
                                    fw_name_constraint=wf.fws[2].name
                                )
                                wf = add_modify_2d_nscf_kpoints(
                                    wf,
                                    modify_kpoints_params={"mode": "uniform"},
                                    fw_name_constraint=wf.fws[3].name
                                )
                                wf = add_modify_incar(wf, {"incar_update": {"METAGGA":"SCAN"}})
                                wf = add_modify_incar(wf)
                                wf = add_additional_fields_to_taskdocs(
                                    wf,
                                    {"pc_from": "symBaseBinaryQubit/scan_relax_pc/SCAN_relax/{}".format(mx2["task_id"]),
                                     "pc_from_id": mx2["task_id"],
                                     "combo": mx2["sym_data"]["combo"],
                                     "class": mx2["c2db_info"]["class"],
                                     "defect_name": defect_data["name"],
                                     "site_info": gen_defect.site_info,
                                     "perturbed": gen_defect.distort}
                                )
                                wf = add_additional_fields_to_taskdocs(
                                    wf,
                                    {"pc_from": "symBaseBinaryQubit/scan_relax_pc/SCAN_relax/{}".format(mx2["task_id"]),
                                     "pc_from_id": mx2["task_id"],
                                     "class": mx2["c2db_info"]["class"],
                                     "defect_name": defect_data["name"],
                                     },
                                    task_name_constraint="IRVSPToDb"
                                )

                                if idx % 2 == 1:
                                    wf = set_execution_options(wf, category=cat, fworker_name="efrc")
                                else:
                                    wf = set_execution_options(wf, category=cat, fworker_name="owls")
                                wf = set_queue_options(wf, "24:00:00")
                                wf = set_queue_options(wf, "00:30:00", fw_name_constraint="irvsp")
                                wf = preserve_fworker(wf)
                                wf.name = wf.name+":dx[{}]".format(gen_defect.distort)
                                wf = bash_scp_files(
                                    wf,
                                    dest="/home/tsai/Research/projects/Scan2dDefect/calc_data/scf",
                                    port=12348,
                                    fw_name_constraint="SCAN_scf",
                                    task_name_constraint="RunVasp"
                                )
                                wf = clean_up_files(wf, files=["WAVECAR*"], task_name_constraint="SCP",
                                                    fw_name_constraint=wf.fws[1].name)
                                wf = clean_up_files(wf, files=["CHGCAR*"], task_name_constraint="ToDb",
                                                    fw_name_constraint=wf.fws[3].name)
                                wf = clean_up_files(wf, files=["WAVECAR"], task_name_constraint="IRVSPToDb",
                                                    fw_name_constraint=wf.fws[4].name)


                                print(wf)
                                lpad.add_wf(wf)

if __name__ == '__main__':
    pc()