from qubitPack.workflow.tool_box import *

from fireworks import LaunchPad, Workflow

from atomate.vasp.database import VaspCalcDb
from atomate.vasp.workflows.jcustom.wf_full import get_wf_full_hse, get_wf_full_scan
from atomate.vasp.fireworks.core import ScanOptimizeFW
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    add_modify_kpoints,
    set_queue_options,
    set_execution_options,
    clean_up_files
)

from pymatgen import Structure
from pymatgen.io.vasp.sets import MPScanRelaxSet, MPHSERelaxSet

from pycdt.core.defectsmaker import ChargedDefectsStructures

from monty.serialization import loadfn
import os


def relax_pc():
    lpad = LaunchPad.from_file("/home/tug03990/atomate/example/config/project/"
                               "symBaseBinaryQubit/scan_relax_pc/my_launchpad.yaml")

    mx2s = loadfn("/home/tug03990/atomate/example/config/project/symBaseBinaryQubit/"
                  "scan_relax_pc/gap_gt1-binary-NM.json")

    for mx2 in mx2s:
        if mx2["irreps"] and mx2["formula"] == "Rh2Br6" and mx2["spacegroup"] == "P3":
            pc = mx2["structure"]
            pc = modify_vacuum(pc, 20)
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


def binary_scan_defect(cat="binary_defect", defect_type=("vacancies", "S"), impurity_on_nn=None): #BN_vac
    lpad = LaunchPad.from_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/defect_db/{}/my_launchpad.yaml".format(cat)))
    col = VaspCalcDb.from_db_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/symBaseBinaryQubit/scan_relax_pc/db.json")).collection

    mx2s = col.find(
        {
            "task_id":60
        }
    )
    geo_spec = None
    aexx = 0.25
    test = []
    for mx2 in mx2s:
        pc = Structure.from_dict(mx2["output"]["structure"])
        if mx2["nsites"] == 2:
            geo_spec = {25*2: [20]}
        if mx2["nsites"] == 3:
            geo_spec = {25*3: [20]}
        if mx2["nsites"] == 4:
            geo_spec = {25*4: [20]}
        if mx2["nsites"] == 8:
            geo_spec = {25*8: [20]}

        defect = ChargedDefectsStructures(pc, antisites_flag=True).defects

        cation, anion = find_cation_anion(pc)
        good_sym_site_symbols = list(dict.fromkeys(mx2["c2db_info"]["irreps"]))
        print(good_sym_site_symbols)
        for good_sym_site_symbol in good_sym_site_symbols:
            defect_type = ("substitutions", "on_{}".format(good_sym_site_symbol))
            for de_idx in range(len(defect[defect_type[0]])):
                print(cation, anion)
                # cation vacancy
                if defect_type[1] == "_".join(defect[defect_type[0]][de_idx]["name"].split("_")[-2:]):
                    for na, thicks in geo_spec.items():
                        for thick in thicks:
                            for dtort in [0]:
                                if not impurity_on_nn:
                                    impurity_on_nn = []
                                gen_defect = GenDefect(pc, [defect_type[0], de_idx], na, thick, distort=dtort,
                                                       sub_on_side=list(impurity_on_nn))
                                # gen_defect.vacancies(dtort, list(impurity_on_nn))
                                enmax = 1.3*max([potcar.enmax for potcar in MPScanRelaxSet(gen_defect.defect_st).potcar])

                                defect_data = gen_defect.defect_entry

                                wf = get_wf_full_scan(
                                    structure=gen_defect.defect_st,
                                    charge_states=[0],
                                    gamma_only=False,
                                    dos=True,
                                    nupdowns=[-1],
                                    encut=enmax,
                                    vasptodb={
                                        "category": cat,
                                        "NN": gen_defect.NN,
                                        "NN_dist": gen_defect.nn_dist,
                                        "defect_entry": defect_data,
                                    },
                                    wf_addition_name="{}:{}".format(gen_defect.defect_st.num_sites, thick)
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
                                     "perturbed": gen_defect.distort}
                                )

                                wf = add_modify_incar(wf, {"incar_update": {"NSW":150}}, "SCAN_relax")
                                wf = add_modify_incar(wf, {"incar_update": {"LWAVE": False}}, "SCAN_scf")
                                wf = set_queue_options(wf, "24:00:00", fw_name_constraint="SCAN_relax")
                                wf = set_queue_options(wf, "24:00:00", fw_name_constraint="SCAN_scf")
                                # related to directory
                                wf = set_execution_options(wf, category=cat)
                                wf = preserve_fworker(wf)
                                wf.name = wf.name+":dx[{}]".format(gen_defect.distort)
                                # task_name_constraint=x meaning after x do powersup
                                wf = clean_up_files(wf, files=["CHG*", "DOS*", "LOCPOT*"], fw_name_constraint="SCAN_scf",
                                                    task_name_constraint="VaspToDb")

                                print(wf)
                                # lpad.add_wf(wf)


if __name__ == '__main__':
    # relax_pc()
    binary_scan_defect()