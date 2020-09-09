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
    set_execution_options
)

from pymatgen import Structure
from pymatgen.io.vasp.sets import MPScanRelaxSet, MPHSERelaxSet

from pycdt.core.defectsmaker import ChargedDefectsStructures

from monty.serialization import loadfn


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


class SymBaseBinaryQubit:
    @staticmethod
    def binary_vacancy(cat="BN_vac"):
        lpad = LaunchPad.from_file("/home/tug03990/config/project/symBaseBinaryQubit/BN_vac/my_launchpad.yaml")
        col = VaspCalcDb.from_db_file('/home/tug03990/atomate/example/config/project/symBaseBinaryQubit/scan_relax_pc/'
                                      'db.json').collection

        mx2s = col.find(
            {
                "task_id": {"$in": [1]}
            }
        )

        geo_spec = None
        aexx = 0.25
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

            for vac in range(len(defect["vacancies"])):
                print(cation, anion)
                # cation vacancy
                if "{}".format("N") not in defect["vacancies"][vac]["name"]:
                    continue
                print(vac)
                for na, thicks in geo_spec.items():
                    for thick in thicks:
                        for dtort in [0]:
                            vacancies = GenDefect(pc, ["vacancies", vac], na, thick)
                            vacancies.vacancies(dtort, "C")
                            enmax = 1.3*max([potcar.enmax for potcar in MPHSERelaxSet(vacancies.defect_st).potcar])
                            wf = get_wf_full_hse(
                                structure=vacancies.defect_st,
                                charge_states=[0],
                                gamma_only=False,
                                dos_hse=True,
                                nupdowns=[2],
                                encut=enmax,
                                include_hse_relax=True,
                                vasptodb={
                                    "category": cat,
                                    "NN": vacancies.NN,
                                    "NN_dist": vacancies.nn_dist,
                                    "defect_entry": vacancies.defect_entry
                                },
                                wf_addition_name="{}:{}".format(vacancies.defect_st.num_sites, thick)
                            )

                            def kpoints(kpts):
                                kpoints_kwarg = {
                                    'comment': "JCustom",
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
                                {"lattice_constant": "SCAN",
                                 "perturbed": vacancies.distort}
                            )

                            wf = add_modify_incar(wf, {"incar_update": {"NSW":150}}, "PBE_relax")
                            wf = add_modify_incar(wf, {"incar_update": {"AEXX":aexx, "NSW":150}}, "HSE_relax")
                            wf = add_modify_incar(wf, {"incar_update": {"LWAVE": False, "AEXX":aexx}}, "HSE_scf")
                            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="PBE_relax")
                            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_relax")
                            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")
                            # related to directory
                            wf = set_execution_options(wf, category=cat)
                            wf = preserve_fworker(wf)
                            wf.name = wf.name+":dx[{}]".format(vacancies.distort)
                            lpad.add_wf(wf)

    @staticmethod
    def binary_vacancy_scan(cat="BN_vac"):
        lpad = LaunchPad.from_file("/home/tug03990/config/project/symBaseBinaryQubit/BN_vac/my_launchpad.yaml")
        col = VaspCalcDb.from_db_file('/home/tug03990/atomate/example/config/project/symBaseBinaryQubit/scan_relax_pc/'
                                      'db.json').collection

        mx2s = col.find(
            {
                "task_id": {"$in": [1]}
            }
        )

        geo_spec = None
        aexx = 0.25
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

            for vac in range(len(defect["vacancies"])):
                print(cation, anion)
                # cation vacancy
                if "{}".format("B") not in defect["vacancies"][vac]["name"]:
                    continue
                print(vac)
                for na, thicks in geo_spec.items():
                    for thick in thicks:
                        for dtort in [0]:
                            vacancies = GenDefect(pc, ["vacancies", vac], na, thick)
                            vacancies.vacancies(dtort, False)
                            enmax = 1.3*max([potcar.enmax for potcar in MPScanRelaxSet(vacancies.defect_st).potcar])
                            wf = get_wf_full_scan(
                                structure=vacancies.defect_st,
                                charge_states=[0],
                                gamma_only=True,
                                dos_hse=True,
                                nupdowns=[-1],
                                encut=enmax,
                                vasptodb={
                                    "category": cat,
                                    "NN": vacancies.NN,
                                    "NN_dist": vacancies.nn_dist,
                                    "defect_entry": vacancies.defect_entry
                                },
                                wf_addition_name="{}:{}".format(vacancies.defect_st.num_sites, thick)
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
                                {"lattice_constant": "SCAN",
                                 "perturbed": vacancies.distort}
                            )

                            wf = add_modify_incar(wf, {"incar_update": {"NSW":150}}, "SCAN_relax")
                            wf = add_modify_incar(wf, {"incar_update": {"LWAVE": False}}, "SCAN_scf")
                            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="SCAN_relax")
                            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="SCAN_scf")
                            # related to directory
                            wf = set_execution_options(wf, category=cat)
                            wf = preserve_fworker(wf)
                            wf.name = wf.name+":dx[{}]".format(vacancies.distort)
                            lpad.add_wf(wf)

if __name__ == '__main__':
    # relax_pc()
    SymBaseBinaryQubit.binary_vacancy_scan()