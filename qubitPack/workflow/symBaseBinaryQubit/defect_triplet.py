from qubitPack.workflow.antisiteQubit.wf_owls import DefectWF
from qubitPack.workflow.tool_box import find_cation_anion

from fireworks import LaunchPad

from atomate.vasp.database import VaspCalcDb
from atomate.vasp.workflows.jcustom.hse_full import get_wf_full_hse
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    preserve_fworker,
    add_modify_incar,
    add_modify_kpoints,
    set_queue_options,
    set_execution_options
)

from pymatgen import Structure

from pycdt.core.defectsmaker import ChargedDefectsStructures






class SymBaseBinaryQubit(DefectWF):
    # def __init__(self, orig_st, natom, defect_type, substitution, distort=0, vacuum_thickness=None):
    #     super().__init__(orig_st, natom, defect_type, substitution, distort, vacuum_thickness)

    @classmethod
    def binary_vacancy(cls, cat="c2dbVac"):
        lpad = LaunchPad.from_file("/home/tug03990/config/project/symBaseBinaryQubit/c2dbVac/my_launchpad.yaml")
        col = VaspCalcDb.from_db_file("/home/tug03990/config/category/mx2_antisite_pc/db.json").collection
        # mx2s = col.find({"task_id":{"$in":[3091, 3083, 3093, 3097, 3094, 3102]}})
        # 3091: S-W, 3083: Se-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te
        mx2s = col.find({"task_id":{"$in":[3091, 3083, 3093, 3097, 3094, 3102]}})

        geo_spec = {5* 5 * 3: [15]}
        aexx = 0.25
        for mx2 in mx2s:
            pc = Structure.from_dict(mx2["output"]["structure"])

            defect = ChargedDefectsStructures(pc, antisites_flag=True).defects

            cation, anion = find_cation_anion(pc)

            for vac in range(len(defect["substitutions"])):
                print(cation, anion)
                # cation vacancy
                if "{}".format(cation) not in defect["vacancies"][vac]["name"]:
                    continue
                for na, thicks in geo_spec.items():
                    for thick in thicks:
                        for dtort in [0, 0.001]:
                            vacancy_insttance = cls(pc, na, ("vacancies", vac), False, dtort, thick)
                            wf = get_wf_full_hse(
                                structure=vacancy_insttance.defect_st,
                                charge_states=[0],
                                gamma_only=False,
                                dos_hse=True,
                                nupdowns=[2],
                                encut=320,
                                include_hse_relax=True,
                                vasptodb={"category": cat, "NN": vacancy_insttance.NN,
                                          "defect_entry": vacancy_insttance.defect_entry},
                                wf_addition_name="{}:{}".format(vacancy_insttance.defect_st.num_sites, thick)
                            )

                            def kpoints(kpts):
                                kpoints_kwarg = {
                                    'comment': "mx2_antisite",
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

                            wf = add_modify_kpoints(wf, {"kpoints_update": kpoints([1,1,1])}, "PBE_relax")
                            wf = add_modify_kpoints(wf, {"kpoints_update": kpoints([1,1,1])}, "HSE_relax")
                            wf = add_modify_kpoints(wf, {"kpoints_update": kpoints([1,1,1])}, "HSE_scf")
                            wf = add_additional_fields_to_taskdocs(
                                wf,
                                {"lattice_constant": "HSE",
                                 "perturbed": vacancy_insttance.distort}
                            )

                            wf = add_modify_incar(wf, {"incar_update": {"NSW":150}}, "PBE_relax")
                            wf = add_modify_incar(wf, {"incar_update": {"AEXX":aexx, "NSW":150}}, "HSE_relax")
                            wf = add_modify_incar(wf, {"incar_update": {"LWAVE": True, "AEXX":aexx}}, "HSE_scf")
                            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="PBE_relax")
                            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_relax")
                            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")
                            # related to directory
                            wf = set_execution_options(wf, category=cat)
                            wf = preserve_fworker(wf)
                            wf.name = wf.name+":dx[{}]".format(vacancy_insttance.distort)
                            lpad.add_wf(wf)


if __name__ == '__main__':
    SymBaseBinaryQubit.binary_vacancy()