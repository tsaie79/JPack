from fireworks import Workflow
from fireworks import LaunchPad

from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import *
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW, ScanOptimizeFW
from atomate.vasp.fireworks.jcustom import *
from atomate.vasp.workflows.jcustom.wf_full import get_wf_full_hse
from atomate.vasp.config import RELAX_MAX_FORCE
from atomate.vasp.database import VaspCalcDb

from pymatgen.io.vasp.inputs import Structure, Kpoints
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet, MPHSERelaxSet, MPHSEBSSet

from VaspBandUnfolding import unfold

import os, copy
from glob import glob
from monty.serialization import loadfn
from unfold import find_K_from_k
import math
from qubitPack.tool_box import *


class HostWF:
    @classmethod
    def hse_scf_wf(cls):
        col = VaspCalcDb.from_db_file("/home/tug03990/config/category/mx2_antisite_pc/db.json").collection
        # 3091: S-W, 3083: Se-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te
        mx2s = col.find({"task_id":{"$in":[3091]}})

        for mx2 in mx2s:
            pc = Structure.from_dict(mx2["output"]["structure"])
            pc = modify_vacuum(pc, 20)
            scf_wf = get_wf_full_hse(
                structure=pc,
                charge_states=[0],
                gamma_only=False,
                gamma_mesh=True,
                scf_dos=True,
                nupdowns=[-1],
                task="hse_relax",
                category="pc",
            )
            print(scf_wf)

            wf = add_modify_incar(scf_wf,
                                  {
                                      "incar_update":
                                          {"LASPH": True, "EDIFF": 1E-7, "EDIFFG": -0.001, "ISIF": 3, "NSW": 250, "LCHARG":False, "LWAVE":False}})
            LPAD.add_wf(wf)

    @classmethod
    def hse_band_structure(cls):
        pc = Structure.from_dict({'@module': 'pymatgen.core.structure',
                                  '@class': 'Structure',
                                  'charge': None,
                                  'lattice': {'matrix': [[3.1529358858637777, 1.0053975e-09, 0.0],
                                                         [-1.5764674411892312, 2.730521508558374, 0.0],
                                                         [0.0, 0.0, 20.0]],
                                              'a': 3.1529358858637777,
                                              'b': 3.1529347125859775,
                                              'c': 20.0,
                                              'alpha': 90.0,
                                              'beta': 90.0,
                                              'gamma': 120.00000176314633,
                                              'volume': 172.18318506083142},
                                  'sites': [{'species': [{'element': 'W', 'occu': 1}],
                                             'abc': [0.0, 0.0, 0.5],
                                             'xyz': [0.0, 0.0, 10.0],
                                             'label': 'W',
                                             'properties': {}},
                                            {'species': [{'element': 'S', 'occu': 1}],
                                             'abc': [0.6666669999999968, 0.3333330000000032, 0.5778428253948604],
                                             'xyz': [1.5764696866472017, 0.9101729266825626, 11.556856507897209],
                                             'label': 'S',
                                             'properties': {}},
                                            {'species': [{'element': 'S', 'occu': 1}],
                                             'abc': [0.6666669999999968, 0.3333330000000032, 0.4221571746051396],
                                             'xyz': [1.5764696866472017, 0.9101729266825626, 8.443143492102791],
                                             'label': 'S',
                                             'properties': {}}]}
                                 )
        bs_wf = get_wf_full_hse(
            structure=pc,
            charge_states=[0],
            gamma_only=False,
            gamma_mesh=True,
            scf_dos=False,
            nupdowns=[-1],
            task="hse_scf-hse_bs",
            category="pc",
        )
        bs_wf = add_modify_incar(bs_wf,{"incar_update": { "EDIFF": 1E-7, "LCHARG":True, "LWAVE":False}}, "HSE_scf")
        bs_wf = add_modify_incar(bs_wf,{"incar_update": { "EDIFF": 1E-7, "LCHARG":False, "LWAVE":False}}, "HSE_bs")
        LPAD.add_wf(bs_wf)

class PBEDefectWF:
    def get_wf_point_defects(self, bulk_structure, defect_structure,
                             charge_states=None, name="point_defects",
                             lepsilon=False, vasp_cmd="vasp", db_file=None, user_kpoints_settings=None,
                             tag=""):
        # defect
        fws = []
        vis_relax = MPRelaxSet(defect_structure)
        nelect = vis_relax.nelect
        for cs in charge_states:
            uis_relax = {"EDIFFG": -0.02, "EDIFF": 1E-4, "ISIF": 2, "NELECT": nelect - cs, "ENCUT": 400, "ISPIN":1}
            v = vis_relax.as_dict()
            v.update({"user_incar_settings": uis_relax})
            vis_relax = vis_relax.__class__.from_dict(v)
            defect_optimize_fw = OptimizeFW(structure=defect_structure, job_type="normal",
                                            vasp_input_set=vis_relax, name="{} {} structure optimization".format(tag, cs))
            uis_static = {"EDIFF": 1E-5, "NELECT": nelect - cs, "ENCUT": 400, "ISPIN":1}
            vis_static = MPStaticSet(defect_structure, force_gamma=False, lepsilon=False,
                                     user_kpoints_settings=user_kpoints_settings,
                                     user_incar_settings=uis_static)

            defect_static_fw = StaticFW(defect_structure, vasp_input_set=vis_static, parents=defect_optimize_fw,
                                        name="{} {} static calculation".format(tag, cs))
            fws.append(defect_optimize_fw)
            fws.append(defect_static_fw)

        # bulk
        # vis_relax = MPRelaxSet(bulk_structure)
        # uis_relax = {"EDIFFG": -0.02, "EDIFF": 1E-4, "ISIF": 2, "ISPIN":1}
        # v = vis_relax.as_dict()
        # v.update({"user_incar_settings": uis_relax})
        # vis_relax = vis_relax.__class__.from_dict(v)
        # bulk_opt_fw = OptimizeFW(structure=bulk_structure, job_type="normal",
        #                          vasp_input_set=vis_relax, name="{} bulk structure optimization".format(tag))
        #
        # uis_static = {"EDIFF": 1E-5, "ISPIN":1}
        # vis_static = MPStaticSet(bulk_structure, force_gamma=False, lepsilon=lepsilon,
        #                          user_kpoints_settings=user_kpoints_settings,
        #                          user_incar_settings=uis_static)
        #
        # bulk_static_fw = StaticFW(bulk_structure, vasp_input_set=vis_static, parents=bulk_opt_fw,
        #                          name="{} bulk static calculation".format(tag))
        # fws.append(bulk_opt_fw)
        # fws.append(bulk_static_fw)

        wfname = "{}:{}".format(defect_structure.composition.reduced_formula, name)

        return Workflow(fws, name=wfname)

    def bn_sub_wf(self, lzs):
        lpad = LaunchPad.auto_load()
        st = Structure.from_file("/gpfs/work/tug03990/2D_formation_energy_corr/host/POSCAR")
        chg = ChargedDefectsStructures(st, substitutions={"N": "C"}, cellmax=9*9*2)

        for lz in lzs:
            bulk = chg.get_ith_supercell_of_defect_type(0, "bulk")
            bulk = modify_vacuum(bulk, lz)
            print(bulk.lattice.c)
            # defect = chg.get_ith_supercell_of_defect_type(0, "substitutions")
            defect = chg.get_ith_supercell_of_defect_type(0, "vacancies")
            defect = modify_vacuum(defect, lz)
            print(defect.lattice.c)
            wf = self.get_wf_point_defects(bulk, defect, [0, -1], name="{} point_defects".format(lz))
            lpad.add_wf(wf)


class DefectWF:

    @classmethod
    def relax_pc(cls):
            # lpad = LaunchPad.from_file("/home/tug03990/config/my_launchpad.efrc.yaml")
            lpad = LaunchPad.from_file("/home/tug03990/config/project/antisiteQubit/scan_opt_test/my_launchpad.yaml")
            col = VaspCalcDb.from_db_file("/home/tug03990/PycharmProjects/my_pycharm_projects/database/db_config/"
                                          "db_dk_local.json").collection
            # defect_st_col = VaspCalcDb.from_db_file("/home/tug03990/PycharmProjects/my_pycharm_projects/database/db_config/"
            #                               "db_c2db_tmdc_bglg1.json").collection
            mx2s = col.find(
                {"class": {"$in": ["TMDC-T", "TMDC-H", "TMDC-T'"]},
                 "gap_hse": {"$gt": 1},
                 "ehull": {"$lt": 0.3},
                 "magstate": "NM",
                 # "formula": {"$in":["WTe2"]}
                 # "spacegroup": "P-6m2",
                 "formula": {"$in": ["MoS2", "MoSe2",  "MoTe2", "WS2", "WSe2", "WTe2"]}
                 # "formula": {"$in": ["WSe2"]}
                 }
            )


            for mx2 in mx2s:
                pc = Structure.from_dict(mx2["structure"])
                # pc = Structure.from_file("/home/tug03990/work/sandwich_BN_mx2/structures/BN_c2db.vasp")
                pc = modify_vacuum(pc, 20)

                def mphserelaxset(aexx):
                    vis_relax = MPHSERelaxSet(pc, force_gamma=True)
                    v = vis_relax.as_dict()
                    v.update({"user_incar_settings":{"AEXX": aexx, "ALGO":"All"}})
                    vis_relax = vis_relax.__class__.from_dict(v)
                    return vis_relax

                scan_opt = ScanOptimizeFW(structure=pc, name="SCAN_relax")

                # pbe_relax = OptimizeFW(structure=pc, name="PBE_relax")
                # hse_relax_25 = OptimizeFW(structure=pc, name="HSE_relax", vasp_input_set=mphserelaxset(0.3),
                #                           parents=pbe_relax)
                # hse_relax_35 = OptimizeFW(structure=pc, name="HSE_relax", vasp_input_set=mphserelaxset(0.35), parents=pbe_relax)
                wf = Workflow([scan_opt], name="{}:SCAN_opt".format(mx2["formula"]))
                # wf = add_modify_incar(wf, {"incar_update":{"NCORE":4, "NSW":100, "ISIF":3, "EDIFF":1E-7,
                #                                            "EDIFFG":-0.001,
                #                                            "LCHARG":False, "LWAVE":False}})
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

                wf = set_execution_options(wf, category="scan_opt_test")
                lpad.add_wf(wf)

    @classmethod
    def MX2_anion_antisite(cls, defect_type="vacancies", category="vac_tmd", distort=0, add_irvsp_fw=False):
        # for analysis and plotting, use script in qubitPak/defect_formation_energy_correction
        lpad = LaunchPad.from_file("/home/tug03990/config/project/antisiteQubit/{}/my_launchpad.yaml".format(category))
        aexx = 0.25

        col = VaspCalcDb.from_db_file("/home/tug03990/config/category/mx2_antisite_pc/db.json").collection
        # 3091: S-W, 3083: Se-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te
        mx2s = col.find({"task_id":{"$in":[3093, 3102]}})
        geo_spec = {5*5*3: [20]}
        for mx2 in mx2s:
            # pc = Structure.from_dict(mx2["structure"])
            pc = Structure.from_dict(mx2["output"]["structure"])
            defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
            cation, anion = find_cation_anion(pc)
            for idx, antisite in enumerate(range(len(defect[defect_type]))):
                print(cation, anion)
                if "vac_1_{}".format(anion) not in defect[defect_type][antisite]["name"]:
                    continue
                for na, thicks in geo_spec.items():
                    sc_size = np.eye(3, dtype=int)*math.sqrt(na/3)

                    for thick in thicks:
                        antisite_st = GenDefect(
                            orig_st=pc,
                            defect_type=(defect_type, antisite),
                            natom=na,
                            vacuum_thickness=thick,
                            sub_on_side=None,
                            distort=distort
                        )

                        wf_antisite = get_wf_full_hse(
                            structure=antisite_st.defect_st,
                            charge_states=[-2, -2],
                            gamma_only=None, #[(0,0,0)],#[list(find_K_from_k(sc_size, [0, 0, 0])[0])],
                            gamma_mesh=True,
                            nupdowns=[0, 2],
                            task="hse_relax-hse_scf",
                            vasptodb={"category": category, "NN": antisite_st.NN,
                                      "defect_entry": antisite_st.defect_entry, "defect_type": "defect"},
                            wf_addition_name="{}:{}".format(na, thick),
                            task_arg={"parse_dos":True}
                        )
                        if add_irvsp_fw:
                            from pytopomat.workflows.fireworks import IrvspFW
                            from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
                            fws = wf_antisite.fws
                            for i in range(2):
                                fw_irvsp = IrvspFW(
                                    structure=antisite_st.defect_st,
                                    parents=fws[i+1],
                                    irvsptodb_kwargs={
                                        "additional_fields":
                                            {
                                                "spg_pymatgen": SpacegroupAnalyzer(antisite_st.defect_st).get_space_group_symbol()
                                            }
                                },

                                )
                                fws.append(fw_irvsp)

                            wf_antisite = Workflow(fws, name=wf_antisite.name)
                            wf_antisite = add_modify_kpoints(wf_antisite, {
                                "kpoints_update" :{
                                "num_kpts": 2,
                                "kpts": [(0,0,0), (0,0,0)],
                                "kpts_weights": [1, 0],
                                "labels": [None, "\Gamma"]}}, fw_name_constraint="HSE_scf")
                            wf_antisite = add_modify_incar(
                                wf_antisite,
                                {"incar_update":{"LWAVE": True}},
                                fw_name_constraint="HSE_scf"
                            )

                        wf_antisite = add_modify_incar(
                            wf_antisite,
                            {"incar_update":{"AEXX": aexx}},
                        )
                        wf_antisite = add_modify_incar(wf_antisite)
                        wf_antisite = set_queue_options(wf_antisite, "24:00:00", fw_name_constraint="HSE_scf")
                        wf_antisite = set_queue_options(wf_antisite, "24:00:00", fw_name_constraint="HSE_relax")
                        # related to directory
                        wf_antisite = set_execution_options(wf_antisite, category=category, fworker_name="owls")
                        wf_antisite = preserve_fworker(wf_antisite)

                        wf = add_additional_fields_to_taskdocs(
                            wf_antisite,
                            {
                                "defect_obj": antisite_st.pmg_obj.as_dict(),
                                "lattice_constant": "HSE",
                                "perturbed": antisite_st.distort,
                                "pc_from": "mx2_antisite_pc/HSE_scf/{}".format(mx2["task_id"])
                            }
                        )
                        # wf = scp_files(
                        #     wf,
                        #     "/home/jengyuantsai/Research/projects/single_photon_emitter/standard_defect",
                        #     fw_name_constraint="HSE_scf",
                        # )
                        # lpad.add_wf(wf)

    @classmethod
    def defct_static(cls):
        from qubitPack.tool_box import get_db
        db = get_db("antisiteQubit", "W_Se_Ef_gamma", port=12345)
        st = db.collection.find_one({"task_id":153})["output"]["structure"]
        category = "standard_defect"
        wf = get_wf_full_hse(
            structure=Structure.from_dict(st),
            charge_states=[0],
            gamma_only=True,
            gamma_mesh=False,
            task_arg={"parse_dos":False},
            nupdowns=[-1],
            task="hse_relax-hse_scf",
            vasptodb={"category": category},
            wf_addition_name=None,
        )

        wf = add_additional_fields_to_taskdocs(wf, {"pc_from": "antisiteQubit/W_Se_Ef_gamma/153"})
        wf = scp_files(
            wf,
            "/home/jengyuantsai/Research/projects/single_photon_emitter/standard_defect",
            fw_name_constraint="HSE_scf",
        )
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")
        wf = add_modify_incar(wf, {"incar_update": {"LWAVE":True, "LCHARG":False}}, "HSE_scf")
        # related to directory
        wf = preserve_fworker(wf)
        wf = set_execution_options(wf, fworker_name="owls", category=category)
        lpad = LaunchPad.from_file("/home/tug03990/config/project/single_photon_emitter/standard_defect/my_launchpad.yaml")
        lpad.add_wf(wf)

    @classmethod
    def MX2_anion_vacancy(cls,cat):

            # lpad = LaunchPad.auto_load()
            lpad = LaunchPad.from_file("/home/tug03990/config/project/antisiteQubit/MxC3vToChDeltaE/my_launchpad.yaml")
            col = VaspCalcDb.from_db_file("/home/tug03990/PycharmProjects/my_pycharm_projects/database/db_config/"
                                          "db_mx2_antisite_pc.json").collection
            # mx2s = col.find({"task_id":{"$in":[3091, 3083, 3093, 3097, 3094, 3102]}})
            # 3091: S-W, 3083: Se-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te
            mx2s = col.find({"task_id":{"$in":[3091, 3083, 3093, 3097, 3094, 3102]}})

            # col = VaspCalcDb.from_db_file("/home/tug03990/config/category/mx2_antisite_basic_aexx0.25_final/db.json").collection
            # # mx2s = col.find({"task_id":{"$in":[3281, 3282, 3291, 3285]}}) #3281, 3282, 3281, 3291, 3285
            # mx2s = col.find({"task_id":{"$in":[3302]}})
            # col = VaspCalcDb.from_db_file("/home/tug03990/config/category/mx2_antisite_basic/db.json").collection
            # mx2s = col.find({"task_id":{"$in":[2365]}})
            # geo_spec = {5 * 5 * 3: [20, 30, 40], 6 * 6 * 3: [20, 30, 40]}
            geo_spec = {5 * 5 * 3: [20]}
            aexx = 0.25
            for mx2 in mx2s:
                # pc = Structure.from_dict(mx2["structure"])
                pc = Structure.from_dict(mx2["output"]["structure"])
                # if "Te" in pc.formula:
                #     pc = special_treatment_to_structure(pc, "selective_dynamics", NN=[74, 55, 49, 54])
                # else:
                #     pc = special_treatment_to_structure(pc, "selective_dynamics", NN=[0, 6, 5, 25])

                defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
                cation, anion = find_cation_anion(pc)

                for vacancy in range(len(defect["vacancies"])):
                    print(cation, anion)
                    # cation vacancy
                    if cation not in defect["vacancies"][vacancy]["name"]:
                        continue
                    for na, thicks in geo_spec.items():
                        for thick in thicks:
                            for dtort in [0, 0.001]:
                                se_antisite = DefectWF(pc, natom=na, vacuum_thickness=thick, substitution=None,
                                                       antisite=False, type_vac=vacancy, bulk=False, distort=dtort)

                                wf = get_wf_full_hse(
                                    structure=se_antisite.defect_st,
                                    charge_states=[0],
                                    gamma_only=True,
                                    dos_hse=True,
                                    nupdowns=[2],
                                    encut=320,
                                    include_hse_relax=True,
                                    vasptodb={"category": cat, "NN": se_antisite.NN},
                                    wf_addition_name="{}:{}".format(na, thick)
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
                                     "perturbed": se_antisite.distort}
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
                                wf.name = wf.name+":dx[{}]".format(se_antisite.distort)
                                lpad.add_wf(wf)
                    # break

    @classmethod
    def MX2_formation_energy(cls, category="formation_energy"):
            # for analysis and plotting, use script in qubitPak/defect_formation_energy_correction
            lpad = LaunchPad.from_file("/home/tug03990/config/project/antisiteQubit/{}/my_launchpad.yaml".format(category))
            aexx = 0.25

            col = VaspCalcDb.from_db_file("/home/tug03990/config/category/mx2_antisite_pc/db.json").collection
            # 4229: S-W, 4239: Se-W, 4236: Te-W, 4237:Mo-S, 4238: Mo-Se, 4235:Mo-Te
            mx2s = col.find({"task_id":{"$in":[4239]}})
            geo_spec = {5*5*3: [20]}
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
                                    gamma_only=[list(find_K_from_k(sc_size, [1/3, 1/3, 0])[0])],
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

    @classmethod
    def other_defects(cls):
        def WS2(charge=-2):
            from qubitPack.tool_box import get_db
            from mpinterfaces.utils import ensure_vacuum
            db = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", user="adminUser", password="qiminyan", port=12345)
            e = db.collection.find_one({"uid": "WS2-MoS2-NM"})
            host = Structure.from_dict(e["structure"])
            host = ensure_vacuum(host, 15)
            from pymatgen.analysis.defects.core import Vacancy, Substitution, Interstitial
            from pymatgen.io.vasp import Poscar
            from pymatgen import MPRester, PeriodicSite
            defect_site = PeriodicSite(host.sites[1].specie, host.sites[1].frac_coords,
                                       host.lattice)
            vac = Vacancy(host, defect_site, charge=charge)
            print("Creating Defect {} in charge state {}:".format(vac.name, vac.charge))
            vac = vac.generate_defect_structure(supercell=(4, 4, 1))
            from qubitPack.tool_box import phonopy_structure
            vac = phonopy_structure(vac)
            vac.set_charge(charge)
            return vac

        category = "vac_tmd"
        lpad = LaunchPad.from_file("/home/tug03990/config/project/antisiteQubit/{}/my_launchpad.yaml".format(category))
        from atomate.vasp.workflows.base.core import get_wf
        from atomate.vasp.config import GAMMA_VASP_CMD
        from atomate.vasp.jpowerups import clear_to_db
        wf = get_wf(WS2(), "/home/tug03990/scripts/JPack/projects/antisiteQubit/wf_yaml/static.yaml")
        wf = add_modify_incar(wf)
        # wf = clear_to_db(wf)
        wf = set_queue_options(wf, "24:00:00")
        wf = set_execution_options(wf, category=category)
        wf = preserve_fworker(wf)
        wf = use_gamma_vasp(wf, GAMMA_VASP_CMD)
        wf = add_modify_kpoints(wf, modify_kpoints_params={"kpoints_update": {"kpts": [[1,1,1]], "style": "G"}})
        lpad.add_wf(wf)

    @classmethod
    def antisite_wse2_singlet_triplet(cls):
            primitive_st_wse2 = Structure.from_file("/gpfs/work/tug03990/research/"
                                                    "BinaryMaterials/OptimizeLatticeConst/noO_hexagonal/"
                                                    "electron_even_odd/even_electron/040752_Se2W1/relax_uc/CONTCAR")
            wf = []
            wf_name = "{}:{}".format(primitive_st_wse2.composition.reduced_formula, "0th_perturbed")
            for nupdown_set in [0]:
                antisite = cls(primitive_st_wse2, 4*4*3, 0, False, True, 21.311627)
                fws = antisite.hse_scf_wflow([0], [0, 0, 0], False, "0th-perturbed", nupdown_set).fws
                for fw in fws:
                    wf.append(fw)
            wf = Workflow(wf, name=wf_name)
            wf = add_namefile(wf)
            lpad = LaunchPad.auto_load()
            lpad.add_wf(wf)

    @classmethod
    def wse2_bulk(cls):
            primitive_st_wse2 = Structure.from_file("/gpfs/work/tug03990/research/"
                                                    "BinaryMaterials/OptimizeLatticeConst/noO_hexagonal/"
                                                    "electron_even_odd/even_electron/040752_Se2W1/relax_uc/CONTCAR")
            wf = []
            wf_name = "{}:{}".format(primitive_st_wse2.composition.reduced_formula, "4*4-21.312")
            for nupdown_set in [-1]:
                antisite = cls(primitive_st_wse2, 4 * 4 * 3, 0, False, False, 21.311627, bulk=True)
                fws = antisite.hse_scf_wflow([0], [0.25, 0.25, 0], True, "bulk", nupdown_set).fws
                for fw in fws:
                    wf.append(fw)
            wf = Workflow(wf, name=wf_name)
            wf = add_namefile(wf)
            lpad = LaunchPad.auto_load()
            lpad.add_wf(wf)

    @classmethod
    def C2DB_vacancies_db_wflow(cls):
            primitive_list = loadfn("/home/tug03990/work/vacancies_db/c2db_vacancy_V1.json")
            primitive_list = primitive_list

            for uid, info in primitive_list.items():
                nsite = info[1]
                st = info[0]
                formula = info[2]
                gap_hse = round(info[3], 3)
                defect = ChargedDefectsStructures(st).defects
                fws = []
                for vac_type in range(len(defect["vacancies"])):
                    print("**"*20, vac_type)
                    vac = cls(st, 4*4*nsite, vac_type, None, False, 20)
                    print(vac.defect_st.lattice.c)
                    if abs(vac.defect_st.lattice.c - 20) > 0.001:
                        return
                    wf_name = "{}:{}".format(st.composition.reduced_formula, gap_hse)
                    fw = vac.hse_scf_no_relax_wf([0], True, True, {"formula": formula, "gap_hse":gap_hse})
                    for fw_sole in fw.fws:
                        fws.append(fw_sole)

                wf = Workflow(fws, name=wf_name)
                wf = add_namefile(wf)

                lpad = LaunchPad.auto_load()
                lpad.add_wf(wf)

    @classmethod
    def qimin_db_even_antisite_lg_bg(cls):
            col = VaspCalcDb.from_db_file("/home/tug03990/PycharmProjects/my_pycharm_projects/database/db_config/"
                                          "db_qimin.json").collection
            even = loadfn('/gpfs/work/tug03990/qimin_db_even_antisite/even_large_bg.js')
            odd = loadfn('/gpfs/work/tug03990/qimin_db_even_antisite/odd_large_bg.js')
            all_lg_bg = [i.split("_")[0] for i in list(even.values()) + list(odd.values())]
            path = '/gpfs/work/tug03990/qimin_db_even_antisite/structures'
            for icsd_id in all_lg_bg:
                distinct_pointgps = col.distinct("output.spacegroup.point_group",
                                                 {"icsd_id":icsd_id, "task_type":"static"})
                print(distinct_pointgps)
                for point_gp in distinct_pointgps:
                    e = col.find_one({"icsd_id":icsd_id, "task_type":"static", "output.spacegroup.point_group": point_gp})
                    nsite = e["nsites"]
                    pc = Structure.from_dict(e["input"]["structure"])
                    defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
                    fws = []

                    for antisite in range(len(defect["substitutions"])):
                        print("**"*20, antisite)
                        vac = cls(pc, 4*4*nsite, antisite, None, True, 20)
                        if MPRelaxSet(vac.defect_st).nelect % 2 == 0 and vac.defect_st.lattice.is_hexagonal():
                            os.makedirs(os.path.join(path, icsd_id, "pt_gp_" + point_gp), exist_ok=True)
                            vac.defect_st.to("POSCAR", "{}/{}-{}.vasp".format(
                                os.path.join(path, icsd_id, "pt_gp_" + point_gp),
                                vac.defect_st.formula, point_gp))
                            print(vac.defect_st.lattice.c)
                            if abs(vac.defect_st.lattice.c - 20) > 0.001:
                                return
                            wf_name = "{}:{}:{}".format(pc.composition.reduced_formula, icsd_id, point_gp)
                            fw = vac.hse_scf_no_relax_wf([0], [0.25, 0.25, 0], True, {"icsd":icsd_id,
                                                                                      "bulk_pointgp":point_gp}, icsd=icsd_id)
                            for fw_sole in fw.fws:
                                fws.append(fw_sole)
                    if fws:
                        wf = Workflow(fws, name=wf_name)
                        wf = add_namefile(wf)

                        lpad = LaunchPad.auto_load()
                        lpad.add_wf(wf)

    @classmethod
    def c2db_even_antisite_lg_bg(cls):
            col = VaspCalcDb.from_db_file("/home/tug03990/PycharmProjects/my_pycharm_projects/database/db_config/"
                                          "db_dk_local.json").collection
            tmdc = col.find(
                {"class":{"$in":["TMDC-T", "TMDC-H", "TMDC-T'"]},
                 "gap_hse":{"$gt":1},
                 "ehull":{"$lt":0.3},
                 "magstate": "NM",
                 # "spacegroup": "P-6m2",
                 "formula":{"$in": [
                                    # "Os2O4",
                                    # "Os2Te4",
                                    # "Ru2Te4"
                                    # "Ti2O4"
                                    "WSe2"
                                    ]}
                 }
            )
            path = '/gpfs/work/tug03990/c2db_TMD_bglg1/structures'

            for st in tmdc:
                pc = Structure.from_dict(st["structure"])
                for idx, el in enumerate(list(dict.fromkeys(pc.species))):
                    if el.is_metal:
                        cation = list(dict.fromkeys(pc.species))[idx].name
                    else:
                        anion = list(dict.fromkeys(pc.species))[idx].name
                defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
                fws = []
                for antisite in range(len(defect["substitutions"])):
                    print(cation, anion)
                    if "{}_on_{}".format(cation, anion) not in defect["substitutions"][antisite]["name"]:
                        continue
                    vac = cls(pc, 4*4*6, antisite, None, True, 20)
                    if MPRelaxSet(vac.defect_st).nelect % 2 == 0:
                        os.makedirs(os.path.join(path, st["formula"]), exist_ok=True)
                        vac.defect_st.to("POSCAR", "{}/{}.vasp".format(
                            os.path.join(path, st["formula"]), vac.defect_st.formula))
                        print(vac.defect_st.lattice.c)
                        if abs(vac.defect_st.lattice.c - 20) > 0.001:
                            print("vacuum thickness reset wrong!")
                            return
                        wf_name = "{}:{}".format(vac.defect_st.formula, st["spacegroup"])
                        for spin in [0, 2]:
                            fw = vac.hse_scf_no_relax_wf([0], True, True, {"antisite": pc.species[1].name},
                                                         icsd=st["formula"], nupdown=spin, encut=520)
                            for fw_sole in fw.fws:
                                fws.append(fw_sole)
                if fws:
                    wf = Workflow(fws, name=wf_name)
                    wf = add_namefile(wf)

                    lpad = LaunchPad.auto_load()
                    lpad.add_wf(wf)


class ZPLWF:
    def __init__(self, prev_calc_dir, spin_config):
        self.prev_calc_dir = prev_calc_dir
        structure = MPHSEBSSet.from_prev_calc(self.prev_calc_dir).structure

        if "magmom" in structure.site_properties.keys():
            structure.remove_site_property("magmom")

        self.structure = structure
        self.nelect = MPHSEBSSet.from_prev_calc(self.prev_calc_dir).nelect
        self.spin_config = spin_config
        self.encut = 1.3*max([potcar.enmax for potcar in MPHSERelaxSet(self.structure).potcar])

    def set_selective_sites(self, center_n, distance):
        selective_dyn = []
        for i in self.structure.get_sites_in_sphere(self.structure[center_n].coords, distance, include_image=True):
            for idx, site in enumerate(self.structure.sites):
                if np.allclose(i[0].frac_coords-i[2], site.frac_coords):
                    selective_dyn.append(idx)
        return selective_dyn

    def wf(self, task, charge, up_occupation, down_occupation, nbands, gamma_only=False, selective_dyn=None,
           specific_structure=None, nonsoc_prev_dir=None):

        user_kpoints_settings = Kpoints.gamma_automatic() if gamma_only else Kpoints.from_dict(
            {
                'comment': 'Automatic kpoint scheme',
                'nkpoints': 1,
                'generation_style': 'Reciprocal',
                'kpoints': [[0.25, 0.25, 0.0]],
                'usershift': (0, 0, 0),
                'kpts_weights': [1.0],
                'coord_type': None,
                'labels': ['None'],
                'tet_number': 0,
                'tet_weight': 0,
                'tet_connections': None,
                '@module': 'pymatgen.io.vasp.inputs',
                '@class': 'Kpoints'
            }
        )

        hse_incar_part = {
            # "AMIX": 0.2,
            # "AMIX_MAG": 0.8,
            # "BMIX": 0.0001,
            # "BMIX_MAG": 0.0001,
        }

        uis_static = {
            "user_incar_settings":
                {
                    "ENCUT": self.encut,
                    "ICHARG": 0,
                    "EDIFF": 1E-5,
                    "LCHARG": False,
                    "LWAVE": False,
                    "ISTART": 1,
                    "NELM": 100,
                    "NELECT": self.nelect-charge,
                    "NSW": 0,
                },
            "user_kpoints_settings": user_kpoints_settings
        }

        uis_relax = {
            "user_incar_settings":
                {
                    "ENCUT": self.encut,
                    # "ICHARG": 0,
                    "ISIF": 2,
                    "EDIFF": 1E-5,
                    "EDIFFG": -0.01,
                    "LCHARG": False,
                    "LWAVE": False,
                    "ISTART": 1,
                    "NELM": 250,
                    "NELECT": self.nelect - charge,
                    "NSW": 240,
                },
            "user_kpoints_settings": user_kpoints_settings
        }

        fws = []

        uis = copy.deepcopy(uis_static)
        uis["user_incar_settings"].update(hse_incar_part)
        cdft_B = JHSEcDFTFW(
            structure=self.structure,
            up_occupation=up_occupation,
            down_occupation=down_occupation,
            nbands=nbands,
            job_type="normal",
            prev_calc_dir=self.prev_calc_dir,
            specific_structure=specific_structure,
            vasp_input_set_params=uis,
            selective_dynamics=None,
            name="CDFT-B-HSE_scf",
            vasptodb_kwargs={
                "additional_fields": {
                    "task_type": "JHSEcDFTFW"
                },
                "parse_eigenvalues": False,
                "parse_dos": False
            }
        )

        uis = copy.deepcopy(uis_relax)
        uis["user_incar_settings"].update(hse_incar_part)
        cdft_C = JHSEcDFTFW(
            structure=self.structure,
            up_occupation=up_occupation,
            down_occupation=down_occupation,
            nbands=nbands,
            job_type="normal",
            vasp_input_set_params=uis,
            prev_calc_dir=self.prev_calc_dir,
            specific_structure=specific_structure,
            selective_dynamics=selective_dyn,
            name="CDFT-C-HSE_relax",
            vasptodb_kwargs={
                "additional_fields": {
                    "task_type": "JHSEcDFTFW"
                },
                "parse_eigenvalues": False,
                "parse_dos": False,
            }
        )

        uis = copy.deepcopy(uis_static)
        uis["user_incar_settings"].update(hse_incar_part)
        if self.spin_config:
            uis["user_incar_settings"].update({"NUPDOWN": 0 if self.spin_config == "singlet" else 2})

        parent = cdft_C
        if task == "D":
            parent = None
        cdft_D = JHSEStaticFW(
            structure=self.structure,
            parents=parent,
            vasp_input_set_params=uis,
            prev_calc_dir=nonsoc_prev_dir,
            name="CDFT-D-HSE_scf",
            vasptodb_kwargs={
                "additional_fields": {
                    "task_type": "JHSEStaticFW"
                },
                "parse_eigenvalues": False,
                "parse_dos": False,
            }
        )

        if task == "B":
            fws.append(cdft_B)
        elif task == "C":
            fws.append(cdft_C)
        elif task == "D":
            fws.append(cdft_D)
        elif task == "B-C":
            fws.append(cdft_B)
            fws.append(cdft_C)
        elif task == "C-D":
            fws.append(cdft_C)
            fws.append(cdft_D)
        elif task == "B-C-D":
            fws.append(cdft_B)
            fws.append(cdft_C)
            fws.append(cdft_D)

        wf = Workflow(fws, name="{}:{}:{}".format(self.structure.composition.reduced_formula,
                                                 self.spin_config, "{}CDFT".format(task)))
        wf = add_namefile(wf)

        try:
            vprun = Vasprun(os.path.join(self.prev_calc_dir, "vasprun.xml"))
        except FileNotFoundError:
            vprun = Vasprun(os.path.join(self.prev_calc_dir, "vasprun.xml.gz"))

        wf = add_additional_fields_to_taskdocs(
            wf, {
                "source":
                    {
                        "structure": self.structure.composition.reduced_formula,
                        "spin_config": self.spin_config,
                        "prev_path": self.prev_calc_dir,
                        "total_energy": vprun.final_energy
                    }
            }
        )
        return wf


class Ef_div:

    def __init__(self, charge, thickness):
        self.structure = Structure.from_file('/gpfs/work/tug03990/Ef_div_problem/bn66.vasp')
        if thickness:
            ase_atom_obj = AseAtomsAdaptor.get_atoms(self.structure)
            ase_atom_obj.center(vacuum=thickness / 2, axis=2)
            self.structure = AseAtomsAdaptor.get_structure(ase_atom_obj)
        self.charge = charge
        self.lpad = LaunchPad.auto_load()

    def wf(self):
        lpad = LaunchPad.auto_load()
        uis = {
            "ENCUT": 400,
            "ISIF": 2,
            "EDIFFG": -0.02,
            "EDIFF": 1E-4,
            "ISPIN": 1
        }
        vis_relax = MPRelaxSet(self.structure, force_gamma=False)
        v = vis_relax.as_dict()
        v.update({"user_kpoints_settings": Kpoints.monkhorst_automatic((3,3,1)), "user_incar_settings": uis})
        vis_relax = vis_relax.__class__.from_dict(v)

        uis = {
            "user_incar_settings": {
                "ENCUT": 400,
                "EDIFF": 1E-5,
                "ISPIN": 1
            },
            "user_kpoints_settings": Kpoints.monkhorst_automatic((3, 3, 1)).as_dict()
        }

        fws = []
        fws.append(OptimizeFW(structure=self.structure, vasp_input_set=vis_relax,
                             name="{}-{}-{}-PBE_relax-{}".format(
                                 LaunchPad.auto_load().get_fw_ids()[-1]+2,
                                 self.structure.composition.reduced_formula, self.charge,
                                 round(self.structure.lattice.c, 2))))
        fws.append(StaticFW(
            structure=self.structure, vasp_input_set_params=uis,
            name="{}-{}-{}-PBE_static".format(LaunchPad.auto_load().get_fw_ids()[-1]+1,
                                              self.structure.composition.reduced_formula,
                                              self.charge, round(self.structure.lattice.c, 2)),
            parents=fws[0],
            vasptodb_kwargs={
                "parse_dos": True,
                "parse_eigenvalues": True,
                "additional_fields": {"fw_id": LaunchPad.auto_load().get_fw_ids()[-1]+1}}
        ))
        wf = Workflow(fws, name="{}:{}:{}".format(self.structure.composition.reduced_formula,
                                                  self.charge, round(self.structure.lattice.c, 2)))
        lpad.add_wf(wf)


class SPWflows:
    @classmethod
    def MX2_bandgap_hse(cls):
        """
        using existing wf from atomate "wf_bandstructure_plus_hse"
        :return:
        """
        from atomate.vasp.workflows.presets.core import get_wf
        from atomate.vasp.powerups import add_modify_incar, add_additional_fields_to_taskdocs
        # col = VaspCalcDb.from_db_file("/home/tug03990/PycharmProjects/my_pycharm_projects/database/db_config/"
        #                               "db_dk_local.json").collection
        # mx2s = col.find(
        #     {"class": {"$in": ["TMDC-T", "TMDC-H", "TMDC-T'"]},
        #      "gap_hse": {"$gt": 1},
        #      "ehull": {"$lt": 0.3},
        #      "magstate": "NM",
        #      # "formula": "MoS2"
        #      # "spacegroup": "P-6m2",
        #      "formula": {"$in": ["MoS2", "MoSe2", "WSe2", "WS2", "WTe2", "MoTe2"]}
        #      }
        # )

        col = VaspCalcDb.from_db_file("/home/tug03990/PycharmProjects/my_pycharm_projects/database/db_config/"
                                      "db_mx2_antisite_pc.json").collection
        # mx2s = col.find({"task_id":{"$in":[3083, 3091, 3093, 3102, 3097, 3094]}})
        mx2s = col.find({"task_id":3083})
        # lpad = LaunchPad.from_file("/home/tug03990/config/my_launchpad.efrc.yaml")
        lpad = LaunchPad.auto_load()
        for mx2 in mx2s:
            pc = Structure.from_dict(mx2["output"]["structure"])
            pc.remove_site_property("magmom")
            # wf = wf_bandstructure_plus_hse(pc, gap_only=True)
            wf = get_wf(pc, '/gpfs/work/tug03990/mx2_antisite_basic_bandgap/bandstructure_hsegap.yaml')
            wf = add_modify_incar(wf, {"incar_update":{"NCORE": 4}})
            wf = add_modify_incar(wf, {"incar_update":{"EDIFF":1E-7}}, fw_name_constraint="static")
            wf = add_modify_incar(wf, {"incar_update":{"AEXX": 0.25, "LVHAR":True, "EDIFF":1E-7, "KPAR":4}},
                                  fw_name_constraint="hse gap")
            wf = add_modify_incar(wf, {"incar_update": {"NSW":0, "KPAR":4}}, fw_name_constraint="optimization")
            magmom = MPRelaxSet(pc).incar.get("MAGMOM", None)
            wf = add_modify_incar(wf, {"incar_update": {"MAGMOM": magmom}})
            wf = set_execution_options(wf, category="mx2_antisite_basic_bandgap")
            wf = add_additional_fields_to_taskdocs(wf, {"defect_type": "bulk", "wf": [fw.name for fw in wf.fws]})
            wf = set_queue_options(wf, "04:00:00", fw_name_constraint="optimization")
            lpad.add_wf(wf)


class Bilayer:
    def __init__(self, primitive_cell):
        self.pc = Structure.from_file(primitive_cell)
        for idx, el in enumerate(list(dict.fromkeys(self.pc.species))):
            if el.is_metal:
                cation = list(dict.fromkeys(self.pc.species))[idx].name
            else:
                anion = list(dict.fromkeys(self.pc.species))[idx].name
        defect = ChargedDefectsStructures(self.pc, antisites_flag=True).defects

        for antisite in range(len(defect["substitutions"])):
            print("{}_on_{}".format(cation, anion) in defect["substitutions"][antisite]["name"])
            print(round(defect["substitutions"][antisite]["unique_site"].c, 2))
            if "{}_on_{}".format(cation, anion) in defect["substitutions"][antisite]["name"] and \
                    round(defect["substitutions"][antisite]["unique_site"].c, 2) in [0.25,0.24, 0.15]:
                print(cation, anion)
                defect_info = DefectWF._defect_from_primitive_cell(
                    orig_st=Structure.from_file(primitive_cell),
                    antisite=True,
                    bulk=False,
                    type_vac=antisite,
                    natom=150
                )
                self.defect_st = defect_info[0]
                self.NN = defect_info[2]
                self.defect_st.to(
                    "poscar",
                    "/gpfs/work/tug03990/mx2_bilayer/bulk_st/{}.vasp".format(self.defect_st.formula.replace(" ", "-")))
                self.magmom = MPRelaxSet(self.defect_st).incar.get("MAGMOM", None)
            else:
                print("OWO")
                continue


    def bilayer_wf(self):
        vis_opt = MPRelaxSet(self.defect_st, user_incar_settings={
            "LCHARG":True,
            "ISIF":2,
            "EDIFF": 1E-4,
            "EDIFFG":-0.01,
            "ENCUT":320,
            "NSW": 150,
            "NCORE":4,
            "METAGGA": "SCAN"
        })
        opt = OptimizeFW(self.defect_st, name="PBE_relax", vasp_input_set=vis_opt)

        uis_static = {
            "ICHARG":1,
            "ISIF":2,
            "EDIFF": 1E-5,
            "ENCUT":320,
            "NEDOS": 9000,
            "EMAX": 10,
            "EMIN": -10,
            "LVHAR": True,
            "NELM": 120,
            "NCORE":4,
            "LAECHG": False
        }
        # static = StaticFW(self.defect_st,
        #                   name="bilayer:PBE_scf",
        #                   vasp_input_set_params={"user_incar_settings": uis_static},
        #                   vasptodb_kwargs=dict(parse_dos=True, parse_eigenvalues=True,
        #                                        additional_fields={"charge_state":0, "NN": self.NN}),
        #                   parents=opt)

        # uis_static_scan ={
        #     "ICHARG":1,
        #     "ISIF":2,
        #     "EDIFF": 1E-5,
        #     "ENCUT":320,
        #     "NEDOS": 9000,
        #     "EMAX": 10,
        #     "EMIN": -10,
        #     "LVHAR": True,
        #     "NELM": 120,
        #     "NCORE": 4,
        #     "METAGGA": "SCAN",
        #     "LAECHG": False
        # }
        uis_hse_scf = {
                "LEPSILON": False,
                "LVHAR": True,
                # "AMIX": 0.2,
                # "AMIX_MAG": 0.8,
                # "BMIX": 0.0001,
                # "BMIX_MAG": 0.0001,
                "EDIFF": 1.e-05,
                "ENCUT": 320,
                "ISMEAR": 0,
                "ICHARG": 1,
                "LWAVE": False,
                "LCHARG": True,
                "NSW": 0,
                "NUPDOWN": -1,
                "NELM": 150,
                "NCORE": 5
            }
        static_scan = StaticFW(self.defect_st,
                               name="HSE_scf",
                               vasp_input_set_params={"user_incar_settings": uis_hse_scf},
                               vasptodb_kwargs=dict(parse_dos=True, parse_eigenvalues=True,
                                                    additional_fields={"task_label":"HSE_scf", "charge_state":0, "NN": self.NN}),
                               parents=opt)

        fws = []
        for fw in [opt, static_scan]:
            fws.append(fw)
        wf = Workflow(fws, name="{}-bilayer".format(self.defect_st.formula))
        k = unfold.find_K_from_k([1/3, 1/3, 0], [[5, 0, 0], [0, 5, 0], [0, 0, 1]])[0]

        kpoints_kwarg = {
            'comment': 'bilayer',
            "style": "G",
            "num_kpts": 0,
            'kpts': [[1,1,1]],
            'kpts_weights': None,
            'kpts_shift': (0, 0, 0),
            'coord_type': None,
            'labels': None,
            'tet_number': 0,
            'tet_weight': 0,
            'tet_connections': None
        }

        wf = add_modify_kpoints(wf, {"kpoints_update": kpoints_kwarg})
        wf = add_modify_incar(wf, {"incar_update":{"NCORE":4, "NSW":100}}, "PBE_relax")
        wf = add_modify_incar(wf, {"incar_update":{"NCORE":4, "LWAVE": True, "LCHARG":False,
                                                   "AEXX":0.25}}, "HSE_scf")
        wf = set_queue_options(wf, "10:00:00", fw_name_constraint="HSE_scf")
        wf = set_queue_options(wf, "10:00:00", fw_name_constraint="PBE_relax")
        wf = set_execution_options(wf, category="mx2_antisite_bilayer")
        lpad = LaunchPad.from_file("/home/tug03990/config/my_launchpad.efrc.yaml")
        wf_singlet = add_modify_incar(
            wf,
            {"incar_update": {"MAGMOM": self.magmom, "ISMEAR":0, "SIGMA":0.05, "NUPDOWN":0}}
        )
        lpad.add_wf(wf_singlet)
        wf_triplet = add_modify_incar(
            wf,
            {"incar_update": {"MAGMOM": self.magmom, "ISMEAR":0, "SIGMA":0.05, "NUPDOWN":2}}
        )
        lpad.add_wf(wf_triplet)

    @classmethod
    def wf(cls):
        for st in glob("/gpfs/work/tug03990/mx2_bilayer/*.vasp"):
            if "Te4-Mo2.vasp" not in st:
                print(st)
                cls(st).bilayer_wf()

class Sandwich:
    def __init__(self, structure):
        if "magmom" in structure.site_properties.keys():
            structure.remove_site_property("magmom")
        self.structure = structure

        self.encut = 1.3*max([potcar.enmax for potcar in MPHSERelaxSet(structure).potcar])

    def wf(self):
        # hse_opt = JHSERelaxFW(self.structure, name="HSE-vdw_relax")
        # hse_opt = JScanOptimizeFW(self.structure, name="SCAN_rvv10_relax",
        #                           job_type="normal", max_force_threshold=False)
        hse_opt = JOptimizeFW(self.structure, job_type="normal", name="PBE_VDW_relax",
                              vasptodb_kwargs={"parse_eigenvalues": False, "parse_dos": False}
                              )

        hse_scf = JHSEStaticFW(
            self.structure,
            vasptodb_kwargs=dict(additional_fields={
                "charge_state":0,
                "pc_from": "mx2_antisite_pc/4237",
            }),
            name="HSE_scf",
            parents=hse_opt
        )

        # vis_kpt = MPRelaxSet(self.structure, vdw="optPBE", force_gamma=True).kpoints.as_dict()
        # vis_kpt = Kpoints.gamma_automatic((2,2,1)).as_dict()
        # vis_kpt.pop('@module')
        # vis_kpt.pop('@class')

        fws = [hse_opt, hse_scf]

        wf = Workflow(fws, name="{}:optPBE-vdw".format(self.structure.formula))
        wf = add_modify_incar(
            wf,
            {
                "incar_update":{
                    'LUSE_VDW': True,
                    'AGGAC': 0.0,
                    'GGA': 'Or',
                    "LASPH": True,
                    # "BPARAM": 15.7,#rvv10
                    "ISMEAR":0,
                    "SIGMA":0.05,
                }
            },
            hse_opt.name
        )

        wf = add_modify_incar(wf, {"incar_update":{
            "EDIFF":1E-4,
            "EDIFFG":-0.02,
            "LCHARG":False,
            "LWAVE":False
        }}, hse_opt.name)

        hse_scf_incar = {
            "EMAX":10,
            "EMIN":-10,
            "NEDOS":9000,
            "LCHARG":False,
            "LWAVE":False,
            "ISMEAR":0,
            "SIGMA":0.05,
            "LASPH":True,
            "LVHAR":True
        }
        wf = add_modify_incar(wf, {"incar_update": hse_scf_incar}, hse_scf.name)

        wf = add_modify_incar(wf, {"incar_update": {"ENCUT": self.encut}})

        wf = add_additional_fields_to_taskdocs(wf, {"fws": [fw.name for fw in wf.fws]})
        wf = add_modify_incar(wf)

        return wf

class Cluster:
    def __init__(self, structure):
        hydrogen_site = [i for i in range(17, 17+13)]
        a = []
        for i in range(len(structure.sites)):
            if i in hydrogen_site:
                a.append([True,True,True])
            else:
                a.append([False,False,False])
        structure.add_site_property("selective_dynamics", a)
        self.structure = structure

    def wf(self):
        lpad = LaunchPad.auto_load()
        opt = OptimizeFW(self.structure, name="HSE_relax", job_type="normal", max_force_threshold=False)
        hse_relax_incar = MPHSERelaxSet(self.structure).incar
        hse_relax_incar.update(
            {"ISIF":2,
             "IBRION":2,
             "EDIFF":1E-4,
             "EDIFFG":-0.01,
             "ENCUT":520,
             "LCHARG":False,
             "LWAVE":False,
             "LASPH": True
             }
        )

        fws = [opt]
        wf = Workflow(fws, name="{}:opt_cluster".format(self.structure.formula))

        wf = add_modify_incar(wf)
        wf = add_modify_incar(wf, {"incar_update": hse_relax_incar}, "relax")
        vis_kpt = Kpoints.gamma_automatic().as_dict()
        vis_kpt.pop('@module')
        vis_kpt.pop('@class')
        wf = add_modify_kpoints(wf, {"kpoints_update":vis_kpt})
        wf = set_execution_options(wf, category="cluster")
        wf = preserve_fworker(wf)
        lpad.add_wf(wf)


class FormationEnergy:
    @classmethod
    def comp_wf(cls):
        from atomate.vasp.workflows.jcustom.wf_full import get_wf_full_hse
        from pymatgen.ext.matproj import MPRester
        from fireworks import LaunchPad
        from qubitPack.tool_box import get_encut

        category = "formation_energy"
        lpad = LaunchPad.from_file("/home/tug03990/config/project/antisiteQubit/{}/my_launchpad.yaml".format(category))
        for mp_id in ["mp-570481"]:
            with MPRester() as mp:
                st = mp.get_structure_by_material_id(mp_id, conventional_unit_cell=True)
                print(st)
                wf = get_wf_full_hse(
                    structure=st,
                    charge_states=[0],
                    gamma_only=None,
                    gamma_mesh=True,
                    nupdowns=[-1],
                    task="hse_scf",
                    vasptodb={"category": category, "mp_id": mp_id},
                    wf_addition_name="",
                    task_arg={"parse_dos": False}
                )

                wf = add_modify_incar(
                    wf,
                    {"incar_update":{"ENCUT":get_encut(st), "ISIF":3, "ISMEAR": -5}},
                )
                wf = add_modify_incar(wf)
                # wf = set_queue_options(wf, "00:30:00", fw_name_constraint="HSE_scf")
                # wf = set_queue_options(wf, "00:30:00", fw_name_constraint="HSE_relax")
                # related to directory
                wf = set_execution_options(wf, category=category, fworker_name="efrc")
                wf = preserve_fworker(wf)
                lpad.add_wf(wf)

    @classmethod
    def calc_Ef(cls):
        from pymatgen.analysis.defects.core import DefectEntry
        from pymatgen.analysis.defects.core import Vacancy, Substitution
        from pymatgen.core.structure import PeriodicSite

        col = get_db("owls", "mx2_antisite_pc").collection
        # 4229: S-W, 4239: Se-W, 4236: Te-W, 4237:Mo-S, 4238: Mo-Se, 4235:Mo-Te
        pc = col.find_one({"task_id":{"$in":[4229]}})
        pc = Structure.from_dict(pc["output"]["structure"])
        col = get_db("antisiteQubit", "formation_energy", port=1234).collection
        host = col.find_one({"defect_type": "host"})
        defect = col.find_one({"defect_entry.name": "as_1_S_on_W"})
        defect_energy = defect["output"]["energy"]
        bulk_energy = host["output"]["energy"]
        parameters = {"vbm": -5.890, "band_gap": 2.333, "scaling_matrix": [5,5,1]}
        corrections = {} #if corrections exist, then one can import them as a dictionary to the DefectEntry
        uncorrected_energy = defect_energy - bulk_energy
        print(pc)
        d = Substitution(pc, PeriodicSite("S", pc.sites[0].frac_coords, pc.lattice))
        dentry = DefectEntry(d, uncorrected_energy, corrections=corrections, parameters=parameters)
        print("Defect Entry created for {} chg {} with uncorrected energy {}".format(dentry.name, dentry.charge,
                                                                                     dentry.uncorrected_energy))
        return dentry

    @classmethod
    def formation(cls, type):
        import pandas as pd
        col = get_db("antisiteQubit", "formation_energy").collection
        mx2_en = col.find_one({"task_id":767})["output"]["energy"]
        m_en = col.find_one({"task_id":760})["output"]["energy"]
        m_num_sites = col.find_one({"task_id":760})["nsites"]
        x_en = col.find_one({"task_id":763})["output"]["energy"]
        x_num_sites = col.find_one({"task_id":763})["nsites"]

        host_form_en = mx2_en/75*3 - m_en/m_num_sites - 2*x_en/x_num_sites
        print(type)
        print("host Ef:{:.3f}".format(host_form_en))

        def chempots(type=type):
            m_chempot, x_chempot= None, None
            if type == "M_rich":
                m_chempot = m_en/m_num_sites
                x_chempot = x_en/x_num_sites + 0.5*host_form_en
            else:
                x_chempot = x_en/x_num_sites
                m_chempot = m_en/m_num_sites + host_form_en
            return m_chempot, x_chempot

        def W_S():
            defect_form_en = col.find_one({"task_id": 771})["output"]["energy"] - mx2_en - (chempots()[0] - chempots()[1])
            return {"type": type, "defect": "W_S", "Ef": defect_form_en}

        def S_W():
            defect_form_en = col.find_one({"task_id": 768})["output"]["energy"] - mx2_en - (chempots()[1] - chempots()[0])
            return {"type": type, "defect": "S_W", "Ef": defect_form_en}

        def V_W():
            defect_form_en = col.find_one({"task_id": 762})["output"]["energy"] - mx2_en + chempots()[0]
            return {"type": type, "defect": "V_W", "Ef": defect_form_en}

        def V_S():
            defect_form_en = col.find_one({"task_id": 770})["output"]["energy"] - mx2_en + chempots()[1]
            # print(col.find_one({"task_id": 720})["output"]["energy"]/75*3 - mx2_en/75*3)
            # print(col.find_one({"task_id": 720})["output"]["energy"] - mx2_en)
            return {"type": type, "defect": "V_S", "Ef": defect_form_en}

        df = pd.DataFrame([W_S(), S_W(), V_W(), V_S()])
        df.to_clipboard()
        print(df)

class COHP:
    @classmethod
    def wf(cls):
        from atomate.vasp.fireworks.lobster import LobsterFW
        from qubitPack.tool_box import get_db
        from fireworks import LaunchPad, Workflow
        category = "cohp"
        source_static_tk_id = 786
        static = get_db("antisiteQubit", "move_z", port=12345)
        static = static.collection.find_one({"task_id": source_static_tk_id})
        static_path = static["calcs_reversed"][0]["dir_name"]
        fw = LobsterFW(structure=Structure.from_dict(static["input"]["structure"]), prev_calc_dir=static_path)
        wf = Workflow([fw], name=fw.name)
        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category, fworker_name="owls")
        wf = preserve_fworker(wf)
        print(wf)
        LPAD = LaunchPad.from_file(
            os.path.expanduser(os.path.join("~", "config/project/antisiteQubit/{}/"
                                                 "my_launchpad.yaml".format(category))))

        LPAD.add_wf(wf)

    @classmethod
    def analysis(cls):
        from pymatgen.electronic_structure.cohp import CompleteCohp
        from pymatgen.electronic_structure.plotter import CohpPlotter
        from pymatgen.electronic_structure.core import Orbital
        import os
        main_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/cohp/"
        main_path += "launcher_2021-05-25-21-22-33-177425" #"WTe2" -0.402076
        # main_path += "launcher_2021-05-25-21-13-08-771930" #"WS2" -2.023164
        # main_path += "launcher_2021-05-26-15-51-50-934362" #V_S in WS2
        # main_path += "launcher_2021-06-02-19-57-42-192796" #W_S -0.5; EF= -1.403767
        # main_path += "launcher_2021-06-02-19-57-49-331235" #W_Te 0.5


        os.chdir(main_path)
        # completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename="COHPCAR.lobster", structure_file="POSCAR")
        completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename="COOPCAR.lobster.gz", structure_file="POSCAR", are_coops=True)
        completecohp.efermi = -0.402076
        def summed(lobster_list, plot_label):
            labelist = lobster_list #["929", "961", "967"]
            cp = CohpPlotter(are_coops=True)
            # get a nicer plot label
            plotlabel = plot_label #"W75-W50+W55+W56 bonds"

            cp.add_cohp(plotlabel, completecohp.get_summed_cohp_by_label_list(
                label_list=labelist, divisor=1))
            x = cp.get_plot(integrated=False)
            x.ylim([-5, 5])
            x.show()

        def orbital_resolved(label, orbitals):
            #search for the number of the COHP you would like to plot in ICOHPLIST.lobster (the numbers in COHPCAR.lobster are different!)
            label = str(label)
            cp = CohpPlotter()
            #get a nicer plot label
            plotlabel=str(completecohp.bonds[label]['sites'][0].species_string)+'-'+str(completecohp.bonds[label]['sites'][1].species_string+"{}".format(orbitals))

            cp.add_cohp(plotlabel,completecohp.get_orbital_resolved_cohp(label=label, orbitals=orbitals))
            #check which COHP you are plotting

            print("This is a COHP between the following sites: "+str(completecohp.bonds[label]['sites'][0])+' and '+ str(completecohp.bonds[label]['sites'][1]))

            x = cp.get_plot(integrated=False)
            x.ylim([-5, 5])
            x.show()

        def sum_lobster_orbitals(lobster_list, orbitals, plot_label):
            cp = CohpPlotter()
            #"W75-W50+W55+W56 bonds"
            cp.add_cohp(plot_label, completecohp.get_summed_cohp_by_label_and_orbital_list(lobster_list, orbitals))
            x = cp.get_plot(integrated=False)
            x.ylim([-5, 5])
            x.show()

        # summed([str(i) for i in [13, 184, 215]], "W26-W1+W6+W7") #WS2
        summed([str(i) for i in [929, 961, 967]], "W75-W50+W55+W56") #WTe2
        # summed([str(i) for i in [3, 4, 173]], "W1+W6+W7") #Vs in WS2
        # orbital_resolved(13, orbitals=[[5, Orbital.dxy], [5, Orbital.dz2]])
        # antisite_o = [5, Orbital.dz2]
        # for o in [[5, Orbital.dx2], [5, Orbital.dxy], [5, Orbital.dz2], [5, Orbital.dxz], [5, Orbital.dyz]]:
        #     sum_lobster_orbitals([str(i) for i in [929, 961, 967]], [[o, antisite_o],
        #                                                             [o, antisite_o],
        #                                                             [o, antisite_o]], "{}-{}".format(o[1].name, antisite_o[1].name))
class MoveZ:
    def __init__(self):
        self.db = get_db("antisiteQubit", "move_z")
        # self.ws2 = [
        #     {"a1": 317, "e1": 325, "e2": 326, "a1e": 330, "e2e": 329, "e1e": 328, "dz": 0.1, "chemsys": "S-W"},
        #     {"a1": 317, "e1": 325, "e2": 326, "a1e": 330, "e2e": 329, "e1e": 328, "dz": 0, "chemsys": "S-W"},
        #     {"a1": 315, "e1": 325, "e2": 326, "a1e": 330, "e2e": 329, "e1e": 328, "dz": -0.1, "chemsys": "S-W"},
        #     {"a1": 315, "e1": 325, "e2": 326, "a1e": 330, "e2e": 329, "e1e": 328, "dz": -0.2, "chemsys": "S-W"},
        #     {"a1": 315, "e1": 325, "e2": 326, "a1e": 330, "e2e": 329, "e1e": 328, "dz": -0.25, "chemsys": "S-W"},
        #     {"a1": 315, "e1": 325, "e2": 326, "a1e": 330, "e2e": 329, "e1e": 328, "dz": -0.275, "chemsys": "S-W"},
        #     {"a1": 327, "e1": 313, "e2": 314, "a1e": 329, "e2e": 330, "e1e": 328, "dz": -0.3, "chemsys": "S-W"},
        #     {"a1": 327, "e1": 313, "e2": 314, "a1e": 328, "e2e": 330, "e1e": 329, "dz": -0.4, "chemsys": "S-W"},
        #     {"a1": 327, "e1": 313, "e2": 315, "a1e": 328, "e2e": 330, "e1e": 329, "dz": -0.5, "chemsys": "S-W"},
        # ]

        self.charged_ws2 = [
            {"a1": 316, "e1": 319, "e2": 320, "a1e": 330, "e2e": 329, "e1e": 328, "dz": 0.1, "chemsys": "S-W"},
            {"a1": 315, "e1": 319, "e2": 320, "a1e": 330, "e2e": 329, "e1e": 328, "dz": 0, "chemsys": "S-W"},
            {"a1": 315, "e1": 319, "e2": 320, "a1e": 330, "e2e": 329, "e1e": 328, "dz": -0.1, "chemsys": "S-W"},
            {"a1": 315, "e1": 319, "e2": 320, "a1e": 330, "e2e": 329, "e1e": 328, "dz": -0.2, "chemsys": "S-W"},
            {"a1": 315, "e1": 316, "e2": 317, "a1e": 330, "e2e": 329, "e1e": 328, "dz": -0.25, "chemsys": "S-W"},
            {"a1": 309, "e1": 307, "e2": 308, "a1e": 328, "e2e": 330, "e1e": 329, "dz": -0.3, "chemsys": "S-W"},
            {"a1": 309, "e1": 307, "e2": 308, "a1e": 328, "e2e": 330, "e1e": 329, "dz": -0.4, "chemsys": "S-W"},
            # {"a1": 327, "a1e": 328, "e1": 313, "e2": 315, "e2e": 330, "e1e": 329, "dz": -0.5, "chemsys": "S-W"},
        ]

        self.wte2 = [
            {"a1": 310, "a1e": 328, "e1": 258, "e2": 259, "e2e": 330, "e1e": 329, "dz": -0.1, "chemsys": "Te-W"},
            {"a1": 310, "a1e": 328, "e1": 269, "e2": 270, "e2e": 330, "e1e": 329, "dz": 0, "chemsys": "Te-W"},
            {"a1": 310, "a1e": 328, "e1": 269, "e2": 270, "e2e": 330, "e1e": 329, "dz": 0.1, "chemsys": "Te-W"},
            {"a1": 310, "a1e": 328, "e1": 269, "e2": 270, "e2e": 330, "e1e": 329, "dz": 0.2, "chemsys": "Te-W"},
            {"a1": 310, "a1e": 328, "e1": 269, "e2": 270, "e2e": 330, "e1e": 329, "dz": 0.3, "chemsys": "Te-W"},
            {"a1": 310, "a1e": 328, "e1": 269, "e2": 271, "e2e": 330, "e1e": 329, "dz": 0.4, "chemsys": "Te-W"},
            {"a1": 302, "a1e": 330, "e1": 325, "e2": 326, "e2e": 329, "e1e": 328, "dz": 0.5, "chemsys": "Te-W"},
        ]
    def get_sheet(self):
        import pandas as pd
        from matplotlib import pyplot as plt
        data = []
        for eigv in self.charged_ws2:
            e = self.db.collection.find_one({"chemsys": eigv["chemsys"], "dz": eigv["dz"], "charge_state":-1})
            locpot = max(e["calcs_reversed"][0]["output"]["locpot"]["2"])
            eig = self.db.get_eigenvals(e["task_id"])
            eig_a1e = eig["1"][0][eigv["a1e"]][0]-locpot
            eig_e1e = eig["1"][0][eigv["e1e"]][0]-locpot
            eig_e2e = eig["1"][0][eigv["e2e"]][0]-locpot
            eig_a1 = eig["1"][0][eigv["a1"]][0]-locpot
            eig_e2 = eig["1"][0][eigv["e2"]][0]-locpot
            eig_e1 = eig["1"][0][eigv["e1"]][0]-locpot
            data.append({"dz": eigv["dz"], "a1*":eig_a1e, "e1*":eig_e1e, "e2*": eig_e2e, "e1":eig_e1, "e2":eig_e2, "a1":eig_a1,
                         "a1_band": eigv["a1"], "a1e_band":eigv["a1e"], "e1_band":eigv["e1"], "e2_band":eigv["e2"], "e1e_band":eigv["e1e"],
                         "e2e_band": eigv["e2e"]
                         })
        df = pd.DataFrame(data).sort_values("dz")
        df.to_clipboard("\t")
        df = df.loc[:, ["dz", "a1", "a1*", "e1", "e1*", "e2", "e2*"]]
        fig = plt.figure(dpi=280)
        df.plot("dz", style=["r*-", "ro-", "b*--", "bo--", "b*-", "bo-"], ax=plt.gca())
        plt.xlabel("displacement in z-direction (Å)")
        plt.ylabel("level energy relative to vacuum (eV)")
        plt.legend(loc="upper center", frameon=True, ncol=6, bbox_to_anchor=(0.5, 1.1))
        plt.show()

def main():
    # DefectWF.MX2_formation_energy()
    # FormationEnergy.formation("M_rich")
    MoveZ().get_sheet()
    # COHP.analysis()

if __name__ == '__main__':
    main()


