#%%
from fireworks import Workflow
from fireworks import LaunchPad

from atomate.vasp.powerups import *
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
from atomate.vasp.jpowerups import scp_files


#%% calculate pristine WS2 and determine band gap
def hse_scf_wf():
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

def hse_band_structure():
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

#%% make standard defect calc
CATEGORY = "standard_defect"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))

def ws2_anion_antisite(defect_type="vacancies", dtort=0): #1e-4
    col = VaspCalcDb.from_db_file("/home/tug03990/config/category/mx2_antisite_pc/db.json").collection
    # 4229: S-W, 4239: Se-W, 4236: Te-W, 4237:Mo-S, 4238: Mo-Se, 4235:Mo-Te

    # col = VaspCalcDb.from_db_file(os.path.expanduser(os.path.join(
    #     "~", "config/project", "single_photon_emitter/pc", "db.json"))).collection

    mx2s = col.find({"task_id":{"$in":[4237]}}) #WS2

    geo_spec = {5* 5 * 3: [20]}
    aexx = 0.25
    for mx2 in mx2s:
        pc = Structure.from_dict(mx2["output"]["structure"])
        defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
        cation, anion = find_cation_anion(pc)

        for sub in range(len(defect[defect_type])):
            print(cation, anion)
            # cation vacancy
            if "vac_2_{}".format(anion) not in defect[defect_type][sub]["name"]:
                continue
            for na, thicks in geo_spec.items():
                for thick in thicks:
                    se_antisite = GenDefect(
                        orig_st=pc,
                        defect_type=(defect_type, sub),
                        natom=na,
                        vacuum_thickness=thick,
                        distort=dtort,
                        sub_on_side=[Element("Re")]
                    )

                    wf = get_wf_full_hse(
                        structure=se_antisite.defect_st,
                        charge_states=[0],
                        gamma_only=False,
                        gamma_mesh=True,
                        scf_dos=True,
                        nupdowns=[-1],
                        task="opt-hse_relax-hse_scf",
                        vasptodb={"category": CATEGORY, "NN": se_antisite.NN,
                                  "defect_entry": se_antisite.defect_entry},
                        wf_addition_name="{}:{}".format(na, thick),
                        category=CATEGORY,
                    )

                    wf = add_additional_fields_to_taskdocs(
                        wf,
                        {
                            "lattice_constant": "HSE",
                            "perturbed": se_antisite.distort,
                            # "pc_from": "single_photon_emitter/pc/HSE_scf/{}".format(mx2["task_id"])
                            "pc_from": "mx2_antisite_pc/HSE_scf/{}".format(mx2["task_id"])
                        }
                    )

                    # wf = add_modify_incar(wf, {"incar_update": {"LWAVE":False, "LCHARG": True}}, fw_name_constraint="HSE_scf")

                    wf = scp_files(
                        wf,
                        "/home/jengyuantsai/Research/projects/single_photon_emitter/standard_defect",
                        fw_name_constraint="HSE_scf",
                    )
                    wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb",
                                        fw_name_constraint="HSE_scf")

                    wf = set_queue_options(wf, "24:00:00", fw_name_constraint="PBE_relax")
                    wf = set_queue_options(wf, "48:00:00", fw_name_constraint="HSE_relax")
                    wf = set_queue_options(wf, "32:00:00", fw_name_constraint="HSE_scf")
                    wf = add_modify_incar(wf)
                    # related to directory
                    wf = preserve_fworker(wf)
                    wf.name = wf.name+":dx[{}]".format(se_antisite.distort)
                    print(wf.name)
                    LPAD.add_wf(wf)
ws2_anion_antisite()

#%% cdft
"""
Run CDFT on MoS2 w/ perturbation
"""
from atomate.vasp.powerups import *
from fireworks import LaunchPad
from projects.antisiteQubit.wf_defect import ZPLWF
import os

CATEGORY = "cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))


def antistie_triplet_ZPL(): #6.5
    #efrc
    # anti_triplet = ZPLWF("/home/tug03990/work/antisiteQubit/perturbed/"
    #                      "block_2020-11-21-04-37-28-746929/launcher_2020-11-23-22-19-31-039641", "triplet")
    anti_triplet = ZPLWF("/home/tug03990/work/single_photon_emitter/cdft/A/launcher_2021-01-13-04-09-57-563577", "triplet")

    # sd_sites = anti_triplet.set_selective_sites(25, 5)
    # sd_sites.sort()
    # print(sd_sites)

    for poscar in glob.glob("/home/tug03990/work/single_photon_emitter/cdft/interpolate_st/*"):
        st = Structure.from_file(poscar)
        wf = anti_triplet.wf(
            "B", 0, up_occupation="328*1.0 1*0.5 1*0.5 1*1.0 144*0.0",#"303*1.0 1*0 1*1.0 175*0",
            down_occupation="328*1.0 147*0.0", nbands=475, gamma_only=True, selective_dyn=None, specific_structure=st
        )

        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

        wf = add_modify_incar(wf)
        wf = clean_up_files(wf, files=["WAVECAR*"], task_name_constraint="VaspToDb")
        wf = add_additional_fields_to_taskdocs(wf, {"cdft_info": "C3V", "interpolate_st_idx":poscar})

        wf = set_execution_options(wf, category=CATEGORY)
        wf = preserve_fworker(wf)

        # if moving_sites:
        #     wf.name = wf.name + ":delta{:.2f}".format(len(moving_sites) / len(anti_triplet.structure.sites))
        print(wf.name)
        LPAD.add_wf(wf)



antistie_triplet_ZPL()

