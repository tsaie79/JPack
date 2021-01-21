#%%
from fireworks import Workflow
from fireworks import LaunchPad

from atomate.vasp.jpowerups import *
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
import glob


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
CATEGORY = "soc_standard_defect"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))

def ws2_anion_antisite(defect_type="substitutions", dtort=0): #1e-4
    col = VaspCalcDb.from_db_file("/home/tug03990/config/category/mx2_antisite_pc/db.json").collection
    # 4229: S-W, 4239: Se-W, 4236: Te-W, 4237:Mo-S, 4238: Mo-Se, 4235:Mo-Te
    #
    # col = VaspCalcDb.from_db_file(os.path.expanduser(os.path.join(
    #     "~", "config/project", "single_photon_emitter/pc", "db.json"))).collection

    mx2s = col.find({"task_id":{"$in":[4229]}}) #WS2

    geo_spec = {5* 5 * 3: [20]}
    for mx2 in mx2s:
        pc = Structure.from_dict(mx2["output"]["structure"])
        defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
        cation, anion = find_cation_anion(pc)

        for sub in range(len(defect[defect_type])):
            print(cation, anion)
            # cation vacancy
            if "as_1_{}_on_{}".format(cation,anion) not in defect[defect_type][sub]["name"]:
                continue
            for na, thicks in geo_spec.items():
                for thick in thicks:
                    se_antisite = GenDefect(
                        orig_st=pc,
                        defect_type=(defect_type, sub),
                        natom=na,
                        vacuum_thickness=thick,
                        distort=dtort,
                        sub_on_side=None
                    )
                    ###
                    scf_dir = "/home/tug03990/work/single_photon_emitter/cdft/A/WS2_enhenced_relax"
                    st = Structure.from_file(os.path.join(scf_dir, "CONTCAR"))


                    wf = get_wf_full_hse(
                        structure=st,
                        task_arg=None,
                        charge_states=[0],
                        gamma_only=False,
                        gamma_mesh=True,
                        scf_dos=True,
                        nupdowns=[-1],
                        task="hse_scf-hse_soc",
                        vasptodb={"category": CATEGORY, "NN": se_antisite.NN,
                                  "defect_entry": se_antisite.defect_entry,
                                  },
                        wf_addition_name="{}:{}".format(na, thick),
                        category=CATEGORY,
                    )
                    wf = add_modify_incar(wf, {"incar_update":{"LAECHG":False}})

                    wf = add_additional_fields_to_taskdocs(
                        wf,
                        {
                            "lattice_constant": "HSE",
                            "perturbed": se_antisite.distort,
                            # "pc_from": "single_photon_emitter/pc/HSE_scf/{}".format(mx2["task_id"])
                            "pc_from": "mx2_antisite_pc/HSE_scf/{}".format(mx2["task_id"]),
                            "sym": "C_3v_in_paper"
                        }
                    )

                    wf = scp_files(
                        wf,
                        "/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(CATEGORY),
                        fw_name_constraint=wf.fws[-1].name,
                    )
                    # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb",
                    #                     fw_name_constraint="HSE_scf")
                    wf = clear_to_db(wf, fw_name_constraint=wf.fws[0].name)

                    # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="PBE_relax")
                    # wf = set_queue_options(wf, "48:00:00", fw_name_constraint="HSE_relax")
                    wf = set_queue_options(wf, "32:00:00", fw_name_constraint="HSE_scf")
                    wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_soc")
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
from pymatgen.io.vasp import Poscar
import os, glob
from monty.serialization import loadfn
from projects.antisiteQubit.wf_defect import ZPLWF


CATEGORY = "cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))


def antistie_triplet_ZPL(sts_path, state="ground_state"): #6.5
    #efrc
    # anti_triplet = ZPLWF("/home/tug03990/work/antisiteQubit/perturbed/"
    #                      "block_2020-11-21-04-37-28-746929/launcher_2020-11-23-22-19-31-039641", "triplet")
    anti_triplet = ZPLWF("/home/tug03990/scratch/single_photon_emitter/cdft/WS2_C3V", "triplet")

    # sd_sites = anti_triplet.set_selective_sites(25, 5)
    # sd_sites.sort()
    # print(sd_sites)
    process = None
    if state == "ground_state":
        process = "D"
    elif state == "excited_state":
        process = "B"
    else:
        return


    info = loadfn(os.path.join(sts_path, "info.json".format(state)))

    for poscar in glob.glob(os.path.join(sts_path, "*010.vasp".format(state))):

        st = Poscar.from_file(poscar)
        """
        WS2_c3v: 
            FERDO = 328*1.0 147*0.0
            FERWE = 328*1.0 1*0.5 1*0.5 1*1.0 144*0.0
            NBANDS = 475
            ENCUT = 336
            NCORE = 4
            A = /home/tug03990/work/single_photon_emitter/cdft/A/WS2_C3V
        WS2_c3v_in_paper:
            ENCUT = 320
            FERDO = 224*1.0 151*0.0
            FERWE = 224*1.0 1*0.5 1*0.5 1*1.0 148*0.0
            NBANDS = 375
            NCORE = 4
            A = "/home/tug03990/work/single_photon_emitter/cdft/A/WS2_enhenced_relax"
        WS2_ch_in_paper:
            ENCUT = 320
            FERDO = 224*1.0 151*0.0
            FERWE = 224*1.0 1*0.5 1*0.5 1*1.0 148*0.0
            NBANDS = 375
            A = "/home/tug03990/work/single_photon_emitter/cdft/A/WS2_ch_in_paper"
        V_s_c3v:
            ENCUT = 336
            FERWE = 321*1.0 1*0 1*0.5 1*0.5 151*0.0
            FERDO = 322*1.0 153*0.0
            NBANDS = 475
            A='/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869/launcher_2021-01-16-02-46-39-229067'
        """

        wf = anti_triplet.wf(
            process, 0, up_occupation="328*1.0 1*0.5 1*0.5 1*1.0 144*0.0",
            down_occupation="328*1.0 147*0.0", nbands=475, gamma_only=True, selective_dyn=None, specific_structure=Poscar.from_file(poscar)
        )
        # wf = jmodify_to_soc(wf, nbands=575, structure=Structure.from_str(st.get_string(), "poscar"))

        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

        wf = add_modify_incar(wf, {"incar_update":{"ENCUT":336}})

        wf = add_modify_incar(wf)
        wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")
        wf = add_additional_fields_to_taskdocs(wf, {"cdft_info": "{}_C3V".format(state),
                                                    "interpolate_st_idx":poscar.split("/")[-1], "sts_info": info})

        wf = set_execution_options(wf, category=CATEGORY)
        wf = preserve_fworker(wf)
        LPAD.add_wf(wf)

antistie_triplet_ZPL("/home/tug03990/scratch/single_photon_emitter/cdft/interpolate_st/ground_state/extra")


#%%
CATEGORY = "soc_standard_defect"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))

scf_dir = "/home/tug03990/work/single_photon_emitter/cdft/A/WS2_C3V"

st = Structure.from_file(os.path.join(scf_dir, "POSCAR.gz"))

magmom = [[0,0,mag_z] for mag_z in MPRelaxSet(st).incar.get("MAGMOM", None)]

wf = get_wf_full_hse(
    structure=st,
    task_arg={"prev_calc_dir": scf_dir},
    charge_states=[0],
    gamma_only=False,
    gamma_mesh=True,
    scf_dos=True,
    nupdowns=[-1],
    task="hse_soc",
    vasptodb={"category": CATEGORY, "NN": se_antisite.NN,
              "defect_entry": se_antisite.defect_entry},
    wf_addition_name="{}:{}".format(na, thick),
    category=CATEGORY,
    )

wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)

wf = add_modify_incar(wf, {"incar_update":{"LCHARG":False, "LWAVE":False}})

wf = add_modify_incar(wf)

wf = clear_to_db(wf)

# wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")

wf = set_execution_options(wf, category=CATEGORY)
wf = preserve_fworker(wf)


LPAD.add_wf(wf)

#%% standard cdft
CATEGORY = "soc_cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))
"""
V_s_c3v:
    ENCUT = 336
    FERWE = 321*1.0 1*0 1*0.5 1*0.5 151*0.0
    FERDO = 322*1.0 153*0.0
    NBANDS = 475
    A='/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869/launcher_2021-01-16-02-46-39-229067'
"""

anti_triplet = ZPLWF("/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869/"
                     "launcher_2021-01-16-02-46-39-229067", "triplet")

wf = anti_triplet.wf(
    "D", 0, up_occupation="321*1.0 1*0 1*0.5 1*0.5 151*0.0",
    down_occupation="322*1.0 153*0.0", nbands=475, gamma_only=True, selective_dyn=None, specific_structure=None
)

wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
# wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
# wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

# wf = add_modify_incar(wf, {"incar_update":{"ENCUT":336}})

wf = add_modify_incar(wf)
wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")

wf = set_execution_options(wf, category=CATEGORY)
wf = preserve_fworker(wf)
LPAD.add_wf(wf)

#%% soc cdft
CATEGORY = "soc_cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))
"""
Ws_c3v:
    ENCUT = 336
    FERWE = 321*1.0 1*0 1*0.5 1*0.5 151*0.0
    FERDO = 322*1.0 153*0.0
    NBANDS = 475
    A='/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869/launcher_2021-01-16-02-46-39-229067'
"""

anti_triplet = ZPLWF("/home/tug03990/work/single_photon_emitter/soc_standard_defect/block_2021-01-17-04-46-01-108825/"
                     "launcher_2021-01-17-05-03-04-839545", None)

wf = anti_triplet.wf(
    "D", 0, up_occupation="657*1.0 1*0 1*1 291*0.0",
    down_occupation=None, nbands=950, gamma_only=True, selective_dyn=None, specific_structure=None
)

wf = jmodify_to_soc(wf, structure=anti_triplet.structure, nbands=950)

wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
# wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
# wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

# wf = add_modify_incar(wf, {"incar_update":{"ENCUT":336}})

wf = add_modify_incar(wf)
wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")
# wf = clear_to_db(wf)

wf = set_execution_options(wf, category=CATEGORY)
wf = preserve_fworker(wf)
LPAD.add_wf(wf)

#%% soc parabola
import glob
from projects.antisiteQubit.wf_defect import ZPLWF


# Ws_c3v_A_owls = "/home/tug03990/scratch/single_photon_emitter/soc_cdft/Ws_c3v/A/launcher_2021-01-17-05-03-04-839545"

CATEGORY = "soc_cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))
"""
soc_Ws_c3v:
    ENCUT = 336
    FERWE = 657*1.0 1*0 1*1 291*0.0
    NBANDS = 950
    A="/home/tug03990/work/single_photon_emitter/soc_standard_defect/block_2021-01-17-04-46-01-108825/launcher_2021-01-17-05-03-04-839545"
"""
def antistie_triplet_ZPL(sts_path, state="excited_state"): #6.5
    #efrc

    anti_triplet = ZPLWF("/home/tug03990/work/single_photon_emitter/soc_standard_defect/"
                         "block_2021-01-17-04-46-01-108825/launcher_2021-01-17-05-03-04-839545", None)

    # sd_sites = anti_triplet.set_selective_sites(25, 5)
    # sd_sites.sort()
    # print(sd_sites)
    process = None
    if state == "ground_state":
        process = "D"
    elif state == "excited_state":
        process = "B"
    else:
        return


    info = loadfn(os.path.join(sts_path, "{}/info.json".format(state)))

    for poscar in glob.glob(os.path.join(sts_path, "{}/*.vasp".format(state))):

        st = Structure.from_file(poscar)

        wf = anti_triplet.wf(
            process, 0, up_occupation="657*1.0 1*0 1*1 301*0.0",
            down_occupation=None, nbands=950, gamma_only=True, selective_dyn=None, specific_structure=Poscar.from_file(poscar)
        )
        wf = jmodify_to_soc(wf, nbands=950, structure=st)

        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

        wf = add_modify_incar(wf)
        # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")
        wf = add_additional_fields_to_taskdocs(wf, {"cdft_info": "{}_Ws_c3v_soc".format(state),
                                                    "interpolate_st_idx":poscar.split("/")[-1], "sts_info": info,
                                                    "poscar_idx": poscar.split("/")[-1]
                                                    })
        wf.name = wf.name+":{}".format(poscar.split("/")[-1])
        wf = set_execution_options(wf, category=CATEGORY)
        wf = preserve_fworker(wf)
        LPAD.add_wf(wf)

antistie_triplet_ZPL(sts_path="/home/tug03990/work/single_photon_emitter/soc_cdft/interpolate_st/")

#%% soc parabola point D (Ws_c3v_soc)
CATEGORY = "soc_cdft"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/my_launchpad.yaml".format(CATEGORY))))

scf_dir = "/home/tug03990/work/single_photon_emitter/cdft/A/WS2_C3V"

state = "ground_state"

sts_path = "/home/tug03990/work/single_photon_emitter/soc_cdft/interpolate_st/"

info = loadfn(os.path.join(sts_path, "{}/info.json".format(state)))

for poscar in glob.glob(os.path.join(sts_path, "{}/*.vasp".format(state))):

    st = Structure.from_file(poscar)

    wf = get_wf_full_hse(
        structure=st,
        task_arg={"prev_calc_dir": scf_dir},
        charge_states=[0],
        gamma_only=False,
        gamma_mesh=True,
        scf_dos=False,
        nupdowns=[-1],
        task="hse_soc",
        vasptodb={"cdft_info": "{}_Ws_c3v_soc".format(state),
                  "interpolate_st_idx":poscar.split("/")[-1], "sts_info": info,
                  "poscar_idx": float(poscar.split("/")[-1].split("_")[-1].split(".")[0])/10
                  },
        wf_addition_name="{}".format(poscar.split("/")[-1]),
        category=CATEGORY,
    )

    wf = write_PMGObjects(wf, dict(poscar=Poscar.from_file(poscar)))

    wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)

    wf = add_modify_incar(wf, {"incar_update":{"LCHARG":False, "LWAVE":False, "LAECHG":False}})

    wf = add_modify_incar(wf)

    wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")

    wf = set_execution_options(wf, category=CATEGORY)
    wf = preserve_fworker(wf)


    LPAD.add_wf(wf)