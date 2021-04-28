#%%
import numpy as np
from pymatgen.io.vasp.sets import MPRelaxSet, MPMetalRelaxSet
from pymatgen.io.vasp.inputs import Kpoints
from atomate.vasp.fireworks.jcustom import JHSEStaticFW
from atomate.vasp.fireworks import OptimizeFW, HSEBSFW
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import cp_vdw_file
from fireworks import LaunchPad, Workflow
import os, glob
from monty.serialization import loadfn
from atomate.vasp.config import LOBSTER_CMD

CATEGORY = "thinner_slab"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/mag_insulator/{}/my_launchpad.yaml".format(CATEGORY))))

tci_hex_slabs = loadfn("/home/tug03990/config/project/mag_insulator/slab/structures/tci_hex_slabs.json")
ti_hex_slabs = loadfn("/home/tug03990/config/project/mag_insulator/slab/structures/ti_hex_slabs.json")
ti_hex_thinner_slabs = loadfn("/home/tug03990/config/project/mag_insulator/slab/structures/ti_hex_thinner_slabs.json")

slabs = ti_hex_slabs.copy()
slabs.update(tci_hex_slabs)


for icsd_id, v in ti_hex_thinner_slabs.items():
    if v["band_gap"] != 0 and v["magnetic_type"] == "NM" and icsd_id == "043512":
        ismear, sigma = 0, 0.01
        print(icsd_id)
        fws = []
        st = Structure.from_str(v['cifs.conventional_standard'], fmt="cif")
        encut = 1.3*max([potcar.enmax for potcar in MPRelaxSet(st).potcar])

        override_default_vasp_params = {
            "user_incar_settings": {"EDIFF":1E-4, "EDIFFG":-0.01, "ISIF":3, "SIGMA":sigma, "ISPIN":1,
                                    "ENCUT":encut, "LASPH":False, "ISMEAR":ismear,
                                    "LCHARG":False},
            # "user_kpoints_settings": Kpoints.automatic_density_by_vol(st, 6**3),
            "vdw": None}

        vis_set = MPMetalRelaxSet(st, force_gamma=True, **override_default_vasp_params)
        opt_vdw = OptimizeFW(st, vasp_input_set=vis_set, vasp_cmd=">>vasp_cmd<<")

        hse_static = JHSEStaticFW(
            st, vasp_input_set_params={
                "user_incar_settings": {"EDIFF":1E-5, "ISMEAR":0, "SIGMA":sigma, "ENCUT":encut,
                                        "LASPH":False, "ISPIN":1, "NEDOS":9000, "EMAX":10, "EMIN":-10,
                                        "LVHAR":True, "LCHARG": False},
            },
            vasptodb_kwargs={"parse_chgcar":True},
            parents=opt_vdw
        )

        hse_bs = HSEBSFW(structure=st, parents=hse_static, mode="line",
                         input_set_overrides={"other_params": {"two_d_kpoints": True}})

        fws.append(opt_vdw)
        fws.append(hse_static)
        fws.append(hse_bs)

        wf = Workflow(fws, name="{}:{}".format(st.formula, icsd_id))
        v.pop("cifs.conventional_standard")
        wf = add_additional_fields_to_taskdocs(wf, {"icsd_id": icsd_id, "mpid": v["task_id"], "mp_info": v})
        opt_vdw = cp_vdw_file(wf, fw_name_constraint=opt_vdw.name)

        wf = preserve_fworker(wf)
        wf = set_execution_options(wf, category=CATEGORY)
        wf = add_modify_incar(wf)

        wf = set_queue_options(wf, "32:00:00", fw_name_constraint=wf.fws[0].name)
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[2].name)

        LPAD.add_wf(wf)

#%% PBE opt+PBE SCF+HSE NSCF
import numpy as np
from pymatgen.io.vasp.sets import MPRelaxSet, MPMetalRelaxSet
from pymatgen.io.vasp.inputs import Kpoints
from atomate.vasp.fireworks.jcustom import JHSEStaticFW
from atomate.vasp.fireworks import OptimizeFW, StaticFW, HSEBSFW
from atomate.vasp.powerups import *
from atomate.vasp.jpowerups import cp_vdw_file
from fireworks import LaunchPad, Workflow
import os, glob
from monty.serialization import loadfn

CATEGORY = "slab"
LPAD = LaunchPad.from_file(
    os.path.expanduser(os.path.join("~", "config/project/mag_insulator/{}/my_launchpad.yaml".format(CATEGORY))))

# tci_hex_slabs = loadfn("/home/tug03990/config/project/mag_insulator/slab/structures/tci_hex_slabs.json")
# ti_hex_slabs = loadfn("/home/tug03990/config/project/mag_insulator/slab/structures/ti_hex_slabs.json")
# slabs = ti_hex_slabs.copy()
# slabs.update(tci_hex_slabs)

ti_hex_thinner_slabs = loadfn("/home/tug03990/config/project/mag_insulator/slab/structures/ti_hex_thinner_slabs.json")
tci_hex_thinner_slabs = loadfn("/home/tug03990/config/project/mag_insulator/slab/structures/tci_hex_thinner_slabs.json")

for icsd_id, v in tci_hex_thinner_slabs.items():
    if v["monolayer_scan_gap"] == 0 and v["magnetic_type"] == "NM" and icsd_id == "059801":
        ismear, sigma = 0, 0.05
        print(icsd_id)
        fws = []
        st = Structure.from_str(v["cifs.conventional_standard"], fmt="cif")
        encut = 1.3*max([potcar.enmax for potcar in MPRelaxSet(st).potcar])

        override_default_vasp_params = {
            "user_incar_settings": {
                "EDIFF":1E-4,
                "EDIFFG":-0.01,
                "ISIF":3,
                "SIGMA":sigma,
                "ISPIN":1,
                "ENCUT":encut,
                "LASPH":False,
                "ISMEAR":ismear,
                "LCHARG":False
            },
            # "user_kpoints_settings": Kpoints.automatic_density_by_vol(st, 6**3),
            "vdw": None}

        opt = OptimizeFW(st, vasp_cmd=">>vasp_cmd<<", override_default_vasp_params=override_default_vasp_params)

        uis_static = {
            "user_incar_settings": {
                "EDIFF":1E-5,
                "ISMEAR":ismear,
                "SIGMA":sigma,
                "ENCUT":encut,
                "LASPH":False,
                "ISPIN":1,
                "LVHAR":False,
                "LCHARG": True,
                "LWAVE":False
            }
        }
        static = StaticFW(st, parents=opt, vasp_input_set_params=uis_static,
                          vasptodb_kwargs={"parse_eigenvalues":False, "parse_dos":False})

        uis_hse_bs = {
            "EDIFF":1E-5,
            "ISMEAR":ismear,
            "SIGMA":sigma,
            "ENCUT":encut,
            "LASPH":False,
            "ISPIN":1,
            "LVHAR":True,
            "LCHARG": False,
            "NELM": 250,
        }

        hse_bs = HSEBSFW(structure=st, parents=static, mode="line", cp_file_from_prev="CHGCAR",
                         input_set_overrides={"other_params": {"two_d_kpoints": True, "user_incar_settings": uis_hse_bs}})

        fws.append(opt)
        fws.append(static)
        fws.append(hse_bs)

        wf = Workflow(fws, name="{}:{}".format(st.formula, icsd_id))
        v.pop("cifs.conventional_standard")
        wf = add_additional_fields_to_taskdocs(wf, {"icsd_id": icsd_id, "mpid": v["task_id"], "mp_info": v})
        opt_vdw = cp_vdw_file(wf, fw_name_constraint=opt.name)

        wf = preserve_fworker(wf)
        wf = set_execution_options(wf, category=CATEGORY)
        wf = add_modify_incar(wf)

        wf = set_queue_options(wf, "32:00:00", fw_name_constraint=wf.fws[0].name)
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[1].name)
        wf = set_queue_options(wf, "32:00:00", fw_name_constraint=wf.fws[2].name)

        LPAD.add_wf(wf)
