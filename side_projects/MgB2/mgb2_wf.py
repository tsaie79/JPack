from atomate.vasp.fireworks import OptimizeFW, StaticFW
from atomate.vasp.powerups import add_modify_incar, set_execution_options

from fireworks import Workflow, LaunchPad

from pymatgen import Structure

import os, glob


def c_terminate_mgb2(cat="b2_mg_b2"):
    lpad = LaunchPad.from_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/mgb2/{}/my_launchpad.yaml".format(cat)))

    for st in glob.glob("/home/tug03990/config/project/mgb2/{}/structure/*.vasp".format(cat)):
        init_st = Structure.from_file(st)

        opt_fw = OptimizeFW(init_st, force_gamma=True, vasptodb_kwargs={"parse_eigenvalues":False}, name="opt")
        scf_fw = StaticFW(init_st, parents=opt_fw, vasptodb_kwargs={"parse_eigenvalues":False}, name="scf")

        fws = [opt_fw, scf_fw]
        wf = Workflow(fws, name="{}_{}_{}".format(init_st.formula, cat, st.split("/")[-1]))
        wf = add_modify_incar(
            wf,
            {
                "incar_update": {
                    "ENCUT": 520,
                    "EDIFFG": -0.01,
                    "EDIFF": 1E-4,
                    "LDIPOL": True,
                    "IDIPOL": 3,
                    "NSW": 150,
                    "ISPIN":1,
                    "ISIF": 2,
                    "LCHARG": False,
                }
            },
            fw_name_constraint="opt"
        )

        wf = add_modify_incar(
            wf,
            {
                "incar_update": {
                    "ENCUT": 520,
                    "EDIFF": 1E-5,
                    "LDIPOL": True,
                    "IDIPOL": 3,
                    "ISPIN":1,
                }
            },
            fw_name_constraint="scf"
        )

        wf = set_execution_options(wf, category=cat)
        print(wf.name)
        lpad.add_wf(wf)


c_terminate_mgb2()