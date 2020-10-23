from atomate.vasp.fireworks import OptimizeFW, StaticFW
from atomate.vasp.powerups import add_modify_incar, set_execution_options

from fireworks import Workflow, LaunchPad

from pymatgen import Structure

import os


def c_terminate_mgb2(cat="c_term"):
    lpad = LaunchPad.from_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/mgb2/{}/my_launchpad.yaml".format(cat)))

    init_st = Structure.from_file("/home/tug03990/config/project/mgb2/c_terminate_4by4.vasp")

    opt_fw = OptimizeFW(init_st, force_gamma=True, name="opt")
    scf_fw = StaticFW(init_st, parents=opt_fw, name="scf")

    fws = [opt_fw, scf_fw]
    wf = Workflow(fws, name="{}_c_term".format(init_st.formula))
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
            }
        },
        fw_name_constraint="scf"
    )

    wf = set_execution_options(wf, category=cat)
    lpad.add_wf(wf)


c_terminate_mgb2()