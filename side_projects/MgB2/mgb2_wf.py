from atomate.vasp.fireworks import OptimizeFW, StaticFW
from atomate.vasp.fireworks.jcustom import JSelectiveOptFW
from atomate.vasp.powerups import add_modify_incar, set_execution_options, add_additional_fields_to_taskdocs, \
    add_modify_kpoints

from fireworks import Workflow, LaunchPad

from atomate.vasp.firetasks import SelectiveDynmaicPoscar

from pymatgen import Structure

import os, glob


def c_terminate_mgb2(cat="slab_dipol"):
    lpad = LaunchPad.from_file(
        os.path.join(
            os.path.expanduser("~"),
            "config/project/mgb2/slab/my_launchpad.yaml".format(cat)))

    for st in glob.glob("/home/tug03990/config/project/mgb2/slab/structure/mv/c_t*.vasp".format(cat)):
        init_st = Structure.from_file(st).get_sorted_structure()
        st_name = st.split("/")[-1].split(".")[0]

        two_layers_from_top = []
        for idx, site in enumerate(init_st.sites):
            if str(site.specie) == "Mg" or str(site.specie) == "B":
                two_layers_from_top.append(idx)

        opt_0 = JSelectiveOptFW(init_st, selective_dynamics=None,
                                 force_gamma=True, vasptodb_kwargs={"parse_eigenvalues":False}, name="{}_opt_0".format(st_name))

        two_layers_from_top = []
        for idx, site in enumerate(init_st.sites):
            if "c" in st_name and site.c > 0.15:
                # print(idx, site.c)
                two_layers_from_top.append(idx)
            elif "si" in st_name and site.c > 0.09:
                two_layers_from_top.append(idx)

        opt_1 = JSelectiveOptFW(init_st,
                                job_type="normal",
                                max_force_threshold=None,
                                selective_dynamics=None,
                                parents=None,
                                force_gamma=True,
                                vasptodb_kwargs={"parse_eigenvalues":False},
                                name="{}_opt_1".format(st_name))


        scf_fw = StaticFW(init_st, parents=opt_1, vasptodb_kwargs={"parse_eigenvalues":False}, name="{}_scf".format(st_name))

        fws = [opt_1, scf_fw]
        wf = Workflow(fws, name="{}:{}:no_sel".format(st_name, cat))
        wf = add_modify_incar(
            wf,
            {
                "incar_update": {
                    "ENCUT": 520,
                    "EDIFFG": 1E-3,
                    "EDIFF": 1E-4,
                    "LDIPOL": True,
                    # "DIPOL": "0.5 0.5 0.7",
                    "IDIPOL": 3,
                    "NSW": 150,
                    "ISPIN":1,
                    "ISIF": 2,
                    "LCHARG": False,
                    "LVHAR": True,
                }
            },
            fw_name_constraint="opt"
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

        # wf = add_modify_kpoints(wf, {"kpoints_update": kpoints((1,1,1))}, fw_name_constraint=opt_0.name)

        wf = add_modify_incar(
            wf,
            {
                "incar_update": {
                    "ENCUT": 520,
                    "EDIFF": 1E-5,
                    "LDIPOL": True,
                    "IDIPOL": 3,
                    # "DIPOL": "0.5 0.5 0.7",
                    "ISPIN":1,
                    "LVHAR": True,
                }
            },
            fw_name_constraint="scf"
        )
        wf = add_additional_fields_to_taskdocs(wf, {"cat": cat, "stack": st_name, "selective":None})
        wf = set_execution_options(wf, category=cat)
        print(st_name)
        print(wf.name)
        lpad.add_wf(wf)


c_terminate_mgb2()