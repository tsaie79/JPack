#%% WF
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

#%% Build Structures
from pymatgen import Structure
import os

from mpinterfaces.interface import Ligand, Interface
from mpinterfaces.old_transformations import *
from mpinterfaces.utils import *

p = "/Users/jeng-yuantsai/Research/project/MgB2/model/new"

sic = Structure.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/SiC_6H/SiC_mp-7631_computed.cif")
sic.make_supercell([2,2,1])

sic_slab = Interface(sic, hkl=[0,0,1], min_thick=10, min_vac=25, primitive=False, from_ase=True)

# mgb2_slab = slab_from_file([0,0,1], "/Users/jeng-yuantsai/Research/project/MgB2/model/MgB2/MgB2_mp-763_computed.cif")
# mgb2_slab = slab_from_file([0,0,1], '/Users/jeng-yuantsai/Research/project/MgB2/model/new/mgb2_112.vasp')
mgb2 = Structure.from_file('/Users/jeng-yuantsai/Research/project/MgB2/model/new/mgb2_112.vasp')
mgb2_slab = Interface(mgb2, hkl=[0,0,1], min_thick=1, min_vac=25, primitive=False, from_ase=True, scell_nmax=1)



# mgb2_slab.to("poscar", os.path.join(p, "mgb2_slab.vasp"))
# sic_slab.to("poscar", os.path.join(p, "sic_slab.vasp"))

substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    sic_slab,
    mgb2_slab,
    max_area=200
)

substrate_slab_aligned.to("poscar", os.path.join(p, "subs.vasp"))
mat2d_slab_aligned.to("poscar", os.path.join(p, "mat2d.vasp"))

hetero_interfaces = generate_all_configs(mat2d_slab_aligned, substrate_slab_aligned,
                                         nlayers_substrate=1, nlayers_2d=1, seperation=3)
for idx, i in enumerate(hetero_interfaces):
    i.to("poscar", os.path.join(p,"structures", "{}.vasp".format(idx)))

#%%

from mpinterfaces.interface import Ligand, Interface
from pymatgen import Molecule, Structure

bk_st = Structure.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/SiC_6H/SiC_mp-7631_computed.cif")
bk_st.make_supercell([4,4,1], to_unit_cell=False)
hkl = [0,0,1]
min_thick = 10
min_vac = 25

surface_coverage = 0.12
mgb2_mol = Molecule.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/new/mgb2_cluster.xyz")
li = Ligand([mgb2_mol])
x,y,z = 0,0,3

adsorb_on_specie = "Si"
adatom_on_ligand = "Mg"

interface = Interface(
    bk_st,
    hkl=hkl,
    min_thick=min_thick,
    min_vac=min_vac,
    ligand=li,
    displacement=z,
    adatom_on_lig=adatom_on_ligand,
    adsorb_on_species=adsorb_on_specie,
    surface_coverage=surface_coverage,
    from_ase=True,
    primitive=False,
    supercell=[1,1,1]
)

interface.sort()
interface.create_interface()
# interface.to("poscar", "/Users/jeng-yuantsai/Research/project/MgB2/model/new/adsorb/0.vasp")
