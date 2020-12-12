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

#%% make mgb2 cluster for vertical and horizontal configurations
from pymatgen import Structure, Molecule
import os

p = "/Users/jeng-yuantsai/Research/project/MgB2/model/new"
s = '/Users/jeng-yuantsai/Research/project/MgB2/model/new/building'

pc = Structure.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/MgB2/MgB2_mp-763_conventional_standard.cif")
sc_size = [4,4,2]
sc_size_name = "".join([str(i) for i in sc_size])
pc.make_supercell(sc_size, to_unit_cell=False)
Molecule.from_sites(pc.extract_cluster(pc.sites)).to("xyz", os.path.join(s, "mgb2_cluster",
                                                                         "mgb2_{}.xyz".format(sc_size_name)))



#%%
from pymatgen import Structure, Molecule
import os
from qubitPack.tool_box import *

mol = Molecule.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/new/building/mgb2_cluster/mgb2_322_pre_hydro.xyz")
bk_st = Structure.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/SiC_6H/SiC_mp-7631_computed.cif")
bk_st.make_supercell([4,4,1], to_unit_cell=False)
bk_st = modify_vacuum(bk_st, 40)

for i in mol.sites:
    coord = i.coords
    coord += [0,0,30-0.5]
    specie = i.specie
    bk_st.append(specie, coord,coords_are_cartesian=True)
    bk_st.sort()
bk_st.to("poscar", '/Users/jeng-yuantsai/Research/project/MgB2/model/new/building/mgb2_cluster/poscar.vasp')

# for idx in [138,139,140, 111,112,113, 135, 136,
#             137,108,109,110,132,133,134,
#             114,115,116, 141,142,143, 123,124,125]:
for idx in [134,135,110,111,132,133,108,109,130,131,112,113,136,137,118,119,142,143,124,125]:
    bk_st.replace(idx, "H")
bk_st.sort()
bk_st.to("poscar", '/Users/jeng-yuantsai/Research/project/MgB2/model/new/building/mgb2_cluster/poscar.vasp')

#%%
from pymatgen import Structure, Molecule, Element, Lattice
import os

hor = Structure.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/new/building/c_term/hor.vasp")
vert = Structure.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/new/building/c_term/vert.vasp")

want_site = []
for site in hor.extract_cluster(hor.sites):
    if not (str(site.specie) == "C" or str(site.specie) == "Si"):
        want_site.append(site)
    else:
        continue

mol = Molecule.from_sites(want_site)
mol.get_
st = Structure.from_sites(want_site)
st.modify_lattice(Lattice([[15, 0 ,0], [0,15,0], [0,0,15]]))
st.clear()
for s in want_site:
    st.append(s.specie, s.coords+[5.5, 1.5, 5.5], coords_are_cartesian=True)
st.to("poscar", "/Users/jeng-yuantsai/Research/project/MgB2/model/new/building/hydrogenate_mgb2/hydro_hor_box.vasp")
#%%
from pymatgen import Structure, Molecule, Element, Lattice
import os
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.transformations.advanced_transformations import SlabTransformation


mol = Molecule.from_file("/Users/jeng-yuantsai/"
                         "Research/project/MgB2/model/new/building/"
                         "hydrogenate_mgb2/hydro_hor.xyz")
slab = Structure.from_file("/Users/jeng-yuantsai/Research/project/MgB2/model/SiC_6H/SiC_mp-7631_computed.cif")
slab.make_supercell([4,4,1])
slab = SlabTransformation([0,0,1], 15, 30, primitive=False).apply_transformation(slab)

ads_slabs = AdsorbateSiteFinder(slab).generate_adsorption_structures(mol, find_args={"put_inside":True,
                                                                                     "distance":3,
                                                                                     "positions": ['bridge']
                                                                                     })
for n, ads_slab in enumerate(ads_slabs):
    ads_slab.to("poscar", "/Users/jeng-yuantsai/Research/project/MgB2/model/new/building/slabs_pmg/{}.vasp".format(n+1))

#%%
from pymatgen import Structure, Molecule, Element, Lattice
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.transformations.advanced_transformations import SlabTransformation
from atomate.vasp.workflows.base.adsorption import get_wf_slab, get_slab_trans_params, get_wf_molecules
from atomate.vasp.powerups import *
from fireworks import LaunchPad
import os


CATEGORY = "ads"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/mgb2/{}/my_launchpad.yaml".format(CATEGORY)))

mol = Molecule.from_file("/home/tug03990/work/mgb2/ads/hydro_hor.xyz")

bulk = Structure.from_file("/home/tug03990/work/mgb2/ads/SiC_mp-7631_computed.cif")
bulk.make_supercell([4,4,1])
slab = SlabTransformation([0,0,1], 15, 30, primitive=False).apply_transformation(bulk)


wf = get_wf_slab(slab,
                 include_bulk_opt=False,
                 adsorbates=[mol],
                 ads_structures_params={"find_args": {"put_inside":True, "distance":3, "positions": ['bridge']}},
                 add_molecules_in_box=False,
                 vasp_cmd=">>vasp_cmd<<",
                 db_file=">>db_file<<"
                 )

wf = Workflow([fw for fw in wf.fws[1:]], name=wf.name)
wf = add_additional_fields_to_taskdocs(wf, {"stack":"horizontal", "fws": [fw.name for fw in wf.fws]})
wf = set_execution_options(wf, category=CATEGORY)
wf = preserve_fworker(wf)

for fw in wf.fws:
    print(fw.name)

wf_mol = get_wf_molecules([mol], box_width=20, vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<", name="hor_mgb2")

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

wf = add_modify_kpoints(wf_mol, {"kpoints_update": kpoints_kwarg})
wf_mol = add_additional_fields_to_taskdocs(wf_mol, {"stack":"horizontal"})
wf_mol = set_execution_options(wf_mol, category=CATEGORY)
wf_mol = preserve_fworker(wf_mol)
LPAD.add_wf(wf_mol)

#%% fws for molecule mgb2
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen import Molecule

from atomate.vasp.powerups import *
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW

from fireworks import LaunchPad, Workflow
import os

CATEGORY = "ads"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/mgb2/{}/my_launchpad.yaml".format(CATEGORY)))

task_name = "hor"

mol = Molecule.from_file("/home/tug03990/work/mgb2/ads/hydro_{}.xyz".format(task_name))
mol_boxed = mol.get_boxed_structure(20, 20, 20, offset=[10,10,10])

opt_fw = OptimizeFW(mol_boxed, job_type="normal")
scf_fw = StaticFW(mol_boxed, parents=opt_fw)
fws = [opt_fw, scf_fw]
wf = Workflow(fws, name=task_name)

encut = 1.3*max([potcar.enmax for potcar in MPRelaxSet(mol_boxed).potcar])
encut = encut if encut <= 520 else 520

wf = add_modify_incar(wf, {"incar_update": {"ENCUT": encut, "ISIF":2,
                                            "EDIFFG": -0.05, "EDIFF": 1E-4, "LCHARG": False, "ISPIN":1}}, opt_fw.name)
wf = add_modify_incar(wf, {"incar_update": {"ENCUT": encut, "EDIFF": 1E-5, "LCHARG": False, "ISPIN":1}}, scf_fw.name)

wf = add_additional_fields_to_taskdocs(wf, {"stack":task_name})
wf = set_execution_options(wf, category=CATEGORY)

LPAD.add_wf(wf)

#%% fws slab
from pymatgen.transformations.advanced_transformations import SlabTransformation
from pymatgen.analysis.adsorption import AdsorbateSiteFinder

from pymatgen.io.vasp.sets import MPRelaxSet, MVLSlabSet
from pymatgen import Molecule
from pymatgen.io.vasp.inputs import Poscar

from atomate.vasp.powerups import *
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW

from fireworks import LaunchPad, Workflow
import os

CATEGORY = "ads"
LPAD = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/mgb2/{}/my_launchpad.yaml".format(CATEGORY)))

task_name = "vert"

mol = Molecule.from_file("/home/tug03990/work/mgb2/ads/hydro_{}.xyz".format(task_name))
slab = Structure.from_file("/home/tug03990/work/mgb2/ads/SiC_mp-7631_computed.cif")
slab.make_supercell([4,4,1])
slab = SlabTransformation([0,0,1], 15, 30, primitive=False).apply_transformation(slab)

ads_slabs = AdsorbateSiteFinder(slab, selective_dynamics=True).generate_adsorption_structures(mol, find_args={"put_inside":True,
                                                                                     "distance":3,
                                                                                     "positions": ['bridge']
                                                                                     })
ads_slab = ads_slabs[0]


opt_fw = OptimizeFW(ads_slab, vasp_input_set=MVLSlabSet(ads_slab, auto_dipole=True), job_type="normal")
scf_fw = StaticFW(ads_slab, parents=opt_fw)
fws = [opt_fw, scf_fw]
wf = Workflow(fws, name="{}:{}".format(task_name, ads_slab.formula))

encut = 1.3*max([potcar.enmax for potcar in MPRelaxSet(ads_slab).potcar])
encut = encut if encut <= 520 else 520


wf = add_modify_incar(wf, {"incar_update": {"ENCUT": encut,
                                            "EDIFFG": -0.05, "EDIFF": 1E-4, "LCHARG": False, "ISPIN":1}}, opt_fw.name)
wf = add_modify_incar(wf, {"incar_update": {"ENCUT": encut, "EDIFF": 1E-5, "LCHARG": False, "ISPIN":1}}, scf_fw.name)

wf = add_additional_fields_to_taskdocs(wf, {"stack":task_name})
wf = set_execution_options(wf, category=CATEGORY)

LPAD.add_wf(wf)


