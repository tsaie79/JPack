from pymatgen import Structure
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.analysis.adsorption import plot_slab
from pymatgen.ext.matproj import MPRester
from matplotlib import pyplot as plt

import os
# os.chdir("/Users/jeng-yuantsai/Research/code/JPack/projects/co2_hydrogenation")
#
# Rh = Structure.from_file("H_ads/structures/Rh_mp-74_conventional_standard.cif")

with MPRester() as m:
    Rh = m.get_structure_by_material_id("mp-47", conventional_unit_cell=True)

Rh.add_oxidation_state_by_guess()
slabs = generate_all_slabs(Rh, 1, 10, 10)
fig = plt.figure()
for n, slab in enumerate(slabs):
    if slab.is_symmetric() and not slab.is_polar():
        print(slab.miller_index, slab.shift)
#         ax = fig.add_subplot(1, len(slabs), n+1)
#         plot_slab(slab, ax)
#         ax.set_title("{}".format(slab.miller_index))
# plt.show()
print(len(slabs))

#%%
from pymatgen.analysis.adsorption import AdsorbateSiteFinder

with MPRester() as m:
    Rh = m.get_structure_by_material_id("mp-74", conventional_unit_cell=True)
print(Rh.lattice.abc)
slab_100 = SlabGenerator(Rh, [1, 0, 0], 10, 15, center_slab=True).get_slab()

asf = AdsorbateSiteFinder(slab_100)
ads_sites = asf.find_adsorption_sites()




#%%
from pymatgen import Molecule
from monty.serialization import dumpfn

H = Molecule(["H", "H"], [[0,0,0], [0.74, 0, 0]])

ads_sts = asf.generate_adsorption_structures(H, find_args={"distance": 1.6}, repeat=[3,3,1])

fig = plt.figure(dpi=500)
for n, ads_st, name in zip(range(len(ads_sts)), ads_sts, ads_sites.keys()):
    print(ads_st.sites)
    # ads_st.to("poscar", "/Users/jeng-yuantsai/Research/code/JPack/projects/co2_hydrogenation/H_ads/H2_{}_100.vasp".format(name))
    ax = fig.add_subplot(1, 4, n+1)
    plot_slab(ads_st, ax, adsorption_sites=False)
    # ax.set_xlim(-10, 10)
    # ax.set_ylim(-10, 10)
plt.show()
#%%
st_json = dict(zip(ads_sites.keys(), ads_sts))
dumpfn(st_json, "ads_sts.json", indent=4)

#%% H ads results

from qubitPack.tool_box import get_db
from pymatgen import Structure
import os
from matplotlib import pyplot as plt

os.chdir("/Users/jeng-yuantsai/Research/code/JPack/projects/co2_hydrogenation/H_ads/structures/relaxed")
d = get_db("CO2_hydrogenation", "H_ads")

fig = plt.figure(figsize=(20,5))

for tk_id, n, m in zip([12,13,14], ["hollow", "top", "bridge"], range(3)):
    e = d.collection.find_one({"task_id":tk_id})
    d.collection.update_one({"task_id":tk_id}, {"$set": {"ads_site":n}})
    print(e["output"]["energy"])
    # Structure.from_dict(e["output"]["structure"]).to("poscar", "{}_tid_{}.vasp".format(n, tk_id))
    # Structure.from_dict(e["input"]["structure"]).to("poscar", "{}_tid_{}_in.vasp".format(n, tk_id))

    locpot = e["calcs_reversed"][0]["output"]["locpot"]["2"]
    ax = fig.add_subplot(1,1,1)
    ax.plot(locpot)
    ax.set_title(n)
plt.show()
#%%
from pymatgen.electronic_structure.cohp import CompleteCohp
os.chdir('/Users/jeng-yuantsai/Research/project/co2_hydrogenation/calculations/H2-Rh_cohp')

hollow, bridge, top = "tkid_18_hollow", "tkid_20_bridge", "tkid_19_top"
hollow_cohp=CompleteCohp.from_file(fmt="LOBSTER",filename=os.path.join(hollow, "COHPCAR.lobster.gz"),
                                   structure_file=os.path.join(hollow, "POSCAR"))
bridge_cohp=CompleteCohp.from_file(fmt="LOBSTER",filename=os.path.join(bridge, "COHPCAR.lobster.gz"),
                                   structure_file=os.path.join(bridge, "POSCAR"))
top_cohp=CompleteCohp.from_file(fmt="LOBSTER",filename=os.path.join(top, "COHPCAR.lobster.gz"),
                                structure_file=os.path.join(top, "POSCAR"))



#%%
import os
from pymatgen.electronic_structure.plotter import CohpPlotter
from pymatgen.electronic_structure.core import Orbital
from pymatgen.electronic_structure.cohp import CompleteCohp
os.chdir('/Users/jeng-yuantsai/Research/project/co2_hydrogenation/calculations/H-Rh_cohp')

hollow_cohp=CompleteCohp.from_file(fmt="LOBSTER",filename=os.path.join(hollow, "COHPCAR.lobster.gz"),
                                   structure_file=os.path.join(hollow, "POSCAR"))
completecohp = hollow_cohp

#search for the number of the COHP you would like to plot in ICOHPLIST.lobster (the numbers in COHPCAR.lobster are different!)
label="5"
cp=CohpPlotter()
#get a nicer plot label
plotlabel=str(completecohp.bonds[label]['sites'][0].species_string)+'-'+str(completecohp.bonds[label]['sites'][1].species_string)

cp.add_cohp(plotlabel,completecohp.get_orbital_resolved_cohp(label=label, orbitals=[[1, Orbital.s], [4, Orbital.dyz]]))
#check which COHP you are plotting

print("This is a COHP between the following sites: "+str(completecohp.bonds[label]['sites'][0])+' and '+ str(completecohp.bonds[label]['sites'][1]))

x = cp.get_plot(integrated=False)
x.ylim([-10, 6])

x.show()
#%%
from pymatgen.electronic_structure.core import Orbital

# labelist = ["7","8"]
labelist = ["1", "2", "5", "6"]
cp = CohpPlotter()
# get a nicer plot label
plotlabel = "H1-Rh4+Rh5 bonds"

cp.add_cohp(plotlabel, completecohp.get_summed_cohp_by_label_list(
    label_list=labelist, divisor=1))
x = cp.get_plot(integrated=False)
x.ylim([-10, 6])
x.show()

#%%
from pymatgen.io.lobster import Icohplist

#get icohp value by label (labelling according to ICOHPLIST.lobster)
#for spin polarized calculations you can also sum the spin channels
icohplist=Icohplist(filename='tkid_20_bridge/ICOHPLIST.lobster.gz')
icohpcollection=icohplist.icohpcollection

print('icohp value for certain bond by label')
label='7'
print(icohpcollection.get_icohp_by_label(label))
print()
#you can get all Icohpvalue objects for certain bond lengths
print('Icohp values for certain bonds with certain bond lengths')
for key,icohp in icohpcollection.get_icohp_dict_by_bondlengths(minbondlength=0.0, maxbondlength=3.0).items():
    print(key+':'+str(icohp.icohp))
print()
#you can get all icohps for a certain site
print('ICOHP values of certain site')
for key,icohp in (icohpcollection.get_icohp_dict_of_site(site=0,minbondlength=0.0, maxbondlength=3.0).items()):
    print(key+':'+str(icohp.icohp))