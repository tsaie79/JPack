#%%
from atomate.vasp.database import VaspCalcDb
import pandas as pd
from pymatgen import Structure, Element, Spin
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter
from pymatgen.io.vasp.outputs import Vasprun, Locpot, VolumetricData
import os
import numpy as np
from math import floor
import matplotlib.pyplot as plt
from IPython.display import display

db_wse2_like = VaspCalcDb.from_db_file("db_WSe2_like_Ef_from_C2DB.json")
db_wse2_Ef = VaspCalcDb.from_db_file("db_wse2_Ef.json")
db_dk = VaspCalcDb.from_db_file("db_dk_local.json")
db_c2db_tmdc_bglg1 = VaspCalcDb.from_db_file("db_c2db_tmdc_bglg1.json")
db_mx2_antisite_basic_aexx035 = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/qubit/'
                                                        'calculations/mx2_antisite_basic_aexx0.35/db.json')
db_mx2_antisite_basic = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/qubit/'
                                                        'calculations/mx2_antisite_basic/db.json')




#%%
filter_labels = ["S-W"]
labels_in_plot = ["${W_{S}}^0$"]
triplet_defect_state_plot(
    VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_basic_aexx0.25_final/db.json"),
    filter_labels,
    labels_in_plot,
    {"S-W":[224,225,226]},
    {"S-W":[None]}
)
#%% BN-WS2-BN defect states
filter_labels = ["B-N-S-W"]
labels_in_plot = ["${W_{S}}^0$", "${W_{Se}}^0$", "${W_{Te}}^0$", "${Mo_{S}}^0$", "${Mo_{Se}}^0$", "${Mo_{Te}}^0$"]
triplet_defect_state_plot(
    VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/sandwich_BN_mx2/db.json"),
    filter_labels,
    labels_in_plot,
    {"B-N-S-W":[367,368,369]},
    {"B-N-S-W":[None]})



#%%
def ef_025():
    x1 = [0,0.689]
    y1 = [-4.925,-4.236]
    plt.plot(x1, y1, label = "charge state +1", color="purple")
    x2 = [0.689,1.119]
    y2 = [-4.236,-4.236]
    plt.plot(x2, y2, label = "charge state 0", color="purple")
    x3 = [1.119,2]
    y3 = [-4.236,-5.998]
    plt.plot(x3, y3, label="charge state -2", color="purple")
    plt.yticks([-6, -5, -4])
    plt.xticks([0, 0.69, 1.12, 2])
    plt.show()
    plt.savefig("/Users/jeng-yuantsai/Research/qubit/plt/wse2_Ef_0.25.eps", format="eps")


#%%
from glob import glob
import os
from pymatgen import Structure
import pandas as pd
def get_st_data(path):
    data = {}
    data1 = {}
    for s in glob(os.path.join(path, "*HSE*")):
        st = Structure.from_file(s)
        print(s)
        if "Te" in s:
            data1[s.split("/")[-1]] = {"M-M1":st.get_distance(74, 55), "M-M2":st.get_distance(74, 49),
                                      "M-M3":st.get_distance(74, 54), "M-X":st.get_distance(74, 24)}
        else:
            data1[s.split("/")[-1]] = {"M-M1":st.get_distance(0, 25), "M-M2":st.get_distance(6, 25),
                                      "M-M3":st.get_distance(5, 25), "M-X":st.get_distance(50, 25)}
            data[s.split("/")[-1]] = {"a": st.lattice.a, "b": st.lattice.b, "c": st.lattice.c}
            data = pd.DataFrame(data).round(3)
            data1 = pd.DataFrame(data1).round(3)
    return data, data1

lattice, bond = get_st_data('/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_basic/structures')
lattice35, bond35 = get_st_data('/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_basic_aexx0.35/structures')
print(lattice, bond)

#%%
from pymatgen import Structure
col = db_wse2_like.collection
data75 = {}
data108 = {}
for i in col.find({"task_id": {"$in":[1947, 1954, 1949, 1953, 1958, 1948, 1957, 1962, 1968]+
                                    [1877, 1881, 1900, 1889, 1886, 1909, 1876, 1895, 1903]+
                                     [2550, 2545, 2549]+[2547, 2551, 2562]}, "charge_state": -2}):
    st = Structure.from_dict(i["calcs_reversed"][0]["output"]["structure"])
    # st.to("poscar", "/Users/jeng-yuantsai/Research/qubit/calculations/WSe2_like_Ef_from_C2DB/structures/{}.{}.vasp".
    #       format(i["formula_pretty"], i["task_id"]))
    if i["nsites"] == 108:
        data108[i["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"]] = {"M-M1":st.get_distance(7, 36), "M-M2":st.get_distance(0, 36),
                         "M-M3":st.get_distance(6, 36), "M-X":st.get_distance(72, 36)}
    elif i["nsites"] == 75:
        data75[i["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"]] = {"M-M1":st.get_distance(0, 25), "M-M2":st.get_distance(6, 25),
                               "M-M3":st.get_distance(5, 25), "M-X":st.get_distance(50, 25)}

a = pd.DataFrame(data108).round(3)
b = pd.DataFrame(data75).round(3)
pd.concat([a,b], axis=0)



#%%75
db = "mx2_antisite_basic_aexx0.25_final"
col = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/projects/qubit/calculations/{}/db.json".
                              format(db)).collection
chemsys = "Te-W"
task_label = "PBE_relax"
encut = [320]
if "Te" in chemsys:
    for i in col.find({"task_label": task_label, "input.incar.ENCUT":{"$in":encut}, "nsites":5*5*3, "chemsys":chemsys}):
        print("=="*5, i["perturbed"], round(i["output"]["energy"],3))
        init = Structure.from_dict(i["orig_inputs"]["poscar"]["structure"])
        final = Structure.from_dict(i["output"]["structure"])
        for site in [55, 49, 54, 24]:
            print("before:{}, {}, after:{}. {}".format(round(init.get_distance(74,site),3), round(init.lattice.a/5,3),
                                                       round(final.get_distance(74, site), 3), round(final.lattice.a/5,3)))


else:
    for i in col.find({"task_label": task_label, "input.incar.ENCUT":{"$in":encut}, "nsites":5*5*3, "chemsys":chemsys}):
        print("=="*5, i["perturbed"], round(i["output"]["energy"],3))
        init = Structure.from_dict(i["orig_inputs"]["poscar"]["structure"])
        final = Structure.from_dict(i["output"]["structure"])
        for site in [0, 6, 5, 50]:
            print("before:{}, {}, after:{}. {}".format(round(init.get_distance(25,site),3), round(init.lattice.a/5,3),
                                                       round(final.get_distance(25, site), 3), round(final.lattice.a/5,3)))
print(i["calcs_reversed"][0]["dir_name"])
print(i["task_id"])
print("encut:{}, ediffg:{}".format(i["input"]["incar"]["ENCUT"], i["input"]["incar"]["EDIFFG"]))
print(i["calcs_reversed"][0]["input"]["kpoints"]["kpoints"])
print(i["calcs_reversed"][0]["input"]["nkpoints"])



#%% 108
mx2_antisite_perturb_pattern = "mx2_antisite_perturb_pattern"
mx2_antisite_sp_aexx25 = "mx2_antisite_sp_aexx0.25"
col = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/{}/db.json".
                              format(mx2_antisite_perturb_pattern))
col = col.collection

for i in col.find({"defect_pos":"side", "nsites":108}):
    print("=="*5, i["perturbed"], round(i["output"]["energy"],3))
    init = Structure.from_dict(i["orig_inputs"]["poscar"]["structure"])
    final = Structure.from_dict(i["output"]["structure"])
    c = 36
    for site in [0,6,7,72]:
        print("before:{}, {}, after:{}. {}".format(round(init.get_distance(c,site),3), round(init.lattice.a/5,3),
                                                   round(final.get_distance(c, site), 3), round(final.lattice.a/5,3)))
    print(i["dir_name"])
    print(i["task_id"])
    print(i["input"]["incar"]["ENCUT"])
    print(i["calcs_reversed"][0]["input"]["kpoints"]["kpoints"])
    print(i["calcs_reversed"][0]["input"]["nkpoints"])

#%% for antisite at center of supercell
col = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_perturb_pattern/db.json")
col = col.collection

for i in col.find({"defect_pos":"center"}):
    print("=="*5, i["perturbed"], round(i["output"]["energy"],3))
    init = Structure.from_dict(i["orig_inputs"]["poscar"]["structure"])
    final = Structure.from_dict(i["output"]["structure"])
    for site in [17,18,12,62]:
        print("before:{}, after:{}".format(round(init.get_distance(25,site),3), round(final.get_distance(25, site), 3)))
    print(i["dir_name"])


#%% make structure
def get_rand_vec(distance): #was 0.001
    # deals with zero vectors.
    vector = np.random.randn(3)
    vnorm = np.linalg.norm(vector)
    return vector / vnorm * distance if vnorm != 0 else get_rand_vec(distance)

col = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_pc/db.json").collection
wse2_HSE06_relax_host = col.find_one({"task_id":3150})
host_st = Structure.from_dict(wse2_HSE06_relax_host["output"]["structure"])
host_st.make_supercell([5,5,1])
host_st.to("poscar", "/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_pc/wse2_75.vasp")
host_st.replace(37, "W")
distort = 0.02
host_st.translate_sites(37, get_rand_vec(distort))
host_st.to("poscar", "/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_pc/w_se_75_{}.vasp".format(distort))

import subprocess
cmd = "scp /Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_pc/w_se_75_{}.vasp" \
      " tug03990@owlsnesttwo.hpc.temple.edu:/gpfs/work/tug03990/mx2_antisite_perturb_pattern/structures".format(distort)
subprocess.call(cmd.split(" "))


#%% symmetrize strucutre
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.analysis.structure_analyzer import RelaxationAnalyzer

mx2_antisite_basic_aexx25_new = "mx2_antisite_basic_aexx0.25_new"
col = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/{}/db.json".
                              format(mx2_antisite_basic_aexx25_new)).collection

init = Structure.from_dict(col.find_one({"task_id":3126})["output"]["structure"])
print(init.get_space_group_info())
sym = SpacegroupAnalyzer(init, 0.9)
final = sym.get_refined_structure()
# final.to("poscar","/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_basic_aexx0.25_new/refined_wse2.vaspl")
for site in [0, 6, 5, 50]:
    print("before:{}, {}, after:{}. {}".format(round(init.get_distance(25,site),3), round(init.lattice.a/5,3),
                                               round(final.get_distance(25, site), 3), round(final.lattice.a/5,3)))
asym = SpacegroupAnalyzer(final, 0.5)
sym = asym


#%% symmetrize structure but no relax
mx2_antisite_symmetrized_st_aexx25 = "mx2_antisite_symmetrized_st_aexx0.25"
mx2_antisite_c3v_mv_z = "mx2_antisite_c3v_mv_z"
sym = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/{}/db.json".
                              format(mx2_antisite_c3v_mv_z)).collection
mx2_antisite_basic_aexx25_new = "mx2_antisite_basic_aexx0.25_new"
asym = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/{}/db.json".
                               format(mx2_antisite_basic_aexx25_new)).collection

chemsys = "Se-W"
if "Te" in chemsys:
    entry = {}
    for s, a, o in zip(sym.find({"chemsys":chemsys, "task_label":"PBE_relax"}),
                    asym.find({"chemsys":chemsys, "task_label":"PBE_relax", "perturbed":0.02}),
                    asym.find({"chemsys":chemsys, "task_label":"PBE_relax", "perturbed":0})):
        print("=="*5, "s: p:{} id:{} | a: p:{} id:{}".format(s["perturbed"], s["task_id"], a["perturbed"], a["task_id"]))
        init = Structure.from_dict(s["output"]["structure"])
        final = Structure.from_dict(a["output"]["structure"])
        orig = Structure.from_dict(o["output"]["structure"])
        entry["symmetrize"] = dict(zip(["X1", "X2", "X3", "Z1"],
                                [init.get_distance(74,site) for site in [55, 49, 54, 24]]))
        entry["symmetrize"].update({"energy": s["output"]["energy"]})
        entry["symmetrize"].update({"task_id": s["task_id"]})

        entry["perturb"] = dict(zip(["X1", "X2", "X3", "Z1"],
                                 [final.get_distance(74,site) for site in [55, 49, 54, 24]]))
        entry["perturb"].update({"energy": a["output"]["energy"]})
        entry["perturb"].update({"task_id": a["task_id"]})

        entry["no_perturb"] = dict(zip(["X1", "X2", "X3", "Z1"],
                                 [orig.get_distance(74,site) for site in [55, 49, 54, 24]]))
        entry["no_perturb"].update({"energy": o["output"]["energy"]})
        entry["no_perturb"].update({"task_id": o["task_id"]})
else:
    entry = {}
    for s, a, o in zip(sym.find({"chemsys":chemsys, "task_label":"PBE_relax"}),
                    asym.find({"chemsys":chemsys, "task_label":"PBE_relax", "perturbed":0.02, "task_id":{"$nin":[3156]}}),
                    asym.find({"chemsys":chemsys, "task_label":"PBE_relax", "perturbed":0, "task_id":{"$nin":[3218]}})):
        print("=="*5, "s: p:{} id:{} | a: p:{} id:{}".format(s["perturbed"], s["task_id"], a["perturbed"], a["task_id"]))
        init = Structure.from_dict(s["output"]["structure"])
        final = Structure.from_dict(a["output"]["structure"])
        orig = Structure.from_dict(o["output"]["structure"])
        entry["symmetrize"] = dict(zip(["X1", "X2", "X3", "Z1"],
                                [init.get_distance(25,site) for site in [0, 6, 5, 50]]))
        entry["symmetrize"].update({"energy": s["output"]["energy"]})
        entry["symmetrize"].update({"task_id": s["task_id"]})

        entry["perturb"] = dict(zip(["X1", "X2", "X3", "Z1"],
                                 [final.get_distance(25,site) for site in [0, 6, 5, 50]]))
        entry["perturb"].update({"energy": a["output"]["energy"]})
        entry["perturb"].update({"task_id": a["task_id"]})
        entry["no_perturb"] = dict(zip(["X1", "X2", "X3", "Z1"],
                                 [orig.get_distance(25,site) for site in [0, 6, 5, 50]]))
        entry["no_perturb"].update({"energy": o["output"]["energy"]})
        entry["no_perturb"].update({"task_id": o["task_id"]})


pd.DataFrame(entry).round(3)






#%%
import matplotlib.pyplot as plt
x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]
fig = plt.figure()
ax = fig.add_subplot(111)
ax1 = fig.add_subplot(212)
ax.bar(1, 2, 0.4, 1)
ax1.bar(2, 4, 0.4, 1)
ax.text(1, 3,"x")

#%%
from atomate.vasp.database import VaspCalcDb
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.inputs import Element

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/qubit/calculations/sandwich_BN_mx2/db.json')
dos = db.get_dos(3951)
dos_plt = DosPlotter(stack=False, sigma=0.05)
# tdos = DosPlotter(stack=True, zero_at_efermi=False)
# tdos.add_dos("tdos", dos)
# fig2 = tdos.get_plot(xlim=[-5,5])
# fig2.show()
dos_plt.add_dos("Total DOS", dos)
for element, dos in dos.get_element_dos().items():
    if element == Element("B"):
        dos_plt.add_dos("{}".format(element), dos)
    elif element == Element("N"):
        dos_plt.add_dos("{}".format(element), dos)
    # elif element == Element("W"):
    #     dos_plt.add_dos("{}".format(element), dos)
    # elif element == Element("S"):
fig = dos_plt.get_plot(xlim=[-5,5])

fig.savefig("/Users/jeng-yuantsai/Research/qubit/plt/sandwich_ws2_bn.eps", format="eps")
fig.show()

#%%
from pymatgen import Structure
from atomate.vasp.database import VaspCalcDb
import pandas as pd
import numpy as np

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/qubit/calculations/antisiteQubit/MxC3vToChDeltaE/db.json')
col = db.collection

data = [
    {
        "task_id": e["task_id"],
        "chemsys": e["chemsys"],
        "perturbed":e["perturbed"],
        "energy":e["output"]["energy"],
        "dx":np.array([round(Structure.from_dict(e["input"]["structure"]).get_distance(e["NN"][-1], i),3)
                       for i in e["NN"][:-1]+[24 if "Te" in e['chemsys'] else 50]])
    }
    for e in col.find({"task_label":"HSE_scf"})]
raw = pd.DataFrame(data).sort_values(["chemsys", "perturbed"])
print(raw)

d = []
for chem in ["Mo-S", "S-W", "Mo-Se", "Se-W", "Mo-Te", "Te-W"]:
    a = raw.loc[raw["chemsys"]==chem]["energy"].diff().iloc[-1]
    b = raw.loc[raw["chemsys"]==chem]["dx"].diff().iloc[-1]
    e = {"chemsys": chem, "dE":round(a,3), "dx":np.linalg.norm(b).round(3)}
    d.append(e)
pd.DataFrame(d).set_index("chemsys")

#%%
from pymatgen import Structure
from atomate.vasp.database import VaspCalcDb
import pandas as pd
import numpy as np

db = VaspCalcDb.from_db_file('/Users/jeng-yuantsai/Research/project/qubit/calculations/antisiteQubit/scan_opt_test/db.json')
col = db.collection

data = [
    {
        "formula": e["formula_pretty"],
        "origin_a": e["input"]["structure"]["lattice"]["a"],
        "a":e["output"]["structure"]["lattice"]["a"],
        "da": (e["input"]["structure"]["lattice"]["a"] - e["output"]["structure"]["lattice"]["a"])/
              e["input"]["structure"]["lattice"]["a"],
        "origin_b": e["input"]["structure"]["lattice"]["b"],
        "b":e["output"]["structure"]["lattice"]["b"],
        "db": (e["input"]["structure"]["lattice"]["b"] - e["output"]["structure"]["lattice"]["b"])/
              e["input"]["structure"]["lattice"]["b"],
        "origin_c": e["input"]["structure"]["lattice"]["c"],
        "c":e["output"]["structure"]["lattice"]["c"],
        "dc": (e["input"]["structure"]["lattice"]["c"] - e["output"]["structure"]["lattice"]["c"])/
              e["input"]["structure"]["lattice"]["c"]
    }
    for e in col.find()]
raw = pd.DataFrame(data).round(3).sort_values(["formula"])
print(raw)

