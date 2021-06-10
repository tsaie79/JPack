#%% Mv_z
from atomate.vasp.database import VaspCalcDb
from pymatgen.io.vasp.inputs import Structure
import pandas as pd
import numpy as np
import os
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt

class MovingZSplitting:
    @classmethod
    def sheet(cls, chemsys):
        db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_splitting/db.json")
        mx2_antisite_mv_z = db.collection

        orig_db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_basic_aexx0.25_final/db.json")


        # wse2_orig = Structure.from_dict(mx2_antisite_mv_z.find_one({"task_id":2646})["output"]["structure"])
        p = "/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_basic_splitting/"
        # e_orig = orig_db.collection.find_one({"chemsys":chemsys, "task_label":"HSE_scf", "nupdown_set":2})
        # orig = Structure.from_dict(e_orig["output"]["structrue"])
        # dz = [orig[site].c for site in e_orig["NN"][:-1]]
        # if "Te" in chemsys:
        #     d0 = round(orig.get_distance(74, 24),3)
        sheet = []
        for e in list(mx2_antisite_mv_z.find({"chemsys":chemsys, "task_label":"HSE_scf", "lattice_constant":{"$nin":["PBE"]}}))+ \
                 list(orig_db.collection.find({"chemsys":chemsys, "task_label":"HSE_scf", "nupdown_set":2})):
            print(e["task_id"])
            try:
                eig = db.get_eigenvals(e["task_id"])
            except TypeError:
                eig = orig_db.get_eigenvals(e["task_id"])
            entry = {}
            st = Structure.from_dict(e["output"]["structure"])
            avg_z = np.average([st[site].z for site in e["NN"][:-1]])
            entry["perturb"] = e["perturbed"]
            entry["formula_pretty"] = e["formula_pretty"]
            entry["bond_length_diff"] = round(st[e["NN"][-1]].z-avg_z,3)
            entry["energy"] = round(e["output"]["energy"],3)
            if "Mo" in chemsys:
                entry["∆E02"] = round(eig["1"][0][304][0]-eig["1"][0][302][0],3)
                entry["∆E12"] = round((eig["1"][0][304][0]-eig["1"][0][303][0]),3 )#/entry["∆E02"],3)
                entry["∆E01"] = round((eig["1"][0][303][0]-eig["1"][0][302][0]), 3)#/entry["∆E02"],3)
            else:
                entry["∆E02"] = round(eig["1"][0][226][0]-eig["1"][0][224][0],3)
                entry["∆E12"] = round((eig["1"][0][226][0]-eig["1"][0][225][0]), 3)#/entry["∆E02"],3)
                entry["∆E01"] = round((eig["1"][0][225][0]-eig["1"][0][224][0]),3) #/entry["∆E02"],3)
            entry["task_id"] = e["task_id"]
            sheet.append(entry)
        df = pd.DataFrame(sheet).set_index("task_id")
        df.sort_values(["bond_length_diff"], inplace=True)
        fig, ax = plt.subplots(3, sharex=True, figsize=(10,8))
        fig.suptitle(df["formula_pretty"].iloc[0])
        ax[0].plot(df["bond_length_diff"], df["∆E02"], "o", color="k", label="∆E02")
        ax[0].set(ylabel="∆E02 (eV)", xlim=[-0.6, 0.6])
        ax[1].plot(df["bond_length_diff"], df["∆E01"], "o", color="blue", label="∆E01")
        ax[1].plot(df["bond_length_diff"], df["∆E12"], "o", color="red", label="∆E12")
        ax[1].set(ylabel="∆E/∆E02", xlim=[-0.6, 0.6])
        ax[2].plot(df["bond_length_diff"], df["energy"], "o", color="k", label="Etot")
        ax[2].set(xlabel="displacement from initial z-position (Å)", ylabel="total energy (eV)", xlim=[-0.6, 0.6])
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        ax[0].legend()
        ax[1].legend()
        ax[2].legend()
        os.makedirs(os.path.join(p, chemsys, "plt"), exist_ok=True)
        fig.savefig(os.path.join(p, chemsys, "plt", "{}.jpg".format(df["formula_pretty"].iloc[0])), dpi=160, format="jpg")
        return df

    @classmethod
    def plot(cls):
        # colors = ["powderblue", "deepskyblue", "red", "orange"]
        # colors = ["lightgreen", "limegreen", "red"]
        colors = ["darkgray", "dimgray", "white", "None"]

        fig = plt.figure(figsize=(3,3))
        ax = fig.add_subplot(1, 1, 1)

        c = [3.15, 3.28, 3.51]
        # c = [3.15, 3.28]

        w = 0.03
        maxy, miny = 1.49+1, 1.671-0.6


        mos2 = [1.48, 1.58, 1.68, 1.78, 1.88, 1.98, 2.08, 2.18, 2.28, 2.38, 2.48]
        ax.bar(c[0]-w/2*1.2,  miny-1.49, w, 1.49, color=colors[0], label="1-2 configuration")
        ax.bar(c[0]-w/2*1.2, maxy-1.49, w, 1.49, color=colors[1], label="2-1 configuration")
        ax.hlines(1.98, c[0]-w/2*1.2-w/2, c[0]-w/2*1.2+w/2, colors[-2])
        ax.hlines(1.49, c[0]-w/2*1.2-w/2, c[0]-w/2*1.2+w/2, colors[-1])

        # ax.plot([c[0]-w/2*1.2, c[0]+w/2*1.2], [0.8, 1.49-0.4], color="silver", linestyle="dashed")

        ws2 = [1.405, 1.505, 1.605, 1.705, 1.805, 1.905, 2.005, 2.105, 2.205, 2.305, 2.405]
        ax.bar(c[0]+w/2*1.2, miny-1.505, w, 1.505, color=colors[0])
        ax.bar(c[0]+w/2*1.2, maxy-1.505, w, 1.505, color=colors[1])
        ax.hlines(1.905, c[0]+w/2*1.2-w/2, c[0]+w/2*1.2+w/2, colors[-2])
        ax.hlines(1.505, c[0]+w/2*1.2-w/2, c[0]+w/2*1.2+w/2, colors[-1])

        # ax.plot([1.125, x2], [0.8, 1.505-0.5], color="silver", linestyle="dashed")

        # ax.plot([c[0]-w/2*1.2+w/2, c[0]+w/2*1.2-w/2], [1.98, 1.905], color="black", linestyle="dashed")
        # ax.plot([c[0]+w/2*1.2+w/2, c[1]-w/2*1.2-w/2], [1.905, 1.924], color="black", linestyle="dashed")
        # ax.plot([c[0]-w/2*1.2+w/2, c[0]+w/2*1.2-w/2], [1.49, 1.505], color="black", linestyle="dashed")
        # ax.plot([c[0]+w/2*1.2+w/2, c[1]-w/2*1.2-w/2], [1.505, 1.524], color="black", linestyle="dashed")

        mose2 = [1.424, 1.524, 1.624, 1.724, 1.824, 1.924, 2.024, 2.124, 2.224, 2.324, 2.424]
        ax.bar(c[1]-w/2*1.2, miny-1.524, w, 1.524, color=colors[0])
        ax.bar(c[1]-w/2*1.2, maxy-1.524, w, 1.524, color=colors[1])
        ax.hlines(1.924, c[1]-w/2*1.2-w/2, c[1]-w/2*1.2+w/2, colors[-2])
        ax.hlines(1.524, c[1]-w/2*1.2-w/2, c[1]-w/2*1.2+w/2, colors[-1])

        # ax.plot([x1+0.1, x1], [0.8, 1.49-0.4], color="silver", linestyle="dashed")

        wse2 = [1.336, 1.436, 1.536, 1.636, 1.736, 1.836, 1.936, 2.036, 2.136, 2.236, 2.336]
        ax.bar(c[1]+w/2*1.2, miny-1.636, w, 1.636, color=colors[0])
        ax.bar(c[1]+w/2*1.2, maxy-1.636, w, 1.636, color=colors[1])
        ax.hlines(1.836, c[1]+w/2*1.2-w/2, c[1]+w/2*1.2+w/2, colors[-2])
        ax.hlines(1.636, c[1]+w/2*1.2-w/2, c[1]+w/2*1.2+w/2, colors[-1])
        # ax.plot([1.125, x2], [0.8, 1.505-0.5], color="silver", linestyle="dashed")

        # ax.plot([c[1]-w/2*1.2+w/2, c[1]+w/2*1.2-w/2], [1.924, 1.836], color="black", linestyle="dashed")
        # ax.plot([c[1]+w/2*1.2+w/2, c[2]-w/2*1.2-w/2], [1.836, 1.402], color="black", linestyle="dashed")
        # ax.plot([c[1]-w/2*1.2+w/2, c[1]+w/2*1.2-w/2], [1.524, 1.636], color="black", linestyle="dashed")
        # ax.plot([c[1]+w/2*1.2+w/2, c[2]-w/2*1.2-w/2], [1.636, 1.902], color="black", linestyle="dashed")

        mote2 = [1.402, 1.502, 1.602, 1.702, 1.802, 1.902]
        ax.bar(c[2]-w/2*1.2, miny-1.902, w, 1.902, color=colors[0])
        ax.bar(c[2]-w/2*1.2, maxy-1.902, w, 1.902, color=colors[1])
        ax.hlines(1.402, c[2]-w/2*1.2-w/2, c[2]-w/2*1.2+w/2, colors[-2])
        ax.hlines(1.902, c[2]-w/2*1.2-w/2, c[2]-w/2*1.2+w/2, colors[-1])
        # ax.plot([2.625, x1], [0.8, 1.902-0.9], color="silver", linestyle="dashed")


        wte2 = [1.271, 1.371, 1.471, 1.571, 1.671, 1.771]
        ax.bar(c[2]+w/2*1.2, miny-1.671, w, 1.671, color=colors[0])
        ax.bar(c[2]+w/2*1.2, maxy-1.671, w, 1.671, color=colors[1])
        ax.hlines(1.271, c[2]+w/2*1.2-w/2, c[2]+w/2*1.2+w/2, colors[-2])
        ax.hlines(1.671, c[2]+w/2*1.2-w/2, c[2]+w/2*1.2+w/2, colors[-1])
        # ax.plot([2.625, x2], [0.8, 1.671-0.6], color="silver", linestyle="dashed")
        # ax.plot([c[2]-w/2*1.2+w/2, c[2]+w/2*1.2-w/2], [1.402, 1.271], color="black", linestyle="dashed")
        # ax.plot([c[2]-w/2*1.2+w/2, c[2]+w/2*1.2-w/2], [1.902, 1.671], color="black", linestyle="dashed")


        # pos = [1, 1.25, 1.75, 2.0, 2.5, 2.75]
        # ax.plot(pos, [1.98, 1.905, 1.924, 1.836, 1.402, 1.271], color="orangered", marker=".", linestyle="dotted")
        # ax.plot(pos, [1.49, 1.505, 1.524, 1.636, 1.902, 1.671], color="royalblue", marker=".", linestyle="dotted")


        # labels = ["${Mo_{S}}^0$", "${W_{S}}^0$",  "${Mo_{Se}}^0$", "${W_{Se}}^0$", "${Mo_{Te}}^0$", "${W_{Te}}^0$"]
        labels = ["3.15", "3.28", "3.51"]
        x1 = 0.25
        # c = []
        # for i in range(3):
        #     x1 += 0.75
        #     x2 = x1+0.25
        #     c.append(x1)
        #     c.append(x2)
        ax.set_xticks(c)
        # ax.set_xticklabels(labels)
        # ax.set_ylabel("Antisite position along z direction (Å)")
        ax.yaxis.set_minor_locator(AutoMinorLocator(4))
        ax.set_ylim([1, 2.6])
        # ax.legend()

        fig.show()
        fig.savefig('/Users/jeng-yuantsai/Research/project/qubit/plt/in_paper/fig_1_c.eps', format="eps")

MovingZSplitting.plot()

#%% MV_Z aggregation
db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_splitting/db.json")

es = db.collection.aggregate(
    [
        {"$match": {"task_label":"HSE_scf", "lattice_constant":{"$nin":["PBE"]}}},
        {"$group": {"_id":"$perturbed",
                    "tid": {"$push": "$task_id"},
                    "pc": {"$push":"$pc_from"},
                    "NN": {"$push":"$NN"},
                    "point_group": {"$push": "$output.spacegroup.point_group"},
                    "mag": {"$push": "$calcs_reversed.output.outcar.total_magnetization"},
                    }
         },
    ]
)


df = pd.DataFrame(es)
df = df.transpose()
#%% get defect states from moving z

from qubitPack.qc_searching.analysis.main import main
proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_splitting"
db_json = os.path.join(proj_path, "db.json")

host_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_pc"
db_host_json = os.path.join(host_path, "db.json")

perturb = 0.5

chem = "S-W"
tot, proj, d_df = main(
    db_json,
    {"task_label":"HSE_scf", "lattice_constant":{"$nin":["PBE"]}, "chemsys":chem, "perturbed":perturb},
    5,-5,
    os.path.join(proj_path, chem),
    True,
    "dist",
    db_host_json,
    0.05
)
tot.to_clipboard()