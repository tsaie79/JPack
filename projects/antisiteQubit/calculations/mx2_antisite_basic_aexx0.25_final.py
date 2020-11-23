#%% Full Table
from atomate.vasp.database import VaspCalcDb
from pymatgen import Structure
import pandas as pd

# path ="/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_final/
path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_latest/"
col = VaspCalcDb.from_db_file(path+"db.json").collection
perturb = {}
unperturb = {}
for i in col.find({
    "perturbed":{"$in":[0, 0.02, 1e-5]},
    "task_label":"HSE_scf",
    "nupdown_set":2
    #"task_id":{"$in":[3187, 3188, 3190, 3189, 3192, 3191]}
}):
    print(i["task_id"])
    if i["perturbed"] == 1e-5:
        st = Structure.from_dict(i["calcs_reversed"][0]["output"]["structure"])
        st.to("poscar", path+"structures/pert_{}_{}.vasp".format(i["formula_pretty"], i["task_id"]))
        if "Te" in i["chemsys"]:
            perturb[i["chemsys"]] = {
                "M-M1":st.get_distance(55, 74),
                "M-M2":st.get_distance(49, 74),
                "M-M3":st.get_distance(54, 74),
                "M-X":st.get_distance(24, 74),
                "X-A-M":st.get_angle(24, 74, 55),
                "energy":i["output"]["energy"]
            }
        else:
            perturb[i["chemsys"]] = {
                "M-M1":st.get_distance(0, 25),
                "M-M2":st.get_distance(6, 25),
                "M-M3":st.get_distance(5, 25),
                "M-X":st.get_distance(50, 25),
                "X-A-M":st.get_angle(50,25,5),
                "energy":i["output"]["energy"]
            }
    elif i["perturbed"] == 0:
        st = Structure.from_dict(i["calcs_reversed"][0]["output"]["structure"])
        st.to("poscar", path+"structures/unpert_{}_{}.vasp".format(i["formula_pretty"], i["task_id"]))
        if "Te" in i["chemsys"]:
            unperturb[i["chemsys"]] = {
                "M-M1":st.get_distance(55, 74),
                "M-M2":st.get_distance(49, 74),
                "M-M3":st.get_distance(54, 74),
                "M-X":st.get_distance(24, 74),
                "X-A-M":st.get_angle(24, 74, 55),
                "energy":i["output"]["energy"],
                "task_id": i["task_id"]
            }
        else:
            unperturb[i["chemsys"]] = {
                "M-M1":st.get_distance(0, 25),
                "M-M2":st.get_distance(6, 25),
                "M-M3":st.get_distance(5, 25),
                "M-X":st.get_distance(50, 25),
                "X-A-M":st.get_angle(50,25,5),
                "energy":i["output"]["energy"],
                "task_id": i["task_id"]
            }
perturb = pd.DataFrame(perturb)
unperturb = pd.DataFrame(unperturb)


perturb = perturb.reindex(columns=["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"]).round(3)
unperturb = unperturb.reindex(columns=["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"]).round(3)

mo_angles_1 = perturb.iloc[4, 3:]
w_angles_1 = perturb.iloc[4, :3]
ax.plot(["S", "Se", "Te"], mo_angles_1, label="Mo_relax", marker="o")
ax.plot(["S", "Se", "Te"], w_angles_1, label="W_relax", marker="x")
ax.legend()
plt.show()

# unperturb.to_clipboard()
# print(perturb)
# print(unperturb)

#%% Defect state plot
from matplotlib.ticker import AutoMinorLocator

def triplet_defect_state_plot(db, filter_label, up, down, save_fig):
    colors = {
        "cbm": "orange",
        "vbm": "deepskyblue"
    }
    evac = []
    for chemsys in filter_label:
        e = db.collection.find_one({"task_label":"HSE_scf", "chemsys":chemsys})
        print(e["task_id"])
        vac = max(e["calcs_reversed"][0]["output"]["locpot"]["2"])
        evac.append(vac)
    evac = np.array(evac)

    vbms_025 = []
    cbms_025 = []
    for chemsys in filter_label:
        filters = {
            "chemsys": chemsys,
            "task_label":"hse gap",
            "input.incar.AEXX": 0.25
        }
        e = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/"
                                    "mx2_antisite_basic_bandgap/db.json").collection.find_one(filters)
        vac = max(e["calcs_reversed"][0]["output"]["locpot"]["2"])
        cbms_025.append(e["output"]["cbm"] - vac)
        vbms_025.append(e["output"]["vbm"] - vac)

    # vbms_025.append(-0.921-evac[0]) for sandwich
    # cbms_025.append(1.145-evac[0]) for sandwich
    cbms_025 = np.array(cbms_025)
    vbms_025 = np.array(vbms_025)

    bandgap = pd.DataFrame({"bandgap":dict(zip(filter_labels, cbms_025-vbms_025))})
    print(bandgap)
    x = np.arange(len(vbms_025))
    emin_025 = min(vbms_025)
    emin = emin_025-0.5


    fig = plt.figure(figsize=(4,6))
    ax = fig.add_subplot(1, 1, 1)

    #up
    up_band = []
    up_band_occ = []
    for chemsys, vac in zip(filter_label, evac):
        e = db.collection.find_one({"task_label":"HSE_scf", "chemsys":chemsys})
        up_ev = []
        up_ev_occ = []
        ev = db.get_eigenvals(e["task_id"])
        for band in up[e['chemsys']]:
            up_ev.append(ev["1"][0][band][0]-vac)
            if ev["1"][0][band][1] == 1:
                up_ev_occ.append(ev["1"][0][band][0]-vac)
        up_band.append(up_ev)
        up_band_occ.append(up_ev_occ)

    for compound_idx in range(len(up_band)):
        ax.hlines(up_band[compound_idx], x[compound_idx]-0.4, x[compound_idx]-0.05)

    # for occ_idx in range(len(up_band_occ)):
    #     for energy, pos in zip(up_band_occ[occ_idx], [0.35, 0.2]):
    #         ax.text(x[occ_idx]-pos, energy, "\u2b06")

    #dn
    dn_band = []
    dn_band_occ = []
    for chemsys, vac in zip(filter_label, evac):
        e = db.collection.find_one({"task_label":"HSE_scf", "chemsys":chemsys})
        ev = db.get_eigenvals(e["task_id"])
        dn_ev = []
        dn_ev_occ = []
        for band in down[e["chemsys"]]:
            print(band)
            if band:
                dn_ev.append(ev["-1"][0][band][0]-vac)
                if ev["-1"][0][band][1] == 1:
                    dn_ev_occ.append(ev["-1"][0][band][0]-vac)
        dn_band.append(dn_ev)
        dn_band_occ.append(dn_ev_occ)
        print(dn_band)
    print(dn_band)
    print(dn_band_occ)

    for compound_idx in range(len(dn_band)):
        ax.hlines(dn_band[compound_idx], x[compound_idx]+0.05, x[compound_idx]+0.4)

    for occ_idx in range(len(dn_band_occ)):
        for energy, pos in zip(dn_band_occ[occ_idx], [-0.1, 0]):
            ax.text(x[occ_idx]-pos, energy, "\u2b07")





    ax.bar(x, np.array(vbms_025) - emin, bottom=emin, color=colors["vbm"])
    ax.bar(x, -np.array(cbms_025), bottom=cbms_025, color=colors["cbm"])

    ax.set_ylim(emin, -3)
    # ax.set_xticks(x)
    # ax.set_xticklabels(labels_in_plot, rotation=0, fontsize=10)
    ax.set_xticks([])
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))


    # plt.ylabel('Energy relative to vacuum (eV)', fontsize=20)
    # plt.tight_layout()
    if save_fig:
        plt.savefig(os.path.join("/Users/jeng-yuantsai/Research/qubit/plt", "{}.eps".format(save_fig)),
                    img_format="eps")
    plt.show()

#%% Defect states
filter_labels = ["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"]
# labels_in_plot = ["${W_{S}}^0$", "${W_{Se}}^0$", "${Mo_{S}}^0$", "${Mo_{Se}}^0$"]
triplet_defect_state_plot(
    VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_basic_aexx0.25_final/db.json"),
    filter_labels,
    {"S-W":[224,225,226],"Se-W":[224,225,226],"Te-W":[224,225,226],
     "Mo-S":[302,303,304],"Mo-Se":[302,303,304], "Mo-Te":[302,303,304]},
    {"S-W":[None],"Se-W":[None],"Te-W":[None], "Mo-S":[None],"Mo-Se":[None], "Mo-Te":[None]},
    "mx2_defect_state_withTe"
)

#%% Trnasition levels figure
from matplotlib.ticker import AutoMinorLocator
from atomate.vasp.database import VaspCalcDb
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os

def transistion_levels(filter_label, tranl):
    colors = {
        "cbm": "orange",
        "vbm": "deepskyblue"
    }
    vbms_025 = []
    cbms_025 = []
    for chemsys in filter_label:
        filters = {
            "chemsys": chemsys,
            "task_label":"hse gap",
            "input.incar.AEXX": 0.25
        }
        e = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/"
                                    "mx2_antisite_basic_bandgap/db.json").collection.find_one(filters)
        vac = max(e["calcs_reversed"][0]["output"]["locpot"]["2"])
        cbms_025.append(e["output"]["cbm"] - vac)
        vbms_025.append(e["output"]["vbm"] - vac)
    cbms_025 = np.array(cbms_025)
    vbms_025 = np.array(vbms_025)

    bandgap = pd.DataFrame({"bandgap":dict(zip(filter_labels, cbms_025-vbms_025))})
    print(bandgap)
    x = np.arange(len(vbms_025))
    emin_025 = min(vbms_025)
    emin = emin_025-0.5

    # ppi = 100
    # figw = 800
    # figh = 800
    # fig = plt.figure(figsize=(figw / ppi, figh / ppi), dpi=ppi)
    fig = plt.figure(figsize=(3,3))
    ax = fig.add_subplot(1, 1, 1)

    for trls, pos, vbm in zip(tranl, x, np.array(vbms_025)):
        text = ["$\epsilon(+/0)$", "$\epsilon(0/-)$"]
        if pos == x[-1]:# for MoTe2
            ax.text(pos-0.4, trls[0]+vbm-0.17, text[0], size=6.5)
            ax.text(pos-0.4, trls[1]+vbm+0.1, text[1], size=6.5)
            ax.hlines(trls[0]+vbm, pos-0.4, pos+0.4)
            ax.hlines(trls[1]+vbm, pos-0.4, pos+0.4)
        elif pos == x[2]:
            ax.text(pos-0.4, trls[0]+vbm-0.15, text[0], size=6.5)
            ax.text(pos-0.4, trls[1]+vbm-0.15, text[1], size=6.5)
            ax.hlines(trls[0]+vbm, pos-0.4, pos+0.4)
            ax.hlines(trls[1]+vbm, pos-0.4, pos+0.4)
        else:
            ax.text(pos-0.4, trls[0]+vbm+0.1, text[0], size=6.5)
            ax.text(pos-0.4, trls[1]+vbm+0.1, text[1], size=6.5)
            ax.hlines(trls[0]+vbm, pos-0.4, pos+0.4)
            ax.hlines(trls[1]+vbm, pos-0.4, pos+0.4)
            # ax.text(pos-0.35, trl+vbm+0.1, text + " {:.2f}".format(trl+vbm))

    # plot bars
    # ax.bar(x, np.array(vbms_c2db) - emin, bottom=emin, label="vbm_c2db")
    ax.bar(x, np.array(vbms_025) - emin, bottom=emin, label="vbm", color=colors["vbm"])
    ax.bar(x, -np.array(cbms_025), bottom=cbms_025, label="cbm", color=colors["cbm"])
    # ax.bar(x, -np.array(cbms_c2db), bottom=cbms_c2db, label="cbm_c2db")


    ax.set_ylim(-6.8, -3)
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    # ax.set_xticks(x)
    # ax.set_xticklabels(labels_in_plot, rotation=0, fontsize=10)
    ax.set_xticks([])


    #         ax.legend(loc="upper right")

    #         plt.title("Defect transition levels (eV)", fontsize=12)
    # plt.ylabel("$\epsilon(q/q')$ Transition levels (eV)", fontsize=20)
    # plt.tight_layout()
    plt.savefig(os.path.join("/Users/jeng-yuantsai/Research/qubit/plt", "mx2_tls_0.25.tempTe.eps"),
                img_format="eps")
    plt.show()



filter_labels = ["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"]
# labels_in_plot = ["${W_{S}}^0$", "${W_{Se}}^0$", "${Mo_{S}}^0$", "${Mo_{Se}}^0$"]
trls = [[2.4945-2.1, 1.82], [2.1977-1.43, 1.65], [1.6603-0.71, 1.43], [2.313-2.01, 1.31], [2.0595-1.93, 1.37], [1.6577-0.88, 0.97]]
transistion_levels(filter_labels, trls)



