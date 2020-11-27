from qubitPack.qc_searching.analysis.dos_plot_from_db import DosPlotDB
from qubitPack.qc_searching.analysis.read_eigen import DetermineDefectState
import matplotlib.pyplot as plt
from qubitPack.qc_searching.analysis.dos_plot_from_db import DB_CONFIG_PATH
from glob import glob
import os
import pandas as pd
import numpy as np
from collections import Counter


def get_defect_state(db, db_filter, cbm, vbm, path_save_fig, plot=True, clipboard="tot", locpot=None, threshold=0.1):
    """
    When one is using "db_cori_tasks_local", one must set ssh-tunnel as following:
    "ssh -f tsaie79@cori.nersc.gov -L 2222:mongodb07.nersc.gov:27017 -N mongo -u 2DmaterialQuantumComputing_admin -p
    tsaie79 localhost:2222/2DmaterialQuantumComputing"
    """

    can = DetermineDefectState(db=db, db_filter=db_filter, cbm=cbm, vbm=vbm, save_fig_path=path_save_fig, locpot=locpot)
    # can.nn = [25, 26, 31, 30, 29, 49, 45]
    tot, proj, d_df = can.get_candidates(
        0,
        threshold=threshold,
        select_up=None,
        select_dn=None
    )
    e = {}
    e.update(
        {
            "up_band_idx": d_df["up_band_idx"][0],
            "up_from_vbm": d_df["up_from_vbm"][0],
            "up_occ": d_df["up_occ"][0],
            "dn_band_idx": d_df["dn_band_idx"][0],
            "dn_from_vbm": d_df["dn_from_vbm"][0],
            "dn_occ": d_df["dn_occ"][0]
        }
    )
    # well-defined in-gap state: energetic difference of occupied states and vbm > 0.1
    sum_occ_up = 0
    for en_up, occ_up in zip(d_df["up_from_vbm"][0], d_df["up_occ"][0]):
        if abs(en_up) > 0.1 and occ_up > 0.2:
            sum_occ_up += occ_up
    sum_occ_dn = 0
    for en_dn, occ_dn in zip(d_df["dn_from_vbm"][0], d_df["dn_occ"][0]):
        if abs(en_dn) > 0.1 and occ_dn > 0.2:
            sum_occ_dn += occ_dn

    e.update({"up_deg": list(Counter(np.around(e["up_from_vbm"], 2)).values()),
              "dn_deg": list(Counter(np.around(e["dn_from_vbm"], 2)).values())})
    if sum_occ_up > sum_occ_dn:
        e.update({"triplet_from": "up"})
    else:
        e.update({"triplet_from": "dn"})

    # Calculate plausible optical transition energy
    for en_up, occ_up, idx in zip(d_df["up_from_vbm"][0], d_df["up_occ"][0], range(len(d_df["up_occ"][0]))):
        if occ_up == 1 and idx != 0:
            h_occ_en = en_up
            l_unocc_en = d_df["up_from_vbm"][0][idx-1]
            e.update({"up_tran_en": round(l_unocc_en - h_occ_en, 3)})
            break
        elif occ_up > 0 and idx != 0 and d_df["up_occ"][0][0] == 0:
            h_occ_en = en_up
            l_unocc_en = d_df["up_from_vbm"][0][idx-1]
            e.update({"up_tran_en": round(l_unocc_en - h_occ_en, 3)})
            break
        elif occ_up > 0 and idx != 0 and d_df["up_occ"][0][0] != 0:
            h_occ_en = d_df["up_from_vbm"][0][idx+1]
            l_unocc_en = en_up
            e.update({"up_tran_en": round(l_unocc_en - h_occ_en, 3)})
            break
        else:
            e.update({"up_tran_en": 0})
    for en_dn, occ_dn, idx in zip(d_df["dn_from_vbm"][0], d_df["dn_occ"][0], range(len(d_df["dn_occ"][0]))):

        if occ_dn == 1 and idx != 0:
            h_occ_en = en_dn
            l_unocc_en = d_df["dn_from_vbm"][0][idx-1]
            e.update({"dn_tran_en": round(l_unocc_en - h_occ_en, 3)})
            break
        elif occ_dn > 0 and idx != 0 and d_df["dn_occ"][0][0] == 0:
            h_occ_en = en_dn
            l_unocc_en = d_df["dn_from_vbm"][0][idx-1]
            e.update({"dn_tran_en": round(l_unocc_en - h_occ_en, 3)})
            break
        elif occ_dn > 0 and idx != 0 and d_df["dn_occ"][0][0] != 0:
            h_occ_en = d_df["dn_from_vbm"][0][idx+1]
            l_unocc_en = en_dn
            e.update({"dn_tran_en": round(l_unocc_en - h_occ_en, 3)})
            break
        else:
            e.update({"dn_tran_en": 0})

    d_df = pd.DataFrame([e]).transpose()
    print(d_df)
    print("=="*20)
    print(proj)
    print("=="*20)
    print(tot)
    # print("=="*20)
    # print(d_df)


    if clipboard == "tot":
        tot.to_clipboard("\t")
    elif clipboard == "proj":
        proj.to_clipboard("\t")
    elif clipboard == "dist":
        d_df.to_clipboard("\t")

    if plot:
        cbm_set, vbm_set = None, None
        if locpot:
            cbm_set = can.cbm + can.vacuum_locpot
            vbm_set = can.vbm + can.vacuum_locpot
        else:
            cbm_set = cbm
            vbm_set = vbm
        dos_plot = DosPlotDB(db=db, db_filter=db_filter, cbm=cbm_set, vbm=vbm_set, path_save_fig=path_save_fig)
        # dos_plot.nn = [25, 26, 31, 30, 29, 49, 45]
        dos_plot.total_dos(energy_upper_bound=2, energy_lower_bound=2)
        print(dos_plot.nn)
        dos_plot.sites_plots(energy_upper_bound=2, energy_lower_bound=2)
        dos_plot.orbital_plot(dos_plot.nn[-1], 2, 2)
        # plt.show()

        if path_save_fig:
            for df, df_name in zip([tot, proj, d_df], ["tot", "proj", "d_state"]):
                path = os.path.join(path_save_fig, "xlsx", "{}_{}_{}_{}.xlsx".format(
                    can.entry["formula_pretty"],
                    can.entry["task_id"],
                    can.entry["task_label"],
                    df_name
                ))
                df.to_excel(path)

    return tot, proj, d_df

if __name__ == '__main__':
    # proj_path = '/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/search_triplet_from_defect_db'
    proj_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_basic_aexx0.25_sigma_test"
    db_json = os.path.join(proj_path, "db.json")

    # host_path = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc"
    host_path = "/Users/jeng-yuantsai/Research/project/qubit/calculations/mx2_antisite_pc"
    db_host_json = os.path.join(host_path, "db.json")

    # for dir_name in ["defect_states", "structures", "xlsx"]:
    #     os.makedirs(os.path.join(save_path, dir_name), exist_ok=True)
    # path = '/Users/jeng-yuantsai/Research/qubit/My_manuscript/mx2_antisite_basic/defect_states'
    # db_json = "/Users/jeng-yuantsai/Research/qubit/My_manuscript/WSe_2/support/c2db_TMDC_search"
    # db_json = '/Users/jeng-yuantsai/Research/qubit/My_manuscript/mx2_antisite_basic/db.json'
    # p1 = os.path.join(db_json, "db_WSe2_like_Ef_from_C2DB.json")
    # p2 = os.path.join(db_json, "db_c2db_tmdc_bglg1.json")

    tot, proj, d_df = get_defect_state(
        db_json,
        {"task_id": 4228},
        1.5,-0.581,
        None,
        True,
        "dist",
        None,
        0.1
    )
    tot.to_clipboard()
