from qubitPack.qc_searching.analysis.dos_plot_from_db import DosPlotDB
from qubitPack.qc_searching.analysis.read_eigen import DetermineDefectState
import matplotlib.pyplot as plt
from qubitPack.qc_searching.analysis.dos_plot_from_db import DB_CONFIG_PATH
from glob import glob
import os
import numpy as np


def main(db, db_filter, cbm, vbm, path_save_fig, plot=True, clipboard="tot", locpot=False, db_host=None):
    """
    When one is using "db_cori_tasks_local", one must set ssh-tunnel as following:
    "ssh -f tsaie79@cori.nersc.gov -L 2222:mongodb07.nersc.gov:27017 -N mongo -u 2DmaterialQuantumComputing_admin -p
    tsaie79 localhost:2222/2DmaterialQuantumComputing"
    """

    can = DetermineDefectState(db=db, db_filter=db_filter, cbm=cbm, vbm=vbm, show_edges="band_edges",
                               save_fig_path=path_save_fig, locpot=locpot, db_host=db_host)

    tot, proj = can.get_candidates(
        0,
        threshold=0,
        select_up=None,
        select_dn=None
    )
    if clipboard == "tot":
        tot.to_clipboard("\t")
    else:
        proj.to_clipboard("\t")

    if plot:
        cbm_set, vbm_set = None, None
        if locpot and db_host:
            cbm_set = can.cbm + can.vacuum_locpot
            vbm_set = can.vbm + can.vacuum_locpot
        else:
            cbm_set = cbm
            vbm_set = vbm
        dos_plot = DosPlotDB(db=db, db_filter=db_filter, cbm=cbm_set, vbm=vbm_set, path_save_fig=path_save_fig)
        dos_plot.sites_plots(energy_upper_bound=2, energy_lower_bound=2)
        dos_plot.total_dos(energy_upper_bound=2, energy_lower_bound=2)
        print(dos_plot.nn)
        dos_plot.orbital_plot(dos_plot.nn[-1], 2, 2)
        plt.show()


if __name__ == '__main__':
    proj_path = '/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/search_triplet_from_defect_db'
    save_path = os.path.join(proj_path)

    host_path = "/Users/jeng-yuantsai/Research/project/symBaseBinaryQubit/calculations/scan_relax_pc"
    db_host_json = os.path.join(host_path, "db.json")

    for dir_name in ["defect_states", "structures", "xlsx"]:
        os.makedirs(os.path.join(save_path, dir_name), exist_ok=True)
    # path = '/Users/jeng-yuantsai/Research/qubit/My_manuscript/mx2_antisite_basic/defect_states'
    # db_json = "/Users/jeng-yuantsai/Research/qubit/My_manuscript/WSe_2/support/c2db_TMDC_search"
    db_json = os.path.join(proj_path, "db.json")
    # db_json = '/Users/jeng-yuantsai/Research/qubit/My_manuscript/mx2_antisite_basic/db.json'
    # p1 = os.path.join(db_json, "db_WSe2_like_Ef_from_C2DB.json")
    # p2 = os.path.join(db_json, "db_c2db_tmdc_bglg1.json")

    main(
        db_json,
        {"task_id": 934},
        0, -2.2,
        save_path,
        False,
        "tot",
        True,
        db_host_json
    )
