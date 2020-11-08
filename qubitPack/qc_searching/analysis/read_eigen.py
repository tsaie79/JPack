from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.electronic_structure.core import Spin
from atomate.vasp.database import VaspCalcDb
import numpy as np
import pandas as pd
import os, json
import matplotlib.pyplot as plt
import plotly
import plotly.tools as tls
from qubitPack.qc_searching.py_energy_diagram.application.defect_levels import EnergyLevel
from collections import defaultdict
from qubitPack.qc_searching.analysis.dos_plot_from_db import DB_CONFIG_LOCAL, FIG_SAVING_FILE_PATH


class EigenParse:
    def __init__(self, path):
        self.path = path
        self.vp = Vasprun(os.path.join(self.path, "vasprun.xml"))

    def eigenvalue(self, spin, k_index):
        global s
        if spin == 1:
            s = Spin.up
        elif spin == 2:
            s = Spin.down

        ev = self.vp.eigenvalues
        print(self.vp.actual_kpoints)
        print("=="*50)
        return ev[s][k_index]

    def usage(self):
        v = []
        for i in range(4):
            for j in [1, 2]:
                v.append(self.eigenvalue(j, i))
        comb = np.concatenate(tuple(v), 1)
        # pd.DataFrame(defect.eigenvalue(1, 0)).to_excel("./pic/defect-1.xlsx", index=False)
        # pd.DataFrame(host.eigenvalue(1, 0)).to_excel("./pic/host-1.xlsx", index=False)

        pd.DataFrame(comb).to_excel(os.path.join(self.path, "ev.xlsx"), index=False)

        return None


class ProcarParse:
    def __init__(self, path):
        self.path = path
        self.procar = Procar(os.path.join(self.path, "PROCAR"))
        self.vasp_run = Vasprun(os.path.join(self.path, "vasprun.xml"))
        self.total_up = None
        self.total_down = None
        self.up_energy = None
        self.down_energy = None

    def band(self, write_out, adjacent, k_index, show_nonadjacent=False):
        band = self.procar
        pro_up = []
        pro_down = []
        total_up = []
        total_down = []
        print(self.path.split("/")[-3])
        for j in range(2, 169+1):
            j = j - 2
            print("---band %d---" % (j+2))
            sum_up = 0
            sum_down = 0
            array_up = []
            array_down = []
            if not show_nonadjacent:
                for i in adjacent:
                    up = band.data[Spin.up][k_index][j][i]
                    sum_up += np.sum(up)
                    up = np.append(up, np.sum(up))
                    down = band.data[Spin.down][k_index][j][i]
                    sum_down += -np.sum(down)
                    down = np.append(down, np.sum(down))
                    print(up, down)
                    array_up.append((i, up.tolist()))
                    array_down.append((i, down.tolist()))
                pro_up.append((j+2, dict(array_up)))
                pro_down.append((j+2, dict(array_down)))
                total_up.append(sum_up)
                total_down.append(sum_down)

            elif show_nonadjacent:
                for i in set(range(band.nions))-set(adjacent):
                    up = band.data[Spin.up][k_index][j][i]
                    up = np.append(up, np.sum(up))
                    down = band.data[Spin.down][k_index][j][i]
                    down = np.append(down, np.sum(down))
                    print(up, down)
                    array_up.append((i, up.tolist()))
                    array_down.append((i, down.tolist()))
                pro_up.append((j+2, dict(array_up)))
                pro_down.append((j+2, dict(array_down)))
            else:
                print("show_nonadjacent = True if one wants to show ions except adjacent ones.")
        print(total_up, total_down)
        if write_out and not show_nonadjacent:
            with open(os.path.join(self.path, self.path.split("/")[-3]+"_up_procar.js"), "w") as f:
                json.dump(dict(pro_up), f, indent=4)
                f.close()
            with open(os.path.join(self.path, self.path.split("/")[-3]+"_down_procar.js"), "w") as f:
                json.dump(dict(pro_down), f, indent=4)
                f.close()

        elif write_out and show_nonadjacent:
            with open(os.path.join(self.path, self.path.split("/")[-3]+"_nonadjacent_up_procar.js"), "w") as f:
                json.dump(dict(pro_up), f, indent=4)
                f.close()
            with open(os.path.join(self.path, self.path.split("/")[-3]+"_nonadjacent_down_procar.js"), "w") as f:
                json.dump(dict(pro_down), f, indent=4)
                f.close()
        else:
            print("Write-out was not processed!")

        def print_total():
            y1 = total_up
            y2 = total_down
            x = range(2, 169+1)
            fig, ax = plt.subplots()
            ax.plot(x, y1, "-", label="up")
            ax.plot(x, y2, "-", label="down")
            ax.set_xlabel('Band index')
            ax.set_ylabel('Total projection')
            ax.set_title(self.path.split("/")[-3] + '-Net Orbital Projection on Adjacent Ions' + " (k index=%d)" % k_index)
            plt.grid(linestyle="--")
            # plt.axvline(x=80, color="k", linestyle="-")

            plotly_plot = tls.mpl_to_plotly(fig)
            plotly_plot['data'][0]['hovertext'] = \
                ["%.3f" % self.vasp_run.eigenvalues[Spin.up][k_index][i-2][0] for i in x]
            plotly_plot['data'][1]['hovertext'] = \
                ["%.3f" % self.vasp_run.eigenvalues[Spin.down][k_index][i-2][0] for i in x]
            plotly.offline.plot(plotly_plot, filename=os.path.join(self.path, "procar.html"))

        print_total()
        return total_up, total_down


class DetermineDefectState:
    def __init__(self, show_edges, db, db_filter, cbm, vbm, save_fig_path, locpot=None):
        db = VaspCalcDb.from_db_file(db)

        self.save_fig_path = save_fig_path

        self.entry = db.collection.find_one(db_filter)

        print("---task_id: %s---" % self.entry["task_id"])

        if locpot:
            self.vacuum_locpot = max(self.entry["calcs_reversed"][0]["output"]["locpot"]["2"])

            db_host = VaspCalcDb.from_db_file(locpot)
            self.entry_host = db_host.collection.find_one({"task_id": int(self.entry["pc_from"].split("/")[-1])})
            self.vacuum_locpot_host = max(self.entry_host["calcs_reversed"][2]["output"]["locpot"]["2"])

            self.cbm = self.entry_host["output"]["cbm"] - self.vacuum_locpot_host #A
            self.vbm = self.entry_host["output"]["vbm"] - self.vacuum_locpot_host
            efermi_to_defect_vac = self.entry["calcs_reversed"][0]["output"]["efermi"] - self.vacuum_locpot

            self.efermi = efermi_to_defect_vac

        else:
            print("No vacuum alignment!")
            self.vacuum_locpot = 0
            self.efermi = self.entry["calcs_reversed"][0]["output"]["efermi"]
            self.cbm = cbm
            self.vbm = vbm

        # if locpot and self.entry["vacuum_locpot"]:
        #     self.vacuum_locpot = self.entry["vacuum_locpot"]["value"]
        # else:
        #     self.vacuum_locpot = 0

        # if locpot and self.entry["info_primitive"]["vacuum_locpot_primitive"]:
        #     self.vacuum_locpot_primitive = self.entry["info_primitive"]["vacuum_locpot_primitive"]["value"]
        # else:
        #     self.vacuum_locpot_primitive = 0

        self.nn = self.entry["NN"]
        self.proj_eigenvals = db.get_proj_eigenvals(self.entry["task_id"])
        self.eigenvals = db.get_eigenvals(self.entry["task_id"])
        self.show_edges = show_edges
        print("total_mag:{:.3f}".format(self.entry["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]))
        print("cbm:{:.3f}, vbm:{:.3f}, efermi:{:.3f}".format(self.cbm+self.vacuum_locpot,
                                                            self.vbm+self.vacuum_locpot,
                                                            self.efermi+self.vacuum_locpot,
                                                            ))

    def get_candidates(self, kpoint=0, threshold=0.2, select_up=None, select_dn=None):
        # find promising eigenstates with (band_index, [energy, occupation])
        eigenvals = defaultdict(list)
        energy_range = None
        if self.show_edges == "band_edges":
            energy_range = [self.vbm-0.1, self.cbm+0.1]
        # elif self.show_edges == "band_edges":
        #     energy_range = [self.show_edges[0], self.show_edges[1]]

        try:
            for band_idx, i in enumerate(self.eigenvals["1"][kpoint]):
                # (band_index, [energy, occupied])
                if energy_range[1] > (i[0] - self.vacuum_locpot) > energy_range[0]:
                    eigenvals["1"].append((band_idx, i))
        except IndexError:
            print("Threshold of projection is too high!")

        try:
            for band_idx, i in enumerate(self.eigenvals["-1"][kpoint]):
                if energy_range[1] > (i[0] - self.vacuum_locpot) > energy_range[0]:
                    eigenvals["-1"].append((band_idx, i))
        except IndexError:
            print("Threshold of projection is too high!")


        # find promising states if total projection of nn is larger than 0.2
        promising_band = defaultdict(list)
        band_up_proj = []
        for band in eigenvals["1"]:
            total_proj = 0
            for ion_idx in self.nn:
                total_proj += sum(self.proj_eigenvals["1"][kpoint][band[0]][ion_idx])
            if total_proj >= threshold:
                # print("band_index: {}".format(band[0]))
                promising_band["1"].append((band[0], band[1], total_proj))
                sheet_procar = defaultdict(list)
                for idx, o in enumerate(['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']):
                    test = {}
                    for ion_idx in self.nn:
                        test.update({ion_idx: self.proj_eigenvals["1"][kpoint][band[0]][ion_idx][idx]})
                    test.update({"band_index": band[0]}) #
                    test.update({"spin":"up"}) #
                    test.update({"orbital":o})
                    band_up_proj.append(test)

        band_dn_proj = []
        for band in eigenvals["-1"]:
            total_proj = 0
            for ion_idx in self.nn:
                total_proj += sum(self.proj_eigenvals["-1"][kpoint][band[0]][ion_idx])
            if total_proj >= threshold:
                # print("band_index: {}".format(band[0]))
                promising_band["-1"].append((band[0], band[1], total_proj))
                # sheet_procar = defaultdict(list)
                for idx, o in enumerate(['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']):
                    test = {}
                    for ion_idx in self.nn:
                        test.update({ion_idx: self.proj_eigenvals["-1"][kpoint][band[0]][ion_idx][idx]})
                    test.update({"band_index": band[0]}) #
                    test.update({"spin":"dn"}) #
                    test.update({"orbital":o})
                    band_dn_proj.append(test)

        # Total defect states
        up_channel, dn_channel = None, None
        up = defaultdict(tuple)
        dn = defaultdict(tuple)
        for i in promising_band["1"]:
            up[i[0]] = (round(i[1][0], 3), (i[1][1] > 0.9, i[1][1], i[2]))
        for i in promising_band["-1"]:
            dn[i[0]] = (round(i[1][0], 3), (i[1][1] > 0.9, i[1][1], i[2]))

        sheet_up = defaultdict(list)
        sheet_up["band_index"] = list(up.keys())
        sheet_up["energy"] = [i[0] for i in up.values()]
        sheet_up["occupied"] = [i[1][0] for i in up.values()]
        sheet_up["tot_proj"] = [i[1][2] for i in up.values()]
        sheet_up["n_occ_e"] = [i[1][1] for i in up.values()]
        sheet_up["spin"] = ["up" for i in range(len(up.keys()))]

        sheet_down = defaultdict(list)
        sheet_down["band_index"] = list(dn.keys())
        sheet_down["energy"] = [i[0] for i in dn.values()]
        sheet_down["occupied"] = [i[1][0] for i in dn.values()]
        sheet_down["tot_proj"] = [i[1][2] for i in dn.values()]
        sheet_down["n_occ_e"] = [i[1][1] for i in dn.values()]
        sheet_down["spin"] = ["dn" for i in range(len(dn.keys()))]

        df_up = pd.DataFrame(sheet_up).set_index(["band_index"])
        df_dn = pd.DataFrame(sheet_down).set_index(["band_index"])
        df_up_proj = pd.DataFrame(band_up_proj).set_index(["band_index"])
        df_dn_proj = pd.DataFrame(band_dn_proj).set_index(["band_index"])
        if select_up and select_dn:
            up_channel = df_up.loc[select_up, :].sort_index(ascending=False)
            dn_channel = df_dn.loc[select_dn, :].sort_index(ascending=False)
            band_up_proj = df_up_proj.loc[select_up, :].sort_index(ascending=False)
            band_dn_proj = df_dn_proj.loc[select_dn, :].sort_index(ascending=False)

        else:
            up_channel = df_up.sort_index(ascending=False)
            dn_channel = df_dn.sort_index(ascending=False)

        if select_up and select_dn:
            up_levels = dict(zip(up_channel.loc[select_up, "energy"], up_channel.loc[select_up, "occupied"]))
            dn_levels = dict(zip(dn_channel.loc[select_dn, "energy"], dn_channel.loc[select_dn, "occupied"]))
            eng = EnergyLevel(up_levels, dn_levels)
            fig = None
        else:
            up_levels = dict(zip(up_channel.loc[:, "energy"], up_channel.loc[:, "occupied"]))
            dn_levels = dict(zip(dn_channel.loc[:, "energy"], dn_channel.loc[:, "occupied"]))
            eng = EnergyLevel(up_levels, dn_levels)
            fig = None

        if self.show_edges == "band_edges":
            fig = eng.plotting(
                round(self.vbm + self.vacuum_locpot, 3),
                round(self.cbm + self.vacuum_locpot, 3)
            )

        if self.save_fig_path:
            fig.savefig(os.path.join(self.save_fig_path, "defect_states", "{}_{}_{}.defect_states.png".format(
                self.entry["task_id"],
                self.entry["formula_pretty"],
                self.entry["task_label"])))

        state_df = pd.concat([up_channel, dn_channel], ignore_index=True)
        proj_state_df = pd.DataFrame(band_up_proj + band_dn_proj)

        # depict defate states
        ## up
        up_states = state_df.loc[state_df["spin"] == "up"]
        up_dist_from_vbm = up_states["energy"] - (self.vbm + self.vacuum_locpot)
        up_dist_from_vbm = up_dist_from_vbm.round(3)
        up_occ = up_states.loc[:, "n_occ_e"]
        # print(list(up_dist_from_vbm))
        ## dn
        dn_states = state_df.loc[state_df["spin"] == "dn"]
        dn_dist_from_vbm = dn_states["energy"] - (self.vbm + self.vacuum_locpot)
        dn_dist_from_vbm = dn_dist_from_vbm.round(3)
        dn_occ = dn_states.loc[:, "n_occ_e"]

        d_df = pd.DataFrame(
            [
                {
                    "up_from_vbm": list(up_dist_from_vbm),
                    "up_occ": list(up_occ),
                    "dn_from_vbm": list(dn_dist_from_vbm),
                    "dn_occ": list(dn_occ),
                }
            ]
        )
        # highest_level_from_cbm = (self.cbm + self.vacuum_locpot) - highest_level
        # lowest_level_from_vbm = lowest_level - (self.vbm + self.vacuum_locpot)
        # highest_occupied_level_from_cbm = (self.cbm + self.vacuum_locpot) - highest_occupied_level
        # highest_occupied_level_from_vbm = highest_occupied_level - (self.vbm + self.vacuum_locpot)
        # d_df = pd.DataFrame([{"highest_from_cbm": highest_level_from_cbm,
        #                       "lowest_from_vbm": lowest_level_from_vbm,
        #                       "occ_level_from_cbm": highest_occupied_level_from_cbm,
        #                       "occ_level_from_vbm": highest_occupied_level_from_vbm,
        #                       }]).round(3)

        return state_df, proj_state_df, d_df

    def get_degeneracy(self, threshold, energy_tolerence, kpoint=0):
        eigenvals = {"1": [], "-1": []}
        for i in self.eigenvals["1"][kpoint]:
            eigenvals["1"].append((self.eigenvals["1"][kpoint].index(i), i))

        for i in self.eigenvals["-1"][kpoint]:
            eigenvals["-1"].append((self.eigenvals["-1"][kpoint].index(i), i))

        # find promising states if total projection of nn is larger than 0.2
        promising_band = {"1": [], "-1": []}
        for band in eigenvals["1"]:
            total_proj = 0
            for ion_idx in self.nn:
                total_proj += sum(self.proj_eigenvals["1"][kpoint][band[0]][ion_idx])
            if total_proj >= threshold:
                promising_band["1"].append((band[0], band[1], total_proj))

        for band in eigenvals["-1"]:
            total_proj = 0
            for ion_idx in self.nn:
                total_proj += sum(self.proj_eigenvals["-1"][kpoint][band[0]][ion_idx])
            if total_proj >= threshold:
                promising_band["-1"].append((band[0], band[1], total_proj))

        candidate_up_band = np.array([[i[0], i[1][0]] for i in promising_band["1"]])
        diff_can_up_band = np.diff(candidate_up_band[:, 1])/np.average(candidate_up_band)
        index = np.where(diff_can_up_band <= np.min(diff_can_up_band)*energy_tolerence[0])
        for i in index[0]: print("up: %s" % candidate_up_band[i:i+2, :])

        candidate_dn_band = np.array([[i[0], i[1][0]] for i in promising_band["-1"]])
        diff_can_dn_band = np.diff(candidate_dn_band[:, 1])/np.average(candidate_dn_band)
        index = np.where(diff_can_dn_band <= np.min(diff_can_dn_band)*energy_tolerence[1])
        for i in index[0]: print("down: %s" % candidate_dn_band[i:i+2, :])


if __name__ == '__main__':
    pass

