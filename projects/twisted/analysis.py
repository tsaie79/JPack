#%% get dos
from qubitPack.tool_box import get_db

from pymatgen.electronic_structure.plotter import DosPlotter

wse2 = get_db("twisted", "Bi2Se3")
dos = wse2.get_dos(2)
se_dos = dos.get_element_spd_dos("Se")
bi_dos = dos.get_element_spd_dos("Bi")

pdos = {"Se_s": list(se_dos.values())[0], "Se_p": list(se_dos.values())[1], "Se_d": list(se_dos.values())[2],
        "Bi_s": list(bi_dos.values())[0], "Bi_p": list(bi_dos.values())[1], "Bi_d": list(bi_dos.values())[2]
        }

dsplotter = DosPlotter()
dsplotter.add_dos_dict(pdos)

plt = dsplotter.get_plot(xlim=[-5,5])
plt.show()

#%%
