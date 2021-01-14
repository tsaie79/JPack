#%% IMPORT
from atomate.vasp.database import VaspCalcDb

db_dict = {
    "host": "localhost",
    "port": 1234,
    "database": "single_photon_emitter",
    "collection": "standard_defect",
    "admin_user": "Jeng",
    "admin_password": "qimin",
    "readonly_user": "Jeng_ro",
    "readonly_password": "qimin",
    "aliases": {}
}
def db(collection):
    return VaspCalcDb(host="localhost", port=1234, database="single_photon_emitter",
                collection=collection, user="Jeng_ro", password="qimin", authsource="single_photon_emitter")


#%% WS2 with C3V
from pymatgen import Structure

col = db("cdft").collection

b = col.find_one({"task_id": 23})
c = col.find_one({"task_id": 24})
d = col.find_one({"task_id": 25})


b_st = Structure.from_dict(b["input"]["structure"])
c_st = Structure.from_dict(c["output"]["structure"])

E_ab = b["output"]["energy"] - b["source"]["total_energy"]
E_bc = c["output"]["energy"] - b["output"]["energy"]
E_cd = d["output"]["energy"] - c["output"]["energy"]
E_da = d["output"]["energy"] - d["source"]["total_energy"]

#%% interpolate bewteen excited eq and ground eq
import numpy as np
import os

def get_init_fin(i_st, f_st, disp_range=np.linspace(0, 1, 11), output_dir=None):
    '''
    '''
    # A. Alkauskas, Q. Yan, and C. G. Van de Walle, Physical Review B 90, 27 (2014)
    struct_i, sorted_symbols = i_st, i_st.symbol_set
    struct_f, sorted_symbols = f_st, f_st.symbol_set
    delta_R = struct_f.frac_coords - struct_i.frac_coords
    delta_R = (delta_R + 0.5) % 1 - 0.5

    lattice = struct_i.lattice.matrix #[None,:,:]
    delta_R = np.dot(delta_R, lattice)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Poscar(struct_i).write_file('disp_dir/POSCAR_i'.format(output_dir))
    # Poscar(struct_f).write_file('disp_dir/POSCAR_f'.format(output_dir))


    masses = np.array([spc.atomic_mass for spc in struct_i.species])
    delta_Q2 = masses[:,None] * delta_R ** 2

    print('Delta_Q^2: {:3}'.format(np.sqrt(delta_Q2.sum())))

    # print(delta_Q2.shape)
    # print(struct_i.lattice, struct_i.species, struct_i.cart_coords, )

    for frac in disp_range:
        disp = frac * delta_R
        struct = Structure(struct_i.lattice, struct_i.species,
                           struct_i.cart_coords + disp,
                           coords_are_cartesian=True)
        struct.to("vasp", '{0}/POSCAR_{1:03d}.vasp'.format(output_dir, int(np.rint(frac*10))))

get_init_fin(b_st, c_st, output_dir='/Users/jeng-yuantsai/Research/project/single_photon/calculations/cdft/WS2_c3v/interpolate_st')
