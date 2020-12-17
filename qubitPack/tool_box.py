from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.analysis.local_env import CrystalNN

import numpy as np

from pycdt.core.defectsmaker import ChargedDefectsStructures

from phonopy.phonon.irreps import character_table



def find_cation_anion(structure):
    cation, anion = None, None
    for idx, el in enumerate(list(dict.fromkeys(structure.species))):
        if el.is_metal:
            cation = list(dict.fromkeys(structure.species))[idx].name
        else:
            anion = list(dict.fromkeys(structure.species))[idx].name
    return cation, anion


def modify_vacuum(orig_st, vacuum):
    if vacuum < orig_st.lattice.c:
        print("Please Set Vacuum > original lattice!!")
    ase_offset = AseAtomsAdaptor.get_atoms(orig_st)
    ase_offset.center(vacuum=0.0, axis=2)
    try:
        offset = AseAtomsAdaptor.get_structure(ase_offset).lattice.c
    except Exception as err:
        print(err)
        offset = 0
    ase_atom_obj = AseAtomsAdaptor.get_atoms(orig_st)
    ase_atom_obj.center(vacuum=(vacuum-offset)/2, axis=2)
    return AseAtomsAdaptor.get_structure(ase_atom_obj)


def get_rand_vec(distance): #was 0.001
    # deals with zero vectors.
    vector = np.random.randn(3)
    vnorm = np.linalg.norm(vector)
    return vector / vnorm * distance if vnorm != 0 else get_rand_vec(distance)


def move_site(structure, tgt_sites_idx, displacement_vector, frac_coords=False):
    modified_structure = structure.translate_sites(tgt_sites_idx, displacement_vector, frac_coords=frac_coords)
    return modified_structure


def selective_dyn_site(structure, tgt_sites_idx):
    where = []
    for site_idx in range(len(structure.sites)):
        if site_idx in tgt_sites_idx:
            where.append([True, True, True])
        else:
            where.append([False, False, False])
    poscar = Poscar(structure)
    poscar.selective_dynamics = where
    modified_structure = poscar.structure
    return modified_structure


def defect_from_primitive_cell(orig_st, defect_type, natom, substitution=None, distort=0.002, vacuum_thickness=None):
    """
    defect_type: ["vacancies", 0] / ["substitutions", 0]
    """
    if vacuum_thickness:
        orig_st = modify_vacuum(orig_st, vacuum_thickness)
    else:
        orig_st = orig_st

    defects = ChargedDefectsStructures(orig_st, cellmax=natom, antisites_flag=True).defects
    bulk_st = defects["bulk"]["supercell"]["structure"]

    if defect_type[0] == "bulk":
        bulk_structure = defects[defect_type[0]]["supercell"]["structure"]
        defect_entry = defects[defect_type[0]]
        defect_entry["supercell"].pop("structure")
        return bulk_structure, None, defect_entry, None

    # find NN in defect structure
    def find_nn(defect=defects, defect_type=defect_type):
        defect_st = defect[defect_type[0]][defect_type[1]]["supercell"]["structure"].get_sorted_structure()
        defect_site_in_bulk = defect[defect_type[0]][defect_type[1]]["bulk_supercell_site"]
        defect_entry = defect[defect_type[0]][defect_type[1]]
        defect_entry["supercell"].pop("structure")
        defect_entry["supercell"]["bulk"] = defects["bulk"]["supercell"]["structure"]

        if defect_type[0] == "substitutions":
            defect_site_in_bulk_index = defect_st.index(defect_site_in_bulk)
            print(defect_site_in_bulk, defect_site_in_bulk.to_unit_cell())
            NN = [defect_st.index(defect_st[nn['site_index']])
                  for nn in CrystalNN().get_nn_info(defect_st, defect_site_in_bulk_index)]
            bond_length = [defect_st.get_distance(defect_site_in_bulk_index, NN_index) for NN_index in NN]
            NN = dict(zip(NN, bond_length))
            print("=="*50, "\nBefore distortion: {}".format(NN))
            return defect_st, NN, defect_entry, defect_site_in_bulk_index

        elif defect_type[0] == "vacancies":
            bulk_st = defect["bulk"]["supercell"]["structure"]
            print("vac coord. = {}-{}".format(defect_site_in_bulk, defect_site_in_bulk.to_unit_cell()))
            try:
                defect_site_in_bulk_index = bulk_st.index(defect_site_in_bulk)
            except ValueError:
                defect_site_in_bulk_index = bulk_st.index(defect_site_in_bulk.to_unit_cell())
            NN = [defect_st.index(bulk_st[nn['site_index']])
                  for nn in CrystalNN().get_nn_info(bulk_st, defect_site_in_bulk_index)]
            bond_length = [bulk_st.get_distance(defect_site_in_bulk_index, NN_index['site_index'])
                           for NN_index in CrystalNN().get_nn_info(bulk_st, defect_site_in_bulk_index)]
            NN = dict(zip(NN, bond_length))
            print("=="*50, "\nBefore distortion: {}".format(NN))
            return defect_st, NN, defect_entry, defect_site_in_bulk_index

    defect_st, NN, defect_entry, defect_site_in_bulk_index = find_nn()

    # Move ions around the vacancy randomly
    # add defect into nn_dist for DOS
    for site in NN.keys():
        perturb = get_rand_vec(distort)
        defect_st.translate_sites(site, perturb)
        bulk_st.translate_sites(site, perturb)
    if defect_type[0] == "substitutions":
        NN[defect_site_in_bulk_index] = 0
        bond_length = [defect_st.get_distance(defect_site_in_bulk_index, NN_index) for NN_index in NN]
        NN = dict(zip(NN.keys(), bond_length))
    elif defect_type[0] == "vacancies":
        bond_length = [bulk_st.get_distance(defect_site_in_bulk_index, NN_index['site_index'])
                       for NN_index in CrystalNN().get_nn_info(bulk_st, defect_site_in_bulk_index)]
        NN = dict(zip(NN.keys(), bond_length))

    print("After distortion: {}\n{}".format(NN, "=="*50))

    # To make a substitution to a nn_dist among given element (default None)
    if substitution:
        defect_st.replace(list(NN.keys())[0], substitution)
        print("substitution coord = {}".format(defect_st[list(NN.keys())[0]]))
        defect_entry["complex"] = {"site": defect_st[list(NN.keys())[0]], "site_specie":  substitution}
        defect_entry["defect_type"] = "complex"
        defect_sites_in_bulk = [defect_st.sites[nn] for nn in NN.keys()]
        defect_st.sort()
        NN = [defect_st.index(nn) for nn in defect_sites_in_bulk]

    if not substitution:
        NN = list(NN.keys())

    print("Nearest neighbors = %s" % NN)
    # static_set.write_input(".")
    return defect_st, NN, defect_entry, defect_site_in_bulk_index


class GenDefect:
    def __init__(self, orig_st, defect_type, natom, vacuum_thickness=None, distort=None, sub_on_side=None):
        for site_property in orig_st.site_properties:
            orig_st.remove_site_property(site_property)
        if vacuum_thickness:
            self.orig_st = modify_vacuum(orig_st, vacuum_thickness)
        else:
            self.orig_st = orig_st
        self.defect_type = defect_type
        self.natom = natom
        self.vacuum_thickness = vacuum_thickness
        self.distort = None
        self.defects = ChargedDefectsStructures(self.orig_st, cellmax=natom, antisites_flag=True).defects
        self.bulk_st = self.defects["bulk"]["supercell"]["structure"]
        self.NN_for_sudo_bulk = []

        if (defect_type[0] == "substitutions") or (defect_type[0] == "vacancies"):
            self.defect_entry = self.defects[self.defect_type[0]][self.defect_type[1]]
            self.defect_st = self.defect_entry["supercell"]["structure"].get_sorted_structure()
            self.defect_site_in_bulk = self.defect_entry["bulk_supercell_site"]
            self.defect_site_in_bulk_index = None
            self.NN = None

            if defect_type[0] == "substitutions":
                try:
                    self.defect_site_in_bulk_index = self.defect_st.index(self.defect_site_in_bulk)
                except ValueError:
                    self.defect_site_in_bulk_index = self.defect_st.index(self.defect_site_in_bulk.to_unit_cell())
                self.NN = [self.defect_st.index(self.defect_st[nn['site_index']])
                           for nn in CrystalNN().get_nn_info(self.defect_st, self.defect_site_in_bulk_index)]

            elif defect_type[0] == "vacancies":
                try:
                    self.defect_site_in_bulk_index = self.bulk_st.index(self.defect_site_in_bulk)
                except ValueError:
                    self.defect_site_in_bulk_index = self.bulk_st.index(self.defect_site_in_bulk.to_unit_cell())
                self.NN = [self.defect_st.index(self.bulk_st[nn['site_index']])
                           for nn in CrystalNN().get_nn_info(self.bulk_st, self.defect_site_in_bulk_index)]

            self.nn_dist = dict(before=None, after=None)
            self.nn_dist["before"] = dict(zip([str(idx) for idx in self.NN], range(len(self.NN))))
            print("defect site coord. = orig: {} unit: {}".format(
                self.defect_site_in_bulk, self.defect_site_in_bulk.to_unit_cell()))

            self.defect_entry["supercell"].pop("structure")
            self.defect_entry["supercell"]["bulk"] = self.bulk_st

        elif defect_type[0] == "bulk":
            self.defect_entry = self.defects[self.defect_type[0]]
            self.defect_entry["supercell"].pop("structure")

        else:
            print("!!!Please insert substitutions, vacancies, or bulk!!!")

        if defect_type[0] == "substitutions":
            self.substitutions(distort, sub_on_side)

        elif defect_type[0] == "vacancies":
            self.vacancies(distort, sub_on_side)

    def substitutions(self, distort, substitution):
        bond_length = [self.defect_st.get_distance(self.defect_site_in_bulk_index, NN_index)
                       for NN_index in self.NN]
        bond_length = np.array(bond_length).round(3)

        self.nn_dist["before"] = dict(zip([str(idx) for idx in self.NN], bond_length))

        if substitution:
            self.make_complex(substitution)

        self.NN.append(self.defect_site_in_bulk_index)
        self.nn_dist["before"][str(self.defect_site_in_bulk_index)] = 0
        print("==" * 50, "\nBefore distortion: {}".format(self.nn_dist["before"]))

        if distort:
            self.move_sites(distort)
            bond_length = [self.defect_st.get_distance(self.defect_site_in_bulk_index, NN_index)
                           for NN_index in self.NN]
            bond_length = np.array(bond_length).round(3)
            self.nn_dist["after"] = dict(zip([str(idx) for idx in self.NN], bond_length))
            print("After distortion: {}\n{}".format(self.nn_dist["after"], "==" * 50))

    def vacancies(self, distort, substitution):
        bond_length = [self.bulk_st.get_distance(self.defect_site_in_bulk_index, NN_index['site_index'])
                       for NN_index in CrystalNN().get_nn_info(self.bulk_st, self.defect_site_in_bulk_index)]
        bond_length = np.array(bond_length).round(3)

        self.nn_dist["before"] = dict(zip([str(idx) for idx in self.NN], bond_length))

        if substitution:
            self.make_complex(substitution)

        print("==" * 50, "\nBefore distortion: {}".format(self.nn_dist["before"]))

        if distort:
            sudo_bulk = self.move_sites(distort)
            bond_length = [sudo_bulk.get_distance(self.defect_site_in_bulk_index, NN_index['site_index'])
                           for NN_index in CrystalNN().get_nn_info(sudo_bulk, self.defect_site_in_bulk_index)]
            bond_length = np.array(bond_length).round(3)
            self.nn_dist["after"] = dict(zip([str(idx) for idx in self.NN], bond_length))
            print("After distortion: {}\n{}".format(self.nn_dist["after"], "==" * 50))

    def move_sites(self, distort):
        self.distort = distort
        sudo_bulk = self.bulk_st.copy()
        if self.NN_for_sudo_bulk:
            NN_tot = zip(self.NN, self.NN_for_sudo_bulk)
        else:
            NN_tot = zip(self.NN, self.NN)
        for site, sudo_bulk_site in NN_tot:
            perturb = get_rand_vec(distort)
            self.defect_st.translate_sites([site], perturb, frac_coords=False)
            sudo_bulk.translate_sites([sudo_bulk_site], perturb, frac_coords=False)
        return sudo_bulk

    def make_complex(self, substitution):
        self.defect_entry["complex"] = {"site": [], "site_specie": []}
        for sub, idx in zip(substitution, range(len(substitution))):
            self.defect_st.replace(self.NN[idx], sub)
            print("substitution coord = {}".format(self.NN[idx]))

            self.NN_for_sudo_bulk.append(self.NN[idx])

            self.defect_entry["complex"]["site"].append(self.defect_st[self.NN[idx]])
            self.defect_entry["complex"]["site_specie"].append(sub)

        self.defect_entry["defect_type"] = "complex"

        defect_sites_in_bulk = [self.defect_st[nn] for nn in self.NN]

        self.defect_st.sort()
        if self.defect_type[0] == "substitutions":
            self.defect_site_in_bulk_index = self.defect_st.index(self.defect_site_in_bulk)

        self.NN = [self.defect_st.index(nn) for nn in defect_sites_in_bulk]
        self.nn_dist["before"] = dict(zip([str(idx) for idx in self.NN], self.nn_dist["before"].values()))


def get_good_ir_sites(species, site_syms):
    good_ir_species = []
    good_ir_syms = []
    for specie, site_sym in zip(species, site_syms):
        print(site_sym)
        site_sym = [x for x in site_sym.split(".") if x][0]
        if site_sym == "-4m2":
            site_sym = "-42m"
        if site_sym == "2mm" or site_sym == "m2m":
            site_sym = "mm2"
        if site_sym == "1":
            continue

        irreps = character_table[site_sym][0]["character_table"]
        for irrep, char_vec in irreps.items():
            print(char_vec)
            if char_vec[0] >= 2:
                good_ir_species.append((str(specie)))
                good_ir_syms.append(site_sym)
                break

    return good_ir_species, good_ir_syms