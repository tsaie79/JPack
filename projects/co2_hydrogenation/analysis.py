from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os


os.environ["PATH"] += os.pathsep + "/Library/TeX/texbin/latex"
plt.style.use(["science", "nature", "light"])



class COHP:
    @classmethod
    def analysis(cls):
        from pymatgen.electronic_structure.cohp import CompleteCohp
        from pymatgen.electronic_structure.plotter import CohpPlotter
        from pymatgen.electronic_structure.core import Orbital
        import os
        main_path = "/Users/jeng-yuantsai/Research/project/Collaboration/co2_hydrogenation/calculations/"
        main_path += "H2-Rh_tk26_cohp" #"WTe2" -0.402076
        # main_path += "launcher_2021-05-25-21-13-08-771930" #"WS2" -2.023164
        # main_path += "launcher_2021-05-26-15-51-50-934362" #V_S in WS2
        # main_path += "launcher_2021-06-02-19-57-42-192796" #W_S -0.5; EF= -1.403767
        # main_path += "launcher_2021-06-02-19-57-49-331235" #W_Te 0.5


        os.chdir(main_path)
        # completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename="COHPCAR.lobster", structure_file="POSCAR")
        completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename="COHPCAR.lobster.gz",
                                              structure_file="POSCAR.gz", are_coops=False)
        def summed(lobster_list, plot_label):
            labelist = lobster_list #["929", "961", "967"]
            cp = CohpPlotter(are_coops=False)
            # get a nicer plot label
            plotlabel = plot_label #"W75-W50+W55+W56 bonds"

            cp.add_cohp(plotlabel, completecohp.get_summed_cohp_by_label_list(
                label_list=labelist, divisor=1))
            x = cp.get_plot(integrated=False, invert_axes=False, xlim=(-5, 5), ylim=(-1, 0.5))
            # x.ylim([-5, 5])
            # x.xlim([-2, 0.5])
            plt.show()

        def orbital_resolved(label, orbitals):
            #search for the number of the COHP you would like to plot in ICOHPLIST.lobster (the numbers in COHPCAR.lobster are different!)
            label = str(label)
            cp = CohpPlotter()
            #get a nicer plot label
            plotlabel=str(completecohp.bonds[label]['sites'][0].species_string)+'-'+str(completecohp.bonds[label]['sites'][1].species_string+"{}".format(orbitals))

            cp.add_cohp(plotlabel,completecohp.get_orbital_resolved_cohp(label=label, orbitals=orbitals))
            #check which COHP you are plotting

            print("This is a COHP between the following sites: "+str(completecohp.bonds[label]['sites'][0])+' and '+ str(completecohp.bonds[label]['sites'][1]))

            x = cp.get_plot(integrated=False)
            x.ylim([-5, 5])
            x.xlim([-1.5, 1.5])
            x.show()

        def sum_lobster_orbitals(lobster_list, orbitals, plot_label):
            cp = CohpPlotter(are_coops=False)
            #"W75-W50+W55+W56 bonds"
            cp.add_cohp(plot_label, completecohp.get_summed_cohp_by_label_and_orbital_list(lobster_list, orbitals))
            x = cp.get_plot(integrated=True)
            x.ylim([-10, 6])
            x.show()

        # summed([str(i) for i in [1,2, 30]], "H1-Rh3+H2-Rh3+H1-H2") #WTe2
        # summed([str(i) for i in [2, 30]], "H1-Rh3+H2-Rh3") #WTe2
        # summed([str(i) for i in [1]], "H1-H2") #WTe2
        summed([str(i) for i in [1,2, 30]], "") #WTe2
        summed([str(i) for i in [2, 30]], "") #WTe2
        summed([str(i) for i in [1]], "") #WTe2

    # summed([str(i) for i in [3, 4, 173]], "W1+W6+W7") #Vs in WS2
        # orbital_resolved(13, orbitals=[[5, Orbital.dxy], [5, Orbital.dz2]])
        # antisite_o = [5, Orbital.dz2]
        # for o in [[5, Orbital.dx2], [5, Orbital.dxy], [5, Orbital.dz2], [5, Orbital.dxz], [5, Orbital.dyz]]:
        #     sum_lobster_orbitals([str(i) for i in [929, 961, 967]], [[o, antisite_o],
        #                                                             [o, antisite_o],
        #                                                             [o, antisite_o]], "{}-{}".format(o[1].name, antisite_o[1].name))

COHP.analysis()