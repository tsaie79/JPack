from fireworks import LaunchPad
from qubitPack.tool_box import *
from my_atomate.workflows.wf_full import get_wf_full_hse
from my_atomate.powerups import *
from atomate.vasp.powerups import *
from atomate.vasp.config import *
from projects.antisiteQubit.wf_defect import ZPLWF


class Defect:
    @classmethod
    def pc(cls):
        pass

    @classmethod
    def standard_defect(cls, defect_type="substitutions", distort=0.0, category="move_z", dz=0.1): #1e-4

        db_name, col_name = "owls", "mx2_antisite_pc"
        col = get_db(db_name, col_name, port=12345).collection
        mx2s = col.find({"task_id":{"$in":[3091]}})
        # 3091: S-W, 3083: Se-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te

        geo_spec = {5* 5 * 3: [20]}
        for mx2 in mx2s:
            pc = Structure.from_dict(mx2["output"]["structure"])
            defect = ChargedDefectsStructures(pc, antisites_flag=True).defects
            cation, anion = find_cation_anion(pc)

            for sub in range(len(defect[defect_type])):
                print(cation, anion)
                # cation vacancy
                #"as_1_{}_on_{}"
                #"vac_1_{}"
                print(defect[defect_type][sub]["name"])
                if "as_1_{}_on_{}".format(cation, anion) not in defect[defect_type][sub]["name"]:
                    continue
                for na, thicks in geo_spec.items():
                    for thick in thicks:
                        gen_defect = GenDefect(
                            orig_st=pc,
                            defect_type=(defect_type, sub),
                            natom=na,
                            vacuum_thickness=thick,
                            distort=distort,
                            sub_on_side=None,
                        )
                        e = get_db("single_photon_emitter", "standard_defect", port=12345).collection.find_one({"task_id": 459})
                        defect_st = Structure.from_dict(e["input"]["structure"])
                        defect_st.translate_sites([25], [0,0, dz], frac_coords=False)
                        wf = get_wf_full_hse(
                            structure=defect_st,
                            task_arg=dict(lcharg=True),
                            charge_states=[-1],
                            gamma_only=False,
                            gamma_mesh=True,
                            nupdowns=[3],
                            task="hse_scf",
                            vasptodb={
                                "dz": dz,
                                "category": category, "NN": gen_defect.NN,
                                "defect_entry": gen_defect.defect_entry,
                                "lattice_constant": "HSE",
                                "perturbed": gen_defect.distort,
                                "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
                            },
                            wf_addition_name="{}:{}".format(na, thick),
                        )

                        # wf = scp_files(
                        #     wf,
                        #     "/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category),
                        #     fw_name_constraint=wf.fws[-1].name,
                        # )
                        # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb",
                        #                     fw_name_constraint="HSE_scf")

                        # wf = remove_todb(wf, fw_name_constraint=wf.fws[0].name)

                        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="PBE_relax")
                        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_relax")
                        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")
                        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_soc")

                        wf = add_modify_incar(wf)
                        wf = set_execution_options(wf, category=category, fworker_name="owls")
                        wf = preserve_fworker(wf)
                        wf.name = wf.name+":dx[{}]".format(gen_defect.distort)
                        print(wf.name)
                        return wf
    @classmethod
    def soc_standard_defect(cls, std_d_tkid, std_d_base_dir=None, category="soc_standard_defect", scp=True, saxis=(0,0,1)):

        std_d_name = "single_photon_emitter"
        std_d_col_name = "standard_defect"
        std_d = get_db(std_d_name, std_d_col_name, port=12345)
        entry = std_d.collection.find_one({"task_id": std_d_tkid})
        std_d_dir = entry["calcs_reversed"][0]["dir_name"].split("/")[-1]
        std_d_base_dir = std_d_base_dir or entry["calcs_reversed"][0]["dir_name"].split("launch")[0]
        std_d_dir = os.path.join(std_d_base_dir, std_d_dir)

        st = Structure.from_dict(entry["output"]["structure"])

        wf = get_wf_full_hse(
            structure=st,
            task_arg={"prev_calc_dir": std_d_dir, "saxis": saxis},
            charge_states=[0],
            gamma_only=False,
            gamma_mesh=True,
            nupdowns=[-1],
            task="hse_soc",
            vasptodb={
                "category": category,
                "NN": entry["NN"],
                "defect_entry": entry["defect_entry"],
                "nonsoc_from": "{}/{}/{}/{}".format(std_d_name, std_d_col_name, entry["task_label"], entry["task_id"])
            },
            wf_addition_name="{}:{}:SOC:{}".format(None, entry["task_id"], saxis), #entry["perturbed"]
        )

        wf = add_modify_incar(wf, {"incar_update": {"LCHARG":False, "LWAVE":True, "LAECHG":False}})
        if scp:
            wf = scp_files(
                wf,
                "/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category),
                fw_name_constraint=wf.fws[-1].name,
            )
        # wf = clear_to_db(wf, fw_name_constraint=wf.fws[0].name)
        # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")

        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category, fworker_name="efrc")
        wf = preserve_fworker(wf)
        return wf

    @classmethod
    def soc_cdft(cls, task, soc_std_d_tkid, nbands, occ=None, soc_std_d_base_dir=None, std_d_base_dir=None,
                 specific_poscar=None, category="soc_cdft", prevent_JT=True):

        soc_std_d = get_db("single_photon_emitter", "soc_standard_defect", port=12345)
        soc_std_d_e = soc_std_d.collection.find_one({"task_id":soc_std_d_tkid})

        if not soc_std_d_base_dir:
            soc_std_d_base_dir = soc_std_d_e["calcs_reversed"][0]["dir_name"].split("launch")[0]

        soc_std_d_path = soc_std_d_e["dir_name"].split("/")[-1]

        std_d_e = soc_std_d_e["nonsoc_from"].split("/")
        std_d = get_db(std_d_e[0], std_d_e[1], port=12345)
        std_d_path = std_d.collection.find_one({"task_id":int(std_d_e[-1])})["dir_name"].split("/")[-1]
        if not std_d_base_dir:
            std_d_base_dir = std_d.collection.find_one({"task_id":int(std_d_e[-1])})["dir_name"].split(":")[-1].split("launch")[0]

        anti_triplet = ZPLWF(os.path.join(soc_std_d_base_dir, soc_std_d_path), None)

        if not occ:
            maj_spin, occ_config = get_lowest_unocc_band_idx(soc_std_d_tkid, soc_std_d, nbands, prevent_JT=prevent_JT,
                                                             second_excite=False)
            occ = occ_config[maj_spin]

        wf = anti_triplet.wf(
            task, 0, up_occupation=occ,
            down_occupation=None, nbands=nbands, gamma_only=True, selective_dyn=None,
            specific_structure=specific_poscar,
            nonsoc_prev_dir=os.path.join(std_d_base_dir, std_d_path)
        )
        print(std_d_base_dir, std_d_path)
        import math
        wf = jmodify_to_soc(wf, nbands=nbands, structure=anti_triplet.structure) #saxis=[1/math.sqrt(2), 1/math.sqrt(2), 0]
        for idx, task in enumerate(task.split("-")):
            if task == "B" or task == "C":
                wf = add_modify_incar(wf, {"incar_update": {"ICHARG":0}}, fw_name_constraint=wf.fws[idx].name)

        wf = add_additional_fields_to_taskdocs(wf, {"cdft_occ": {"up":occ, "dn":None},
                                                    "source_entry":"{}/{}/{}".format(soc_std_d.db_name,
                                                                                          soc_std_d.collection.name,
                                                                                          soc_std_d_tkid)})
        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category)
        wf = preserve_fworker(wf)
        wf.name = wf.name+":SOC"
        return wf

    @classmethod
    def cdft(cls, task, std_d_tkid, nbands, occ=None, std_d_base_dir=None, category="cdft", specific_poscar=None,
             prevent_JT=True):

        std_d = get_db("single_photon_emitter", "standard_defect", port=12345)
        std_d_e = std_d.collection.find_one({"task_id":std_d_tkid})
        std_d_path = std_d_e["dir_name"].split("/")[-1]

        if not std_d_base_dir:
            std_d_base_dir = std_d_e["calcs_reversed"][0]["dir_name"].split("launch")[0]

        cdft = ZPLWF(os.path.join(std_d_base_dir, std_d_path), None)

        up_occ, dn_occ = None, None
        if not occ:
            maj_spin, occ_config = get_lowest_unocc_band_idx(std_d_tkid, std_d, nbands, prevent_JT=prevent_JT,
                                                             second_excite=False)
            if maj_spin == "1" and len(occ_config.keys()) == 1:
                up_occ = occ_config[maj_spin]
                dn_occ = None
            else:
                up_occ = occ_config["1"]
                dn_occ = occ_config["-1"]
        else:
            up_occ = occ["1"]
            dn_occ = occ["-1"]
        print("up_occ: {}, dn_occ:{}".format(up_occ, dn_occ))
        wf = cdft.wf(
            task, 0, up_occupation=up_occ,
            down_occupation=dn_occ, nbands=nbands, gamma_only=True, selective_dyn=None,
            specific_structure=specific_poscar,
            nonsoc_prev_dir=None
        )
        # wf = scp_files(
        #     wf,
        #     dest="/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category),
        # )

        wf = add_additional_fields_to_taskdocs(wf, {"cdft_occ": {"up":up_occ, "dn":dn_occ},
                                                    "source_entry":"{}/{}/{}".format(std_d.db_name,
                                                                                     std_d.collection.name,
                                                                                     std_d)
                                                    })
        print(up_occ, dn_occ)
        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category)
        wf = preserve_fworker(wf)
        return wf

    @classmethod
    def aeps(cls):

        def excited_state_interpolation(soc_d_tkid_a, soc_std_d_base_dir, std_d_base_dir, category):
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                     "my_launchpad.yaml".format(category))))

            soc_d = get_db("single_photon_emitter", "soc_standard_defect", port=12345)
            soc_d_a = soc_d.collection.find_one({"task_id":soc_d_tkid_a})

            soc_cdft = get_db("single_photon_emitter", "soc_cdft", port=12345)
            b_st = Structure.from_dict(soc_cdft.collection.find_one(
                {"chemsys":soc_d_a["chemsys"], "task_label": {"$regex": "B-HSE_scf_soc"}})["input"]["structure"])
            c_st = Structure.from_dict(soc_cdft.collection.find_one(
                {"chemsys":soc_d_a["chemsys"], "task_label": {"$regex": "C-HSE_relax_soc"}})["output"]["structure"])

            dx = np.linspace(0, 2,11)
            sts, info = get_interpolate_sts(
                b_st,
                c_st,
                disp_range=dx,
                output_dir=None
            )

            for frac, st in sts.items():
                wf = cls.soc_cdft("B", soc_d_tkid_a, soc_d_a["input"]["parameters"]["NBANDS"], soc_std_d_base_dir,
                             std_d_base_dir, Poscar(st), category="soc_cdft")

                wf.name += wf.name+"{}".format(frac)
                wf = add_additional_fields_to_taskdocs(wf, {"info":info, "dx": frac})

                wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
                LPAD.add_wf(wf)


        def ground_state_interpolation(soc_d_tkid_a, std_d_base_dir, category="soc_cdft"):
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                     "my_launchpad.yaml".format(category))))

            soc_d = get_db("single_photon_emitter", "soc_standard_defect", port=12345)
            soc_d_a = soc_d.collection.find_one({"task_id":soc_d_tkid_a})

            std_d_e = soc_d_a["nonsoc_from"].split("/")

            soc_cdft = get_db("single_photon_emitter", "soc_cdft", port=12345)
            b_st = Structure.from_dict(soc_cdft.collection.find_one(
                {"chemsys":soc_d_a["chemsys"], "task_label": {"$regex": "B-HSE_scf_soc"}})["input"]["structure"])
            c_st = Structure.from_dict(soc_cdft.collection.find_one(
                {"chemsys":soc_d_a["chemsys"], "task_label": {"$regex": "C-HSE_relax_soc"}})["output"]["structure"])

            dx = np.linspace(-1, 1,11)
            sts, info = get_interpolate_sts(
                b_st,
                c_st,
                disp_range=dx,
                output_dir=None
            )

            for frac, st in sts.items():
                wf = cls.soc_standard_defect(int(std_d_e[-1]), std_d_base_dir, category)
                wf = write_PMGObjects(wf, dict(poscar=Poscar(st)))

                wf = add_additional_fields_to_taskdocs(wf, {"info":info, "dx": frac})
                wf = add_modify_incar(wf, {"incar_update":{"LCHARG":False, "LWAVE":False, "LAECHG":False}})

                wf = set_queue_options(wf, "24:00:00", fw_name_constraint=wf.fws[0].name)
                LPAD.add_wf(wf)

def main():
    def soc_cdft():
        category = "soc_cdft"
        wf = Defect.soc_cdft(
            "B-C-D",
            538,
            960,
            occ=None,
            prevent_JT=False,
            category=category,
        )
        # wf = add_additional_fields_to_taskdocs(wf, {"PS": "prevent_JT"})
        wf = set_execution_options(wf, fworker_name="owls")
        LPAD = LaunchPad.from_file(
            os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                 "my_launchpad.yaml".format(category))))
        LPAD.add_wf(wf)
    def cdft():
        category = "cdft"
        wf = Defect.cdft(
            "B-C-D",
            460,
            475,
            occ=None,
            category=category
        )
        wf = add_additional_fields_to_taskdocs(wf, {"PS": "prevent_JT"})
        wf = set_execution_options(wf, fworker_name="efrc")
        LPAD = LaunchPad.from_file(
            os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                 "my_launchpad.yaml".format(category))))
        LPAD.add_wf(wf)
    def std_defect():
        wf = Defect.standard_defect(distort=0)
        LPAD = LaunchPad.from_file(
            os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                             "my_launchpad.yaml".format("standard_defect"))))
        LPAD.add_wf(wf)

    def soc_std_defect():
        wf = Defect.soc_standard_defect(
            553,
        )
        LPAD = LaunchPad.from_file(
            os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                 "my_launchpad.yaml".format("soc_standard_defect"))))

        LPAD.add_wf(wf)

    def move_z(dz):
        import numpy as np
        wf = Defect.standard_defect(distort=0, dz=dz)
        LPAD = LaunchPad.from_file(
            os.path.expanduser(os.path.join("~", "config/project/antisiteQubit/{}/"
                                                 "my_launchpad.yaml".format("move_z"))))
        LPAD.add_wf(wf)
    for i in [-0.25]:
        move_z(dz=i)
if __name__ == '__main__':
    main()