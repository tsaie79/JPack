from fireworks import LaunchPad
from qubitPack.tool_box import *
from projects.antisiteQubit.wf_defect import get_wf_full_hse
from atomate.vasp.jpowerups import *
from atomate.vasp.powerups import *
from projects.antisiteQubit.wf_defect import ZPLWF

class Defect:
    @classmethod
    def pc(cls):
        pass

    @classmethod
    def standard_defect(cls, defect_type="substitutions", dtort=0, category="standard_defect"): #1e-4

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
                        se_antisite = GenDefect(
                            orig_st=pc,
                            defect_type=(defect_type, sub),
                            natom=na,
                            vacuum_thickness=thick,
                            distort=dtort,
                            sub_on_side=None,
                            move_origin=True
                        )
                        wf = get_wf_full_hse(
                            structure=se_antisite.defect_st,
                            task_arg=dict(lcharg=True),
                            charge_states=[0],
                            gamma_only=False,
                            gamma_mesh=True,
                            nupdowns=[-1],
                            task="hse_relax-hse_scf",
                            vasptodb={
                                "category": category, "NN": se_antisite.NN,
                                "defect_entry": se_antisite.defect_entry,
                                "lattice_constant": "HSE",
                                "perturbed": se_antisite.distort,
                                "PS": "encut320",
                                "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
                            },
                            wf_addition_name="{}:{}".format(na, thick),
                        )

                        wf = scp_files(
                            wf,
                            "/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category),
                            fw_name_constraint=wf.fws[-1].name,
                        )
                        # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb",
                        #                     fw_name_constraint="HSE_scf")

                        wf = clear_to_db(wf, fw_name_constraint=wf.fws[0].name)
                        # wf = clear_to_db(wf, fw_name_constraint=wf.fws[1].name)

                        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="PBE_relax")
                        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_relax")
                        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")
                        # wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_soc")

                        wf = use_gamma_vasp(wf, gamma_vasp_cmd=GAMMA_VASP_CMD)
                        wf = add_modify_incar(wf)
                        wf = preserve_fworker(wf)
                        wf.name = wf.name+":dx[{}]".format(se_antisite.distort)
                        print(wf.name)
                        return wf
    @classmethod
    def soc_standard_defect(cls, std_d_tkid, std_d_base_dir, category="soc_standard_defect", scp=False, saxis=(0,0,1)):

        std_d_name = "single_photon_emitter"
        std_d_col_name = "standard_defect"
        std_d = get_db(std_d_name, std_d_col_name, port=12345)
        entry = std_d.collection.find_one({"task_id": std_d_tkid})
        std_d_dir = entry["calcs_reversed"][0]["dir_name"].split("/")[-1]
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
            wf_addition_name="{}:{}:SOC:{}".format(entry["perturbed"], entry["task_id"], saxis),
        )

        wf = add_modify_incar(wf, {"incar_update":{"LCHARG":False, "LWAVE":True, "LAECHG":False}})
        if scp:
            wf = scp_files(
                wf,
                "/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category),
                fw_name_constraint=wf.fws[-1].name,
            )
        # wf = clear_to_db(wf, fw_name_constraint=wf.fws[0].name)
        # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")

        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category)
        wf = preserve_fworker(wf)
        return wf

    @classmethod
    def soc_cdft(cls, task, soc_std_d_tkid, nbands, soc_std_d_base_dir, std_d_base_dir, specific_poscar=None,
                 category="soc_cdft", secondary=False):


        soc_std_d = get_db("single_photon_emitter", "soc_standard_defect", port=12345)
        soc_std_d_e = soc_std_d.collection.find_one({"task_id":soc_std_d_tkid})
        soc_std_d_path = soc_std_d_e["dir_name"].split("/")[-1]

        std_d_e = soc_std_d_e["nonsoc_from"].split("/")
        std_d = get_db(std_d_e[0], std_d_e[1], port=12345)
        std_d_path = std_d.collection.find_one({"task_id":int(std_d_e[-1])})["dir_name"].split("/")[-1]


        anti_triplet = ZPLWF(os.path.join(soc_std_d_base_dir, soc_std_d_path), None)

        wf = anti_triplet.wf(
            task, 0, up_occupation=get_lowest_unocc_band_idx(soc_std_d_tkid, soc_std_d, nbands, secondary=secondary),
            down_occupation=None, nbands=nbands, gamma_only=True, selective_dyn=None,
            specific_structure=specific_poscar,
            nonsoc_prev_dir=os.path.join(std_d_base_dir, std_d_path)
        )

        wf = jmodify_to_soc(wf, nbands=nbands, structure=anti_triplet.structure)
        for idx, task in enumerate(task.split("-")):
            if task == "B" or task == "C":
                wf = add_modify_incar(wf, {"incar_update": {"ICHARG":0}}, fw_name_constraint=wf.fws[idx].name)

        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category)
        wf = preserve_fworker(wf)
        return wf

    @classmethod
    def cdft(cls, task, std_d_tkid, std_d_base_dir, nbands, category="cdft", specific_poscar=None, secondary=False):

        std_d = get_db("single_photon_emitter", "standard_defect", port=12345)
        std_d_e = std_d.collection.find_one({"task_id":std_d_tkid})
        std_d_path = std_d_e["dir_name"].split("/")[-1]

        anti_triplet = ZPLWF(os.path.join(std_d_base_dir, std_d_path), None)

        maj_occ, minor_occ = get_lowest_unocc_band_idx(std_d_tkid, std_d, nbands, secondary)

        wf = anti_triplet.wf(
            task, 0, up_occupation="328*1 1*0.5 1*0.5 1*1 619*0",
            down_occupation=minor_occ, nbands=nbands, gamma_only=True, selective_dyn=None,
            specific_structure=specific_poscar,
            nonsoc_prev_dir=None
        )
        wf = scp_files(
            wf,
            dest="/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category),
        )

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

if __name__ == '__main__':
    category = "standard_defect"
    wf = Defect.soc_cdft(
        "B-C-D",
        350,
        950,
        "/home/tug03990/work/single_photon_emitter/soc_standard_defect/block_2021-01-17-04-46-01-108825",
        "/home/tug03990/work/single_photon_emitter/standard_defect/block_2021-01-12-17-28-08-335869,"
    )
    wf = add_additional_fields_to_taskdocs(wf, {"PS": "WS2_c3v"})
    LPAD = LaunchPad.from_file(
        os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                             "my_launchpad.yaml".format(category))))
    LPAD.add_wf(wf)