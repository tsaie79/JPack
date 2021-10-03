from fireworks import LaunchPad
from qubitPack.tool_box import *
from my_atomate_jyt.vasp.workflows.wf_full import get_wf_full_hse
from my_atomate_jyt.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    remove_todb,
    bash_scp_files,
    jmodify_to_soc,
    write_PMGObjects
)
from atomate.vasp.powerups import (
    add_tags,
    add_modify_incar,
    add_modify_kpoints,
    add_modify_potcar,
    set_queue_options,
    set_execution_options,
    preserve_fworker,
    clean_up_files
)

from atomate.vasp.config import *
import sys
sys.path.append("/home/tug03990//scripts/JPack")
# sys.path.append("/home/tug03990/code/JPack/")
from projects.antisiteQubit.wf_defect import ZPLWF

def get_st_p(bottom=50, NN=[5,6,0,25]):
    poscar = Structure.from_file("POSCAR.gz")
    contcar =Structure.from_file("CONTCAR.gz")
    data = []
    for st in [poscar, contcar]:
        if not NN:
            NN = e["NN"]
        print(np.round(np.array(st.lattice.abc)/5, 4).tolist())
        print(np.around(np.array([st.get_distance(NN[-1], x) for x in NN+[bottom] ]), 4))
        d = dict(
            zip(
                ["d1", "d2", "d3", "d4", "z1", "z2", "z3", "z4", "a", "b", "c"],
                np.around(
                    np.array(
                        [st.get_distance(NN[-1], x) for x in NN[:-1]+[bottom] ]+
                        [st.sites[NN[-1]].coords[-1] - st.sites[x].coords[-1] for x in NN[:-1]+[bottom]]+
                        np.round(np.array(st.lattice.abc)/5, 4).tolist()
                    ), 3
                )
            )
        )
        data.append(d)
    import  pandas as pd
    df = pd.DataFrame(data)
    return df

def get_st(st, bottom=50, NN=[5,6,0,25]):
    data = []
    for st in [st]:
        if not NN:
            NN = e["NN"]
        print(np.round(np.array(st.lattice.abc)/5, 4).tolist())
        print(np.around(np.array([st.get_distance(NN[-1], x) for x in NN+[bottom] ]), 4))
        d = dict(
            zip(
                ["d1", "d2", "d3", "d4", "z1", "z2", "z3", "z4", "a", "b", "c"],
                np.around(
                    np.array(
                        [st.get_distance(NN[-1], x) for x in NN[:-1]+[bottom] ]+
                        [st.sites[NN[-1]].coords[-1] - st.sites[x].coords[-1] for x in NN[:-1]+[bottom]]+
                        np.round(np.array(st.lattice.abc)/5, 4).tolist()
                    ), 3
                )
            )
        )
        data.append(d)
    import  pandas as pd
    df = pd.DataFrame(data)
    return df

class Host:
    @classmethod
    def standard_pc(cls):
        category = "pc"
        db_name, col_name = "owls", "mx2_antisite_pc"
        col = get_db(db_name, col_name, port=12345).collection
        mx2 = col.find_one({"task_id":{"$in":[3083]}})
        # 3091: S-W, 3083: Se-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te

        wf = get_wf_full_hse(
            structure=Structure.from_dict(mx2["output"]["structure"]),
            task_arg=dict(lcharg=True),
            charge_states=[0],
            gamma_only=False,
            gamma_mesh=True,
            nupdowns=[-1],
            task="hse_relax-hse_scf",
            vasptodb={
                "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
            })

        wf = add_modify_kpoints(wf, {
            "kpoints_update": {
                "kpts": [[9, 9, 1]]
            }
        })

        wf = remove_todb(wf, fw_name_constraint=wf.fws[0].name)

        wf = set_queue_options(wf, "32:00:00", fw_name_constraint="HSE_relax")
        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")

        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category, fworker_name="efrc")
        wf = preserve_fworker(wf)

        LPAD = LaunchPad.from_file(
            os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                 "my_launchpad.yaml".format(category))))
        LPAD.add_wf(wf)

    @classmethod
    def soc_standard_pc(cls):
        std_pc_tkid = 592
        std_pc_base_dir = None
        saxis = (0, 0, 1)
        category = "soc_pc"

        std_pc_name = "single_photon_emitter"
        std_pc_col_name = "pc"
        std_pc = get_db(std_pc_name, std_pc_col_name, port=12345)
        entry = std_pc.collection.find_one({"task_id": std_pc_tkid})
        std_pc_dir = entry["calcs_reversed"][0]["dir_name"].split("/")[-1]
        std_pc_base_dir = std_pc_base_dir or entry["calcs_reversed"][0]["dir_name"].split("launch")[0]
        std_pc_dir = os.path.join(std_pc_base_dir, std_pc_dir)

        st = Structure.from_dict(entry["output"]["structure"])

        wf = get_wf_full_hse(
            structure=st,
            task_arg={"prev_calc_dir": std_pc_dir, "saxis": saxis, "read_wavecar": False},
            charge_states=[0],
            gamma_only=False,
            gamma_mesh=True,
            nupdowns=[-1],
            task="hse_soc",
            vasptodb={
                "category": category,
                "nonsoc_from": "{}/{}/{}/{}".format(std_pc_name, std_pc_col_name, entry["task_label"], entry["task_id"])
            },
            wf_addition_name="{}:{}:SOC:{}".format(None, entry["task_id"], saxis), #entry["perturbed"]
        )

        wf = add_modify_kpoints(wf, {
            "kpoints_update": {
                "kpts": [[9, 9, 1]]
            }
        })
        wf = add_modify_incar(wf, {"incar_update": {"LCHARG":False, "LWAVE":True, "LAECHG":False, "ISYM": 3}})

        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")

        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category, fworker_name="efrc")
        wf = preserve_fworker(wf)

        LPAD = LaunchPad.from_file(
            os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                 "my_launchpad.yaml".format(category))))
        LPAD.add_wf(wf)

    @classmethod
    def std_and_soc(cls, category, scp=False, fworker="owls"):

        db_name, col_name = "owls", "mx2_antisite_pc"
        # db_name, col_name = "single_photon_emitter", "pc"
        col = get_db(db_name, col_name, port=12345).collection
        mx2s = col.find({"task_id":{"$in":[3093, 3102]}})
        # 3091: S-W, 3083: Se-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te

        wfs = []
        for mx2 in mx2s:
            pc = Structure.from_dict(mx2["output"]["structure"])
            wf = get_wf_full_hse(
                structure=pc,
                task_arg=dict(
                    metadata_to_pass={
                        "pc_from": "pc_from",
                        "nosoc_dir": "calcs_reversed.0.dir_name",
                        "nosoc_tkid": "task_id",
                        "nosoc_db": "db",
                        "nosoc_collection": "collection"
                    }),
                charge_states=[0],
                gamma_only=False,
                gamma_mesh=True,
                nupdowns=[-1],
                task="hse_relax-hse_scf-hse_soc",
            )
            wf = add_modify_incar(wf, {"incar_update": {"LCHARG":False, "LWAVE":True, "LAECHG":False}},
                                  fw_name_constraint=wf.fws[2].name)

            wf = add_modify_kpoints(wf, {
                "kpoints_update": {
                    "kpts": [[9, 9, 1]]
                }
            })

            wf = add_modify_incar(wf)

            # wf = remove_todb(wf, fw_name_constraint=wf.fws[0].name)

            wf = add_additional_fields_to_taskdocs(
                wf,
                {
                    #"category": category["scf"],
                    "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
                },
            )

            wf = add_tags(
                wf,
                [
                    {
                        "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
                    }
                ],
            )

            # wf = add_additional_fields_to_taskdocs(
            #     wf,
            #     {
            #         "category": category["soc"]
            #     },
            #     fw_name_constraint=wf.fws[2].name
            # )

            if scp:
                wf = bash_scp_files(
                    wf,
                    dest="/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category[
                                                                                                    "scf"]),
                    port=12346,
                    fw_name_constraint=wf.fws[-1].name,
                    task_name_constraint="VaspToDb"
                )

                wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb",
                                    fw_name_constraint=wf.fws[1].name)

                wf = bash_scp_files(
                    wf,
                    dest="/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category[
                                                                                                    "soc"]),
                    port=12346,
                    fw_name_constraint=wf.fws[2].name,
                    task_name_constraint="VaspToDb"
                )


            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_relax")
            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")
            wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_soc")


            wf = set_execution_options(wf, category=category["scf"], fworker_name=fworker,
                                       fw_name_constraint=wf.fws[0].name)
            wf = set_execution_options(wf, category=category["scf"], fworker_name=fworker,
                                       fw_name_constraint=wf.fws[1].name)
            wf = set_execution_options(wf, category=category["soc"], fworker_name=fworker,
                                       fw_name_constraint=wf.fws[2].name)

            wf = preserve_fworker(wf)
            print(wf.name)

            # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")
            wfs.append(wf)
        return wfs

    @classmethod
    def run(cls):
        def std_and_soc():
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                     "my_launchpad.yaml".format("pc"))))
            wfs = cls.std_and_soc(category={"scf": "pc", "soc": "soc_pc"})
            for wf in wfs:
                print(wf)
                LPAD.add_wf(wf)

        std_and_soc()
class Defect:
    @classmethod
    def standard_defect(cls, defect_type="substitutions", distort=0.0, category="standard_defect", dz=0.0): #1e-4

        db_name, col_name = "owls", "mx2_antisite_pc"
        # db_name, col_name = "single_photon_emitter", "pc"
        col = get_db(db_name, col_name, port=12345).collection
        mx2s = col.find({"task_id":{"$in":[3091, 3097, 3094, 3093, 3102]}})
        # 3091: S-W, 3083: Se-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te

        for mx2 in mx2s:
            pc = Structure.from_dict(mx2["output"]["structure"])
            scaling = find_scaling_for_2d_defect(pc, 15)[0]
            area = scaling[0]*scaling[1]
            geo_spec = {area*pc.num_sites: [20]}

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
                        # special_st = Structure.from_file("/home/tug03990/wse0_456.vasp")
                        # move_site(special_st, [25], dz)
                        wf = get_wf_full_hse(
                            structure=gen_defect.defect_st,
                            task_arg=dict(lcharg=True),
                            charge_states=[0],
                            gamma_only=False,
                            gamma_mesh=True,
                            nupdowns=[2],
                            task="hse_relax-hse_scf",
                            vasptodb={
                                "category": category, "NN": gen_defect.NN,
                                "defect_entry": gen_defect.defect_entry,
                                "lattice_constant": "HSE",
                                # "perturbed": {"move_z": dz},
                                "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
                                "PS": "full_HSE_from_pc_to_defect",
                            },
                            wf_addition_name="{}:{}".format(na, thick),
                        )

                        wf = add_tags(wf,
                            [
                                {
                                    "PS": "full_HSE_from_pc_to_defect",
                                    "category": category, "NN": gen_defect.NN,
                                    "lattice_constant": "HSE",
                                    # "perturbed": None,
                                    "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
                                }
                            ]
                        )

                        # wf = bash_scp_files(
                        #     wf,
                        #     dest="/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category),
                        #     port=12346,
                        #     fw_name_constraint=wf.fws[-1].name,
                        #     task_name_constraint="VaspToDb"
                        # )

                        # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb",
                        #                     fw_name_constraint="HSE_scf")

                        wf = remove_todb(wf, fw_name_constraint=wf.fws[0].name)

                        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_relax")
                        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")

                        wf = add_modify_incar(wf)
                        wf = set_execution_options(wf, category=category, fworker_name="owls")
                        wf = preserve_fworker(wf)
                        wf.name = wf.name+":dx[{}]".format(gen_defect.distort)
                        print(wf.name)
                        return wf, category


    @classmethod
    def soc_standard_defect(cls, std_d_tkid, std_d_base_dir=None,
                            category="soc_standard_defect", scp=False, saxis=(0,0,1), fworker="owls"):

        std_d_name = "owls"
        std_d_col_name = "mx2_antisite_basic_aexx0.25_final"
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
                # "defect_entry": entry["defect_entry"],
                "nonsoc_from": "{}/{}/{}/{}".format(std_d_name, std_d_col_name, entry["task_label"], entry["task_id"])
            },
            wf_addition_name="{}:{}:SOC:{}".format(None, entry["task_id"], saxis), #entry["perturbed"]
        )

        wf = add_modify_incar(wf, {"incar_update": {"LCHARG":False, "LWAVE":True, "LAECHG":False}})
        if scp:
            wf = bash_scp_files(
                wf,
                dest="/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category),
                port=12346,
                fw_name_constraint=wf.fws[-1].name,
                task_name_constraint="VaspToDb"
            )
        # wf = clear_to_db(wf, fw_name_constraint=wf.fws[0].name)
        # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")

        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category, fworker_name=fworker)
        wf = preserve_fworker(wf)
        return wf, category

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
        # std_d = get_db(soc_std_d_e["prev_fw_db"], soc_std_d_e["prev_fw_collection"], port=12345)
        std_d_path = std_d.collection.find_one({"task_id":int(std_d_e[-1])})["dir_name"].split("/")[-1]
        # std_d_path = std_d.collection.find_one({"task_id":int(soc_std_d_e["prev_fw_taskid"])})["dir_name"].split(
        #     "/")[-1]
        if not std_d_base_dir:
            std_d_base_dir = std_d.collection.find_one({"task_id":int(std_d_e[-1])})["dir_name"].split(":")[-1].split("launch")[0]
            #std_d_base_dir = std_d.collection.find_one({"task_id":int(std_d_e[-1])})["dir_name"].split(":")[
            # -1].split("launch")[0]

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
        return wf, category


    @classmethod
    def new_soc_cdft(cls, task, soc_std_taskid, nbands=None, occ=None, soc_std_base_dir=None, std_base_dir=None,
                 specific_poscar=None, category="soc_cdft", prevent_JT=True):

        soc_std_db = get_db("single_photon_emitter", "soc_standard_defect", port=12345)
        soc_std_entry = soc_std_db.collection.find_one({"task_id":soc_std_taskid})

        soc_std_path = soc_std_entry["calcs_reversed"][0]["dir_name"]
        if not soc_std_base_dir:
            soc_std_base_dir = soc_std_path.split("launch")[0]
        soc_std_launchdir = soc_std_path.split("/")[-1]

        std_taskid, std_db_name, std_collection = \
            soc_std_entry["prev_fw_taskid"], soc_std_entry["prev_fw_db"], soc_std_entry["prev_fw_collection"]
        std_db = get_db(std_db_name, std_collection, port=12345)
        std_entry = std_db.collection.find_one({"task_id": std_taskid})
        std_path = std_entry["calcs_reversed"][0]["dir_name"]
        if not std_base_dir:
            std_base_dir = std_path.split("launch")[0]
        std_launchdir = std_path.split("/")[-1]

        anti_triplet = ZPLWF(os.path.join(soc_std_base_dir, soc_std_launchdir), None)
        if not nbands:
            nbands = soc_std_entry["input"]["parameters"]["NBANDS"]
        if not occ:
            maj_spin, occ_config = get_lowest_unocc_band_idx(
                soc_std_taskid, soc_std_db, nbands, prevent_JT=prevent_JT, second_excite=False)
            occ = occ_config[maj_spin]

        wf = anti_triplet.wf(
            task, 0, up_occupation=occ,
            down_occupation=None, nbands=nbands, gamma_only=True, selective_dyn=None,
            specific_structure=specific_poscar,
            nonsoc_prev_dir=os.path.join(std_base_dir, std_launchdir)
        )

        import math
        wf = jmodify_to_soc(wf, nbands=nbands, structure=anti_triplet.structure) #saxis=[1/math.sqrt(2), 1/math.sqrt(2), 0]
        for idx, task in enumerate(task.split("-")):
            if task == "B" or task == "C":
                wf = add_modify_incar(wf, {"incar_update": {"ICHARG":0}}, fw_name_constraint=wf.fws[idx].name)

        wf = add_additional_fields_to_taskdocs(
            wf,
            {
                "cdft_occ": {"up":occ, "dn":None},
                "source_entry":"{}/{}/{}".format(soc_std_entry["db"],
                                                 soc_std_entry["collection"],
                                                 soc_std_entry["task_id"]),
                "prev_fw_taskid": soc_std_entry["task_id"],
                "prev_fw_db": soc_std_entry["db"],
                "prev_fw_collection": soc_std_entry["collection"]
            }
        )
        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category)
        wf = preserve_fworker(wf)
        wf.name = wf.name+":SOC"
        return wf, category

    @classmethod
    def std_and_soc(cls, category, defect_type="substitutions", distort=0.0, dz=0.0, scp=False, fworker="owls"):

        db_name, col_name = "owls", "mx2_antisite_pc"
        # db_name, col_name = "single_photon_emitter", "pc"
        col = get_db(db_name, col_name, port=12345).collection
        mx2s = col.find({"task_id":{"$in":[3091]}})
        # 3091: S-W, 3083: Se-W, 3093: Te-W, 3097:Mo-S, 3094: Mo-Se, 3102:Mo-Te

        wfs = []
        for mx2 in mx2s:
            pc = Structure.from_dict(mx2["output"]["structure"])
            scaling = find_scaling_for_2d_defect(pc, 15)[0]
            area = scaling[0]*scaling[1]
            geo_spec = {area*pc.num_sites: [20]}

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
                            standardize_st=True
                        )

                        a = get_st(gen_defect.defect_st)
                        print(a)
                        # special_st = Structure.from_file("/home/tug03990/wse0_456.vasp")
                        # move_site(special_st, [25], dz)
                        # gen_defect.defect_st = Structure.from_file(
                        #     os.path.join("/gpfs/scratch/tug03990/single_photon_emitter/standard_defect/"
                        #                  "block_2021-01-11-17-39-35-861559/launcher_2021-09-23-05-15-51-873176",
                        #                  "CONTCAR.gz"))

                        wf = get_wf_full_hse(
                            structure=gen_defect.defect_st,
                            task_arg=dict(
                                metadata_to_pass={
                                    "pc_from": "pc_from",
                                    "NN": "NN",
                                    "defect_entry": "defect_entry",
                                    "nosoc_dir": "calcs_reversed.0.dir_name",
                                    "nosoc_tkid": "task_id",
                                    "nosoc_db": "db",
                                    "nosoc_collection": "collection"
                                }),
                            charge_states=[0],
                            gamma_only=False,
                            gamma_mesh=True,
                            nupdowns=[0],
                            task="hse_relax-hse_scf-hse_soc",
                            wf_addition_name="{}:{}".format(na, thick),
                        )

                        wf = add_modify_incar(wf, {"incar_update": {"EDIFF": 1e-5, "EDIFFG": -0.02}},
                                              fw_name_constraint=wf.fws[0].name)
                        wf = add_modify_incar(wf, {"incar_update": {"LCHARG":False, "LWAVE":True, "LAECHG":False}},
                                              fw_name_constraint=wf.fws[2].name)
                        wf = add_modify_incar(wf)

                        # wf = remove_todb(wf, fw_name_constraint=wf.fws[0].name)

                        wf = add_additional_fields_to_taskdocs(
                            wf,
                            {
                                #"category": category["scf"],
                                "NN": gen_defect.NN,
                                "defect_entry": gen_defect.defect_entry,
                                "lattice_constant": "HSE",
                                "perturbed": {"distort": distort},
                                "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
                                "PS": "test_conv",
                            },
                            fw_name_constraint=wf.fws[1].name
                        )

                        wf = add_tags(
                            wf,
                            [
                                {
                                    # "PS": "full_HSE_from_pc_to_defect",
                                    # "category": category, "NN": gen_defect.NN,
                                    # "lattice_constant": "HSE",
                                      # "perturbed": None,
                                    "pc_from": "{}/{}/{}".format(db_name, col_name, mx2["task_id"]),
                                }
                            ],
                        )

                        # wf = add_additional_fields_to_taskdocs(
                        #     wf,
                        #     {
                        #         "category": category["soc"]
                        #     },
                        #     fw_name_constraint=wf.fws[2].name
                        # )

                        if scp:
                            wf = bash_scp_files(
                                wf,
                                dest="/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category[
                                                                                                                "scf"]),
                                port=12346,
                                fw_name_constraint=wf.fws[-2].name,
                                task_name_constraint="VaspToDb"
                            )

                            wf = bash_scp_files(
                                wf,
                                dest="/home/jengyuantsai/Research/projects/single_photon_emitter/{}".format(category[
                                                                                                                "soc"]),
                                port=12346,
                                fw_name_constraint=wf.fws[-1].name,
                                task_name_constraint="VaspToDb"
                            )


                        wf = set_queue_options(wf, "32:00:00", fw_name_constraint="HSE_relax")
                        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_scf")
                        wf = set_queue_options(wf, "24:00:00", fw_name_constraint="HSE_soc")


                        wf = set_execution_options(wf, category=category["scf"], fworker_name=fworker,
                                                   fw_name_constraint=wf.fws[0].name)
                        wf = set_execution_options(wf, category=category["scf"], fworker_name=fworker,
                                                   fw_name_constraint=wf.fws[1].name)
                        wf = set_execution_options(wf, category=category["soc"], fworker_name=fworker,
                                                   fw_name_constraint=wf.fws[2].name)

                        wf = preserve_fworker(wf)
                        wf.name = wf.name+":dx[{}]".format(gen_defect.distort)
                        print(wf.name)

                        # wf = clean_up_files(wf, files=["*"], task_name_constraint="VaspToDb")
                        wfs.append(wf)
        return wfs



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
                                                                                     std_d_tkid)
                                                    })
        print(up_occ, dn_occ)
        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category)
        wf = preserve_fworker(wf)
        return wf, category


    @classmethod
    def new_cdft(cls, task, std_taskid, nbands=None, occ=None, std_base_dir=None, category="cdft", specific_poscar=None,
             prevent_JT=True):

        std_db = get_db("single_photon_emitter", "standard_defect", port=12345)
        # std_db = get_db("owls", "mx2_antisite_basic_aexx0.25_final", port=12345)
        std_entry = std_db.collection.find_one({"task_id":std_taskid})

        std_path = std_entry["calcs_reversed"][0]["dir_name"]
        std_launchdir = std_path.split("/")[-1]
        if not std_base_dir:
            std_base_dir = std_path.split("launch")[0]

        cdft = ZPLWF(os.path.join(std_base_dir, std_launchdir), None)

        if not nbands:
            nbands = std_entry["input"]["parameters"]["NBANDS"]

        up_occ, dn_occ = None, None
        if not occ:
            maj_spin, occ_config = get_lowest_unocc_band_idx(std_taskid, std_db, nbands, prevent_JT=prevent_JT,
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

        wf = add_additional_fields_to_taskdocs(
            wf,
            {
                "cdft_occ": {"up":up_occ, "dn":dn_occ},
                "source_entry":"{}/{}/{}".format(std_entry["db"],
                                                 std_entry["collection"],
                                                 std_entry["task_id"]),
                "prev_fw_taskid": std_entry["task_id"],
                "prev_fw_db": std_entry["db"],
                "prev_fw_collection": std_entry["collection"]
        })
        print(up_occ, dn_occ)
        wf = add_modify_incar(wf)
        wf = set_execution_options(wf, category=category)
        wf = preserve_fworker(wf)
        return wf, category

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
    @classmethod
    def run(cls):
        def soc_cdft():
            wf, cat = cls.soc_cdft(
                "B-C-D",
                604,
                960,
                occ=None, #"656*1 1*0 1*1 1*1 291*0",
                prevent_JT=True,
                category="soc_cdft",
            )
            wf = add_additional_fields_to_taskdocs(wf, {"PS": "0.5-0.5-1"})
            wf = set_execution_options(wf, fworker_name="owls")
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                     "my_launchpad.yaml".format(cat))))
            LPAD.add_wf(wf)

        def new_soc_cdft():
            wf, cat = cls.new_soc_cdft(
                "B-C-D",
                633, #[645, 646, 633]
                None,
                occ="656*1 1*1 1*0 1*1 301*0", #"604*1 1*0 1*1 1*1 303*0",
                prevent_JT=True,
                category="soc_cdft",
                soc_std_base_dir="/home/tug03990/work/single_photon_emitter/soc_standard_defect/block_2021-01-17-04"
                                 "-46-01-108825",
                std_base_dir="/home/tug03990/work/single_photon_emitter/standard_defect/"
                             "block_2021-01-12-17-28-08-335869"
            )

            occ = "1-0-1"
            # wf = add_modify_potcar(wf, {"potcar_symbols": {"Mo": "Mo_pv"}})
            wf = add_additional_fields_to_taskdocs(wf, {"PS": occ})
            wf.name += ":{}".format(occ)
            # wf = add_modify_incar(wf, {"incar_update": {"EDIFF": 1E-4}})
            wf = set_execution_options(wf, fworker_name="efrc")
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                     "my_launchpad.yaml".format(cat))))
            print(wf)
            LPAD.add_wf(wf)

        def cdft():
            wf, cat = cls.cdft(
                "B-C-D",
                264,
                480,
                occ=None,
                category="cdft",
                prevent_JT=False
            )
            wf = add_additional_fields_to_taskdocs(wf, {"PS": "perturbed"})
            wf = set_execution_options(wf, fworker_name="owls")
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                     "my_launchpad.yaml".format(cat))))
            LPAD.add_wf(wf)

        def new_cdft():
            for i in [642, 640, 643, 632]:
                wf, cat = cls.new_cdft(
                    "B-C-D",
                    i,
                    None,
                    occ=None,
                    category="cdft",
                    prevent_JT=True
                )
                # wf = add_modify_potcar(wf, {"potcar_symbols": {"Mo": "Mo_pv"}})
                wf = add_additional_fields_to_taskdocs(wf, {"PS": "0.5-0.5-1"})
                wf = set_execution_options(wf, fworker_name="owls")
                LPAD = LaunchPad.from_file(
                    os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                         "my_launchpad.yaml".format(cat))))
                LPAD.add_wf(wf)

        def std_defect():
            wf, cat = cls.standard_defect(distort=0)
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                 "my_launchpad.yaml".format(cat))))

            # wf = add_modify_potcar(wf, {"potcar_symbols": {"W": "W"}})
            LPAD.add_wf(wf)

        def soc_std_defect(tkid=3285):
            wf, cat = cls.soc_standard_defect(
                tkid,
                scp=False,
                fworker="efrc"
            )
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                     "my_launchpad.yaml".format(cat))))

            wf = add_modify_potcar(wf, {"potcar_symbols": {"Mo": "Mo_pv"}})
            wf = add_modify_incar(wf, {"incar_update": {"ENCUT": 320}})
            wf = add_additional_fields_to_taskdocs(wf, {
                "prev_fw_taskid": tkid,
                "prev_fw_db": "owls",
                "prev_fw_collection": "mx2_antisite_basic_aexx0.25_final"
            })
            wf.name += ":owls/3281"
            LPAD.add_wf(wf)

        def std_soc_defect():
            cat = "standard_defect"
            wfs = cls.std_and_soc(distort=0.001, fworker="efrc", scp=False,
                                 category=dict(scf="standard_defect", soc="soc_standard_defect"))
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/single_photon_emitter/{}/"
                                                     "my_launchpad.yaml".format(cat))))
            for wf in wfs:
                print(wf)
                LPAD.add_wf(wf)

        def move_z():
            import numpy as np
            wf, cat = cls.standard_defect(distort=0, dz=0.2)
            LPAD = LaunchPad.from_file(
                os.path.expanduser(os.path.join("~", "config/project/antisiteQubit/{}/"
                                                     "my_launchpad.yaml".format(cat))))
            LPAD.add_wf(wf)

        new_soc_cdft()


if __name__ == '__main__':
    Defect.run()


