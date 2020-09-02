from qubitPack.workflow.tool_box import defect_from_primitive_cell


class ToolBox:
    @classmethod
    def _find_cation_anion(cls, structure):
        cation, anion = None, None
        for idx, el in enumerate(list(dict.fromkeys(structure.species))):
            if el.is_metal:
                cation = list(dict.fromkeys(structure.species))[idx].name
            else:
                anion = list(dict.fromkeys(structure.species))[idx].name
        return cation, anion







def MX2_anion_vacancy():
    def special_treatment_to_structure(st, option, tgt_st=None, mv=None, nn=None):
        # mv antisite in z
        if option == "mv_z":
            # mv_z = -0.5 #Te
            # mv_z = 0.2 #Se
            # mv_z = 0.286 #S
            st.translate_sites(nn[-1], [0, 0, mv], frac_coords=False)
            return st

        # selective dynamics
        elif option == "selective_dynamics":
            nn.pop(-1)
            where = []
            print(nn)
            for i in range(len(st.sites)):
                if i in nn:
                    where.append([True, True, True])
                else:
                    where.append([False, False, False])
            poscar = Poscar(st)
            poscar.selective_dynamics = where
            st = poscar.structure
            return st

        elif option == "modify_lattice":
            tgt_lattice = tgt_st.lattice
            print(tgt_lattice)
            print(st.lattice)
            st.modify_lattice(tgt_lattice)
            print(st.lattice)
            return st

    lpad = LaunchPad.from_file("/home/tug03990/config/my_launchpad.efrc.yaml")
    col = VaspCalcDb.from_db_file("/home/tug03990/PycharmProjects/my_pycharm_projects/database/db_config/"
                                  "db_dk_local.json").collection

    mx2s = col.find(
        {"class": {"$in": ["TMDC-T", "TMDC-H", "TMDC-T'"]},
         "gap_hse": {"$gt": 1},
         "ehull": {"$lt": 0.3},
         "magstate": "NM",
         # "formula": {"$in":["WTe2"]}
         # "spacegroup": "P-6m2",
         # "formula": {"$in": ["MoS2", "MoSe2",  "MoTe2", "WS2", "WSe2", "WTe2"]}
         "formula": {"$in": ["MoTe2"]}
         }
    )


    # geo_spec = {5 * 5 * 3: [20, 30, 40], 6 * 6 * 3: [20, 30, 40]}
    geo_spec = {5 * 5 * 3: [20]}
    aexx = 0.35
    for mx2 in mx2s:
        pc = Structure.from_dict(mx2["structure"])
        defect = ChargedDefectsStructures(pc, antisites_flag=False).defects
        cation, anion = _find_cation_anion(pc)

        for vac in range(len(defect["vacancies"])):
            print(cation, anion)
            if "{}".format(anion) not in defect["vacancies"][vac]["name"]:
                continue
            for na, thicks in geo_spec.items():
                for thick in thicks:
                    se_antisite = DefectWF(pc, natom=na, vacuum_thickness=thick, substitution=None,
                                           antisite=False, type_vac=vac, bulk=False, distort=0.02)
                    # se_antisite.defect_st = Structure.from_file(os.path.join("/gpfs/work/tug03990/mx2_antisite_basic/"
                    #                                      "block_2020-05-28-17-35-55-920859/"
                    #                                      "launcher_2020-05-28-17-35-58-858406", "CONTCAR.relax2.gz"))
                    # st = special_treatment_to_structure(
                    #     se_antisite.defect_st,
                    #     "mv_z",
                    #     nn=se_antisite.NN,
                    #     mv=mv_z
                    # )
                    # se_antisite.defect_st = special_treatment_to_structure(st, "selective_dynamics",
                    #                                                        nn=se_antisite.NN)

                    wf = se_antisite.hse_scf_wf(
                        charge_states=[0,0], gamma_only=True, dos_hse=False, nupdown_set=[2,0],
                        defect_type={"bulk": pc.composition, "c2db_uid": mx2["uid"]},
                        task_info="{}:{}".format(na, thick),
                        encut=320,
                        include_hse_relax=True
                    )
                    kpoints_kwarg = {
                        'comment': 'ws2',
                        "style": "G",
                        "num_kpts": 0,
                        'kpts': [[2, 2, 1]],
                        'kpts_weights': None,
                        'kpts_shift': (0, 0, 0),
                        'coord_type': None,
                        'labels': None,
                        'tet_number': 0,
                        'tet_weight': 0,
                        'tet_connections': None
                    }

                    # wf = add_modify_kpoints(wf, {"kpoints_update": kpoints_kwarg})
                    wf = add_additional_fields_to_taskdocs(
                        wf,
                        {"perturbed": se_antisite.distort, "wf":[fw.name for fw in wf.fws]}
                    )
                    wf = add_modify_incar(wf, {"incar_update":{"NCORE":4, "NSW":100}}, "PBE_relax")
                    wf = add_modify_incar(wf, {"incar_update":{"NCORE":5, "AEXX":aexx}}, "HSE_relax")
                    wf = add_modify_incar(wf, {"incar_update":{"NCORE":4, "LWAVE": True, "LCHARG":False,
                                                               "AEXX":aexx}}, "HSE_scf")
                    wf = set_queue_options(wf, "02:00:00", fw_name_constraint="HSE_scf")
                    wf = set_queue_options(wf, "10:00:00", fw_name_constraint="HSE_relax")
                    wf = set_queue_options(wf, "02:00:00", fw_name_constraint="PBE_relax")
                    # related to directory
                    wf = set_execution_options(wf, category="mx2_anion_vacancy")
                    wf.name = wf.name+":dx[{}]".format(se_antisite.distort)
                    lpad.add_wf(wf)