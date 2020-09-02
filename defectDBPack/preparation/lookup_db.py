#%%
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from atomate.vasp.database import VaspCalcDb
import pandas as pd
from phonopy.phonon.irreps import character_table

db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/code/qubitPack/database_profile/db_dk_local.json")
col = db.collection

filter = {
    "gap_hse_nosoc":{"$gte":1},
    # "ehull":{"$lt":0.4},
    "nkinds":2,
    "magstate": "NM"
}
show = {
    "formula":1,
    "class":1,
    "ehull":1,
    "spacegroup":1,
    "gap_hse_nosoc":1
}

data = []
for e in col.find(filter):
    st = Structure.from_dict(e["structure"])
    sites = SpacegroupAnalyzer(st, symprec=1e-4).get_symmetry_dataset()["std_positions"]
    site_syms = SpacegroupAnalyzer(st, symprec=1e-4).get_symmetry_dataset()["site_symmetry_symbols"]
    point_gp = SpacegroupAnalyzer(st, symprec=1e-4).get_symmetry_dataset()["pointgroup"]
    site_irreps = []
    for specie, site, site_sym in zip(st.species, sites, site_syms):
        print(site_sym)
        site = list(site.round(4))
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
                # site_irreps.append((str(specie), site, irrep, char_vec[0]))
                site_irreps.append(str(specie))
                break




    print(e["uid"])
    print(site_irreps)
    data.append(
        {
            "formula": e["formula"],
            "irreps": site_irreps,
            "spacegroup": e["spacegroup"],
            "pmg_point_gp": point_gp,
            "gap_hse": e["gap_hse"],
            "gap_hse_nosoc": e["gap_hse_nosoc"],
            "gap_nosoc": e["gap_nosoc"],
            "soc_intensity": e["gap_hse_nosoc"] - e["gap_hse"],
            "ehull": e["ehull"],
            "uid":e["uid"]
        }
    )



df = pd.DataFrame(data).round(3).sort_values(["soc_intensity", "ehull"], ascending=True)
df.to_json("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/gap_gt1-binary-NM.json", orient="records", indent=4)
df.to_excel("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/gap_gt1-binary-NM.xlsx", index=False)

#%%
from atomate.vasp.database import VaspCalcDb
import pandas as pd

db = VaspCalcDb.from_db_file("/global/homes/t/tsaie79/code/qubitPack/database_profile/db_mongo_cori.json")
col = db.collection
data = []
for e in col.find():
    data.append(
        {
            "formula": e["formula_pretty"],
            "spacegroup": e["output"]["spacegroup"]["symbol"],
            "bandgap": e["output"]["bandgap"],
            "task_id": e["task_id"]
        }
    )

df = pd.DataFrame(data).sort_values(["bandgap"], ascending=False)
df.to_json("cori_relax_lc.json", orient="records")

#%%
from  monty.serialization import loadfn
import pandas as pd

src = loadfn("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/cori_relax_lc.json")
tgt = loadfn("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/gap_gt1-binary-NM.json")
data = []
for e in src:
    for idx, ee in enumerate(tgt):
        if e["formula"] == ee["formula"] and e["spacegroup"] == ee["spacegroup"]:
            data.append(ee)

df1 = pd.DataFrame(data).to_excel("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/duplicate_c2db.xlsx", index=False)
# df = pd.DataFrame(tgt).round(3).sort_values(["soc_intensity", "ehull"], ascending=True)
# df.to_excel("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/gap_gt1-binary-NM-no-cori_duplicate.xlsx", index=False)
# df.to_json("/Users/jeng-yuantsai/Research/project/defectDB/xlsx/"
#            "gap_gt1-binary-NM-no-cori_duplicate.json", orient="records", indent=4)

