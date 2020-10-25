#%% Check structure
from atomate.vasp.database import VaspCalcDb
from pymatgen import Structure
import os


cat = "b2_mg_b2"
path = "/Users/jeng-yuantsai/Research/project/MgB2/calculations/mgb2/{}".format(cat)
db = VaspCalcDb.from_db_file(os.path.join(path, "db.json"))
e = db.collection.find_one({"task_id":16})
st = Structure.from_dict(e["input"]["structure"])
st.to("poscar", os.path.join(path, "structure", "{}_{}.vasp".format(e["formula_pretty"], e["task_id"])))