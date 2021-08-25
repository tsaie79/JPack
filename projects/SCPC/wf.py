from qubitPack.tool_box import get_db

from pymatgen.io.vasp.inputs import Structure

from atomate.vasp.fireworks import LepsFW
from atomate.vasp.workflows import get_wf
from atomate.vasp.powerups import *

from fireworks import LaunchPad


DB_SCAN_2D_MAT = get_db("Scan2dMat", "calc_data", port=12347).collection
DB_C2DB = get_db("2dMat_from_cmr_fysik", "2dMaterial_v1", port=12347)
DB_SCPC = get_db("SCPC", "test", port=12347)
print(DB_C2DB.db_name, DB_C2DB.host)

class GetEpsilon:
    def __init__(self, db, entry_filter, lpad):
        self.entry = db.collection.find_one(entry_filter)
        self.lpad = LaunchPad(lpad.host, lpad.port, lpad.db_name, lpad.user, lpad.password)
        self.structure = Structure.from_dict(self.entry["structure"])

    def workflow(self):
        wf = get_wf(self.structure, "wf_yaml/dielectric.yaml")
        wf = preserve_fworker(wf)
        wf = set_execution_options(wf, category="test", fworker_name="owls")
        self.lpad.add_wf(wf)





def main():
    lpad = DB_SCPC
    test = GetEpsilon(DB_C2DB, {"uid": "BN-BN-NM"}, lpad)
