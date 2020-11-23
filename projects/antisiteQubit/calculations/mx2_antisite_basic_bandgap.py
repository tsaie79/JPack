#%% Band Gap
db = VaspCalcDb.from_db_file("/Users/jeng-yuantsai/Research/qubit/calculations/mx2_antisite_basic_bandgap/db.json")
data = {}
for e in db.collection.find({"task_label": "hse gap"}):
    data[e["chemsys"]] = {
        "bandgap": e["output"]["bandgap"],
        "direct": e["output"]["is_gap_direct"],
        "task_id": e["task_id"],
        "db": "mx2_antisite_basic_bandgap"
    }
# WSe2 doesnt have band gap with AEXX = 0.25 (only 0.3 and 0.35 exist)
df = pd.DataFrame(data)
df = df[["S-W", "Se-W", "Te-W", "Mo-S", "Mo-Se", "Mo-Te"]]
df.round(3)