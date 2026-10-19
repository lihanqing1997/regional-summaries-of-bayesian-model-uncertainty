from pathlib import Path
import csv,json,collections,datetime
root=Path(__file__).resolve().parent
def readjson(path):
    try:return json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError,ValueError):return {}
def atomic(path,data):
    temp=path.with_suffix(path.suffix+".writing")
    temp.write_text(json.dumps(data,indent=2,allow_nan=False),encoding="utf-8")
    temp.replace(path)
jobs=list(csv.DictReader(open(root/"config/jobs.csv",newline="")))
conditions=list(csv.DictReader(open(root/"config/conditions.csv",newline="")))
outcomes={};failures={};completed=0
case_for={j["job"]:j["case_id"] for j in jobs}
for job in jobs:
    folder=root/"results"/job["job"]
    for row in readjson(folder/"outcomes.json").get("outcomes",[]):outcomes[row["condition"]]=row
    if (folder/"complete.json").exists():completed+=1
    for kind in ["worker_failure.json","supervisor_failure.json"]:
        failure=readjson(folder/kind)
        if failure:failures[job["job"]]=failure
rows=[];groups=collections.defaultdict(list)
for c in conditions:
    row={**c,"case_id":case_for[c["job"]],**outcomes.get(c["condition"],{})}
    if "status" not in row:
        row["status"]="not-run-dataset-failure" if c["job"] in failures else "pending"
        row["audit_status"]="unavailable" if c["job"] in failures else "pending"
    rows.append(row)
    if row["family"]=="primary":groups[(row["case_id"],str(row["epsilon"]))].append(row)
keys=list(dict.fromkeys(k for r in rows for k in r))
tmp=root/"report/all_planned_outcomes.csv.writing"
with open(tmp,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
tmp.replace(root/"report/all_planned_outcomes.csv")
summary=[]
for (case,eps),rr in groups.items():
    attempted=[r for r in rr if r["condition"] in outcomes]
    audited=[r for r in attempted if r.get("audit_status")=="passed"]
    summary.append(dict(case=case,epsilon=eps,planned=len(rr),attempted=len(attempted),
        audited=len(audited),verified_feasible=sum(r.get("verified_feasible") is True for r in audited),
        verified_optimal=sum(r.get("verified_optimal") is True for r in audited),
        pending=sum(r["status"]=="pending" for r in rr),
        unavailable_or_failure=sum(r["status"]!="pending" and r not in audited for r in rr)))
with open(root/"report/primary_counts.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
progress=dict(updated=datetime.datetime.now().astimezone().isoformat(),
    planned_datasets=len(jobs),completed_datasets=completed,failed_datasets=len(failures),
    planned_searches=len(conditions),saved_searches=len(outcomes),
    audited_searches=sum(r.get("audit_status")=="passed" for r in outcomes.values()),
    search_status_counts=dict(collections.Counter(r.get("status") for r in outcomes.values())),
    reference_or_audit_unavailable=sum(r.get("audit_status") in {"reference-unavailable","audit-unavailable"} for r in outcomes.values()),
    denominator_note="All 10,400 retained searches remain in all_planned_outcomes.csv. Repeated searches are not independent datasets.")
atomic(root/"report/progress.json",progress)
print(json.dumps(progress))
