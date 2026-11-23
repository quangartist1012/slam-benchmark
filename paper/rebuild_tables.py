#!/usr/bin/env python3
"""Recompute paper tables from retained measurements, using portable relative paths."""
from pathlib import Path
import csv,hashlib,json
import numpy as np

O=Path(__file__).resolve().parent;D=O/'data';OUT=O/'generated'
def read(p):return json.loads(p.read_text())
def csvread(p):
    with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def save(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    OUT.mkdir(exist_ok=True)
    for line in (O/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ',1);assert hashlib.sha256((O/name).read_bytes()).hexdigest()==digest,name
    geo=read(D/'geometry_summary.json');names=geo['configurations'];ss={x['name']:x for x in geo['summaries']};raw=csvread(D/'primary_lengths.csv');permap=[]
    for n in names:
        for b in geo['bag_ids']:
            rr=[x for x in raw if x['name']==n and x['bag_id']==b];assert len(rr)==2
            if any(x['measured_m']=='' for x in rr):continue
            errors=np.array([float(x['measured_m'])-float(x['reference_m']) for x in rr]);mape=np.mean([100*abs(float(x['measured_m'])-float(x['reference_m']))/float(x['reference_m']) for x in rr]);mae=float(abs(errors).mean());rmse=float(np.sqrt(np.mean(errors**2)))
            expected=next(x for x in geo['bag_metrics'] if x['name']==n and x['bag_id']==b);assert abs(mae-expected['core_dimension_mae_m'])<1e-10;assert abs(rmse-expected['core_dimension_rmse_m'])<1e-10
            permap.append(dict(name=n,bag_id=b,mae_m=mae,rmse_m=rmse,mape_percent=float(mape)))
    table=[]
    for n in names:
        rr=[x for x in permap if x['name']==n];z=ss[n]['metrics'];a=np.array([r['mae_m'] for r in rr]);row=dict(name=n,n=len(rr),mae_cm=float(a.mean()*100),mae_sd_cm=float(a.std(ddof=1)*100),rmse_cm=float(np.mean([r['rmse_m'] for r in rr])*100),primary_mape_percent=float(np.mean([r['mape_percent'] for r in rr])),mpv_cm=100*z['mpv_m']['mean'],surface_rmse_cm=100*z['core_surface_rmse_m']['mean']);assert abs(row['mae_cm']-100*z['core_dimension_mae_m']['mean'])<1e-8;table.append(row)
    save('condensed_geometry.csv',table);save('geometry_per_map.csv',permap)
    paired=[]
    for a,b in [('LIO-EKF','DLO'),('LeGO-LiDAR','LeGO-IMU'),('KISS-ICP','LIO-SAM')]:
        for bag in geo['bag_ids']:
            x=next((x for x in permap if x['name']==a and x['bag_id']==bag),None);y=next((x for x in permap if x['name']==b and x['bag_id']==bag),None)
            if x is not None and y is not None:paired.append(dict(first=a,second=b,bag_id=bag,first_mae_cm=100*x['mae_m'],second_mae_cm=100*y['mae_m'],first_minus_second_cm=100*(x['mae_m']-y['mae_m'])))
    save('paired_dimension_differences.csv',paired)
    cpu={x['name']:x for x in read(D/'resources_common_recording.json')};lat={x['name']:x for x in read(D/'latency_summary.json')['results']};resource=[]
    for n in names:resource.append(dict(name=n,cpu_mean_percent_one_core=cpu[n]['cpu_mean_percent_one_core'],rss_mean_mib=cpu[n]['rss_mean_mib'],mean_latency_median_ms=lat[n]['mean_ms']['median'],p95_latency_median_ms=lat[n]['p95_ms']['median'],resource_scope='mapping core only' if n.startswith('VoxelMap-') else 'selected launch descendants',timing_endpoint=lat[n]['endpoint']))
    save('condensed_computation.csv',resource)
    report=dict(passed=True,configurations=len(names),valid_geometry_maps=len(permap),paired_values=len(paired),same_primary_lengths=True,scale_change=False)
    (OUT/'tables_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
