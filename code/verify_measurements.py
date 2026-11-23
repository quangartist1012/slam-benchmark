#!/usr/bin/env python3
"""Independently check recorded metric arithmetic, rigid alignment and plane fits."""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.spatial import cKDTree
from run_geometry import CORE_DIMS, CORE_PATCHES, dump, sha, now

O=Path(__file__).resolve().parent


def patch_points(p,key):
    z=(p[:,2]>=.35)&(p[:,2]<=2.15)
    if key.startswith('S'):
        station,j=map(int,key[1:].split('_'));x=[3,7,20][station-1];y=[0.,2.74][j]
        return p[z&(abs(p[:,0]-x)<=.5)&(abs(p[:,1]-y)<=.6)]
    if key=='branch_left':
        mask=(abs(p[:,0]-21.99)<.70)&(p[:,1]>-2.65)&(p[:,1]<-.55)
    elif key=='branch_right':
        mask=(abs(p[:,0]-27.69)<.70)&(p[:,1]>-2.65)&(p[:,1]<-.55)
    else:
        mask=(abs(p[:,1]+3.10)<.70)&(p[:,0]>22.59)&(p[:,0]<27.09)
    # Branch extraction uses strict z bounds, unlike station boxes.
    return p[mask&(p[:,2]>.35)&(p[:,2]<2.15)]


def main():
    from method_controls import control_checks
    controls=control_checks()
    report=json.loads((O/'measurements.json').read_text())
    summary=json.loads((O/'summary.json').read_text())
    inventory=json.loads((O/'source_inventory.json').read_text())
    rows=report['results']
    assert len(rows)==180 and len({(r['name'],r['bag_id']) for r in rows})==180
    assert len({r['bag_id'] for r in rows})==10 and len({r['name'] for r in rows})==18
    assert all(sha(O/f)==sig for f,sig in report['analysis_signature'].items())
    count=dict(maps=0,dimensions=0,planes=0,local_neighbourhoods=0,rigid_distances=0)
    import rosbag
    source_bags={}
    for bid in sorted({r['bag_id'] for r in rows}):
        source=O.parent.parent/'bags'/bid/'common.bag'
        with rosbag.Bag(str(source)) as bag:
            source_bags[bid]=dict(start=bag.get_start_time(),end=bag.get_end_time(),
                                 scans=bag.get_type_and_topic_info().topics['/velodyne_points'].message_count,
                                 sha256=sha(source))
    failures=[]
    by_source={(r['name'],r['bag_id']):r for r in inventory['runs']}
    max_plane_delta=0.
    for row in rows:
        source=by_source[row['name'],row['bag_id']]
        stat=Path(source['source_map']).stat()
        assert stat.st_size==source['source_stat']['size_bytes'] and stat.st_mtime_ns==source['source_stat']['mtime_ns']
        if row['status']!='measured':
            failures.append(dict(name=row['name'],bag_id=row['bag_id'],reason=row['reason']))
            continue
        count['maps']+=1
        d=O/'runs'/row['bag_id']/row['name']
        p=np.load(d/'cloud.npz')['cloud'].astype(float)
        q=np.load(d/'aligned.npz')['cloud'].astype(float)
        stamps=np.load(d/'own_poses.npz')['stamps']
        original=source_bags[row['bag_id']]
        assert stamps.min()>=original['start']-1 and stamps.max()<=original['end']+1, 'Output stamps disagree with source session'
        R=np.array(row['rotation']);t=np.array(row['translation'])
        assert abs(np.linalg.det(R)-1)<1e-7 and np.max(abs(R@R.T-np.eye(3)))<1e-7
        ix=np.linspace(0,len(p)-1,min(64,len(p)),dtype=int)
        assert np.max(abs(p[ix]@R.T+t-q[ix]))<1e-4
        assert np.max(abs(np.linalg.norm(np.diff(p[ix],axis=0),axis=1)-np.linalg.norm(np.diff(q[ix],axis=0),axis=1)))<1e-4
        count['rigid_distances']+=max(0,len(ix)-1)
        for dim in row['dimensions']:
            if dim['measured_m'] is None:
                assert dim['signed_error_m'] is None and dim['abs_error_m'] is None
                continue
            error=dim['measured_m']-dim['reference_m']
            assert abs(error-dim['signed_error_m'])<1e-12
            assert abs(abs(error)-dim['abs_error_m'])<1e-12
            assert abs(100*abs(error)/dim['reference_m']-dim['percent_error'])<1e-10
            count['dimensions']+=1
        vals=[d['signed_error_m'] for d in row['dimensions'] if d['key'] in CORE_DIMS]
        if all(v is not None for v in vals):
            assert len(vals)==2 and abs(np.mean(np.abs(vals))-row['core_dimension_mae_m'])<1e-12
            assert abs(np.sqrt(np.mean(np.square(vals)))-row['core_dimension_rmse_m'])<1e-12
        else:
            assert row['core_dimension_mae_m'] is None
        # Reconstruct in float64 to avoid a point changing ROI at float32 rounding.
        transformed=p@R.T+t
        for key,plane in row['surface_details'].items():
            if not plane['valid']:
                continue
            points=patch_points(transformed,key)
            expected=plane.get('n',plane.get('point_count'))
            assert len(points)==expected,(row['name'],row['bag_id'],key,len(points),expected)
            centred=points-points.mean(axis=0)
            # SVD is independent of the eigendecomposition used by measurement.
            _,_,vh=np.linalg.svd(centred,full_matrices=False)
            residual=centred@vh[-1]
            rmse=float(np.linalg.norm(residual)/np.sqrt(len(points)))
            delta=abs(rmse-plane['rmse_m'])
            max_plane_delta=max(max_plane_delta,delta)
            assert delta<1e-8
            count['planes']+=1
        if all(row['surface_details'][k]['valid'] for k in CORE_PATCHES):
            expected=np.sqrt(np.mean([row['surface_details'][k]['rmse_m']**2 for k in CORE_PATCHES]))
            assert abs(expected-row['core_surface_rmse_m'])<1e-12
        tree=cKDTree(transformed)
        for metric in row['local_quality']['metrics']:
            values=np.load(d/f'neighbourhoods_{int(metric["radius_m"]*100)}cm.npz')['values']
            if not len(values):
                assert metric['mpv_m'] is None and metric['mme_nats'] is None
                continue
            assert metric['valid_queries']==len(values)
            assert abs(np.mean(values[:,3])-metric['mpv_m'])<1e-12
            assert abs(np.mean(values[:,4])-metric['mme_nats'])<1e-12
            assert abs(100*len(values)/metric['query_count']-metric['support_percent'])<1e-10
            for k in np.linspace(0,len(values)-1,min(3,len(values)),dtype=int):
                v=values[k]
                indices=tree.query_ball_point(v[:3],metric['radius_m'])
                local=transformed[indices]
                centred=local-local.mean(axis=0)
                _,_,vh=np.linalg.svd(centred,full_matrices=False)
                mpv=np.percentile(abs(centred@vh[-1]),75)
                # slogdet on a covariance eigenvalue-regularized matrix.
                cov=np.cov(local,rowvar=False)
                eig,vec=np.linalg.eigh(cov)
                regularized=vec@np.diag(np.maximum(eig,1e-12))@vec.T
                sign,logdet=np.linalg.slogdet(2*np.pi*np.e*regularized)
                mme=.5*logdet
                assert len(local)==int(v[5]) and abs(mpv-v[3])<1e-8
                assert sign>0 and abs(mme-v[4])<1e-4
                count['local_neighbourhoods']+=1
        print('VERIFIED',row['name'],row['bag_id'],flush=True)
    assert count['dimensions']==summary['valid_dimensions']
    with (O/'bag_metrics.csv').open(encoding='utf-8-sig') as f:
        csv_rows=list(csv.DictReader(f))
    assert len(csv_rows)==180
    result=dict(status='numeric_verification_passed',created_utc=now(),checks=count,
                synthetic_software_controls=controls,
                original_sensor_bags=source_bags,
                invalid_maps_retained=failures,maximum_independent_plane_rmse_difference_m=max_plane_delta,
                plane_crosscheck='Independent SVD, all valid fixed patches, reconstruct original float64 rigid transform',
                local_crosscheck='Saved per-query values independently checked by SVD and covariance slogdet; up to three neighbourhoods per radius per map',
                limitations='Numeric verification does not certify physical reference dimensions or semantic ROI correspondence. Visual alignment review is separate.')
    dump(O/'verification.json',result)
    print('NUMERIC_VERIFICATION_PASSED',count,flush=True)


if __name__=='__main__':
    main()
