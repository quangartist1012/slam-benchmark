#!/usr/bin/env python3
"""Resume all-ten-bag geometry analysis, after acquiring the shared timing lock.

No ROS nodes are started. Original bags, maps and previous pilot are read-only.
NumPy/SciPy/ROS imports happen only after the timing gate.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import time
import traceback

O = Path(__file__).resolve().parent
R = O.parent
PILOT = R / 'paper_closeout_20261120/geometry_3d'
CORE_DIMS = ['straight_length', 'branch_width']
CORE_PATCHES = ['S3_1', 'branch_left', 'branch_right', 'branch_back']
SCAN_NAMES = {'CT-ICP', 'VoxelMap-KISS', 'FAST-LIO2', 'Faster-LIO',
              'Point-LIO', 'LOG-LIO', 'VoxelMap-FAST', 'KISS-ICP', 'LIO-EKF'}


def now():
    return datetime.now(timezone.utc).isoformat()


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    os.replace(temp, path)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(4194304), b''):
            h.update(chunk)
    return h.hexdigest()


def state(status, **extra):
    dump(O / 'progress.json', dict(status=status, scope='criterion 1 only',
         updated_utc=now(), timing_campaign_active=False,
         cpu_heavy_analysis_active=status == 'running',
         word_write_owner='criterion 1 chat after coordination', **extra))


def inventory():
    rows = json.loads((R / 'paper_unified_20261120/unified_180.json').read_text())
    base = {x['name']: x for x in json.loads((R / 'performance_20261119/analysis.json').read_text())}
    common = {x['name']: x for x in json.loads((R / 'paper_unified_20261120/common_geometry.json').read_text())}
    assert len(rows) == 180 and len({(x['name'], x['bag_id']) for x in rows}) == 180
    for row in rows:
        name = row['name']
        directory = Path(row['output_dir'])
        candidates = [directory / 'output.bag', directory / 'slam_output.bag', directory / 'mapping/output.bag']
        found = [p for p in candidates if p.is_file()]
        assert len(found) == 1, (name, row['bag_id'], found)
        row['output_bag'] = str(found[0])
        row['pose_topic'] = base[name]['pose_topic']
        row['map_topic'] = common[name]['cloud_source']
        if name == 'KISS-ICP':
            row['map_topic'] = '/kiss/frame'
        if name == 'LIO-EKF':
            row['map_topic'] = '/frame'
        native = directory / 'finalCloud.pcd' if name.startswith('LeGO-') else (
            directory / 'map/GlobalMap.pcd' if name == 'LIO-SAM' else None)
        assert native is None or native.is_file(), (name, str(native))
        row['native_map'] = str(native) if native else None
        row['source_map'] = str(native or found[0])
        st = Path(row['source_map']).stat()
        row['source_stat'] = dict(size_bytes=st.st_size, mtime_ns=st.st_mtime_ns)
    bags = sorted({x['bag_id'] for x in rows})
    assert len(bags) == 10 and len(base) == 18
    dump(O / 'source_inventory.json', dict(created_utc=now(),
         manifest_sha256=sha(R / 'paper_unified_20261120/unified_180.json'),
         counts=dict(bags=10, configurations=18, pairs=180), bag_ids=bags,
         configurations=list(base), runs=rows))
    return rows


def gate(wait):
    lock = (R / '.benchmark_timing.lock').open('a+')
    last = None
    while True:
        progress = R / 'criterion_3_20261123/progress.json'
        other = json.loads(progress.read_text()) if progress.is_file() else {}
        busy = other.get('timing_campaign_active') or other.get('cpu_heavy_analysis_active')
        acquired = False
        if not busy:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = True
            except BlockingIOError:
                pass
        if acquired:
            # Second read after locking prevents a startup race with the timer.
            other = json.loads(progress.read_text()) if progress.is_file() else {}
            if not (other.get('timing_campaign_active') or other.get('cpu_heavy_analysis_active')):
                return lock
            fcntl.flock(lock, fcntl.LOCK_UN)
        current = (other.get('completed'), other.get('current'), other.get('repeat'))
        if current != last:
            state('waiting_for_timing', other_progress=other,
                  coordination='No NumPy, ROS bag reads, extraction or plotting while timing campaign owns the machine.')
            print('WAITING_FOR_TIMING', current, flush=True)
            last = current
        if not wait:
            raise SystemExit('Criterion 3 is active; re-run with --wait-for-timing or after it finishes.')
        time.sleep(10)


def read_pose(row):
    import numpy as np
    import rosbag
    stamps, pos, quat, frames = [], [], [], set()
    with rosbag.Bag(row['output_bag']) as bag:
        for _, m, _ in bag.read_messages(topics=[row['pose_topic']]):
            p, q = m.pose.pose.position, m.pose.pose.orientation
            v = np.array([p.x, p.y, p.z, q.x, q.y, q.z, q.w])
            if m.header.stamp.to_sec() <= 0 or not np.isfinite(v).all() or np.linalg.norm(v[3:]) < 1e-9:
                continue
            stamps.append(m.header.stamp.to_sec())
            pos.append(v[:3])
            quat.append(v[3:] / np.linalg.norm(v[3:]))
            frames.add(m.header.frame_id)
    assert len(stamps) >= 20, 'insufficient valid own odometry'
    tt = np.array(stamps)
    ix = np.argsort(tt, kind='stable')
    return tt[ix], np.array(pos)[ix], np.array(quat)[ix], sorted(frames)


def extract(row, directory, pose):
    import numpy as np
    import rosbag
    from scipy.spatial.transform import Rotation
    import geometry_core as g
    file = directory / 'cloud.npz'
    meta = directory / 'extraction.json'
    signature = dict(source=row['source_map'], **row['source_stat'],
                     voxel_m=.05, revision=1, geometry_core_sha256=sha(O / 'geometry_core.py'))
    if file.exists() and meta.exists():
        old = json.loads(meta.read_text())
        if old.get('signature') == signature:
            return np.load(file)['cloud'].astype(float), old
    name = row['name']
    topic = row['map_topic']
    frames, stamps, total, used, max_dt = set(), [], 0, 0, 0.
    if row['native_map']:
        points = g.pcd(Path(row['native_map']))
        total, used, method = len(points), 1, 'native final map'
        frames = set(pose[3])
    else:
        parts = []
        with rosbag.Bag(row['output_bag']) as bag:
            info = bag.get_type_and_topic_info().topics
            assert topic in info, ('map topic absent', topic, list(info))
            for i, (_, m, _) in enumerate(bag.read_messages(topics=[topic])):
                if name not in SCAN_NAMES and i != info[topic].message_count - 1:
                    continue
                a = g.xyz(m)
                if not len(a):
                    continue
                ts = m.header.stamp.to_sec()
                frames.add(m.header.frame_id)
                stamps.append(ts)
                total += len(a)
                used += 1
                if name == 'LIO-EKF':
                    tt, pos, quat, _ = pose
                    j = int(np.argmin(abs(tt-ts)))
                    dt = float(abs(tt[j]-ts))
                    max_dt = max(max_dt, dt)
                    assert dt < .001 and m.header.frame_id == 'pointcloud_frame', ('body scan pose match', dt, m.header.frame_id)
                    a = Rotation.from_quat(quat[j]).apply(a) + pos[j]
                if name == 'KISS-ICP':
                    assert m.header.frame_id == 'odom', ('unexpected KISS scan frame', m.header.frame_id)
                parts.append(g.voxel(a))
                if len(parts) >= 50:
                    parts = [g.voxel(np.vstack(parts))]
        assert parts, 'no nonempty map messages'
        points = np.vstack(parts)
        method = 'all registered scans' if name in SCAN_NAMES else 'last full published map'
    if name.startswith('LeGO-'):
        points = points[:, [2, 0, 1]]
    points = g.voxel(points).astype(np.float32)
    assert len(points) >= 100, 'insufficient extracted map points'
    np.savez_compressed(file, cloud=points)
    value = dict(signature=signature, source_sha256=sha(Path(row['source_map'])),
                 cloud_npz_sha256=sha(file), extraction=method, topic=None if row['native_map'] else topic,
                 source_frames=sorted(frames), point_count=len(points),
                 point_count_before_common_voxel=total, messages_used=used,
                 axis_permutation=[2, 0, 1] if name.startswith('LeGO-') else [0, 1, 2],
                 source_stamp_span=[min(stamps), max(stamps)] if stamps else None,
                 body_pose_max_difference_s=max_dt if name == 'LIO-EKF' else None,
                 quantiles_xyz=np.quantile(points, [.01, .5, .99], axis=0).tolist(),
                 source_native_filtering_different=True)
    dump(meta, value)
    return points.astype(float), value


def patches(metrics):
    result = {}
    for i, station in enumerate(metrics.get('stations', []), 1):
        for j, wall in enumerate(station['walls']):
            result[f'S{i}_{j}'] = wall
    for key, value in metrics.get('branch_planes', {}).items():
        result['branch_' + key] = value
    return result


def aggregates(metrics):
    import numpy as np
    dims = {x['key']: x for x in metrics['dimensions']}
    valid = all(dims[k]['measured_m'] is not None for k in CORE_DIMS)
    errors = np.array([dims[k]['signed_error_m'] for k in CORE_DIMS]) if valid else None
    pp = patches(metrics)
    surface_valid = all(pp.get(k, {}).get('valid') for k in CORE_PATCHES)
    return dict(core_dimension_keys=CORE_DIMS, core_surface_keys=CORE_PATCHES,
                core_dimensions_valid=bool(valid),
                core_dimension_mae_m=float(np.mean(abs(errors))) if valid else None,
                core_dimension_rmse_m=float(np.sqrt(np.mean(errors**2))) if valid else None,
                core_dimension_max_m=float(max(abs(errors))) if valid else None,
                core_surfaces_valid=bool(surface_valid),
                core_surface_rmse_m=float(np.sqrt(np.mean([pp[k]['rmse_m']**2 for k in CORE_PATCHES]))) if surface_valid else None,
                surface_details=pp)


def local_quality(p, directory):
    import numpy as np
    from scipy.spatial import cKDTree
    mask = (p[:, 0] >= .5) & (p[:, 0] <= 21) & (p[:, 1] >= -.4) & (p[:, 1] <= 3.14) & (p[:, 2] >= -.3) & (p[:, 2] <= 2.7)
    roi = p[mask]
    _, ids = np.unique(np.floor(roi / .2).astype(int), axis=0, return_index=True)
    query = roi[ids]
    if len(query) > 4000:
        query = query[np.linspace(0, len(query)-1, 4000, dtype=int)]
    tree = cKDTree(p)
    results = []
    for radius in [.3, .6]:
        values = []
        for q, indices in zip(query, tree.query_ball_point(query, radius)):
            if len(indices) < 10:
                continue
            points = p[indices]
            v = points - points.mean(axis=0)
            cov = v.T @ v / (len(v)-1)
            eig, vec = np.linalg.eigh(cov)
            mpv = float(np.quantile(abs(v @ vec[:, 0]), .75))
            mme = .5 * float(3*np.log(2*np.pi*np.e) + np.log(np.maximum(eig, 1e-12)).sum())
            values.append([*q, mpv, mme, len(v)])
        array = np.asarray(values, dtype=float).reshape(-1, 6)
        np.savez_compressed(directory / f'neighbourhoods_{int(100*radius)}cm.npz', values=array)
        results.append(dict(radius_m=radius, query_count=len(query), valid_queries=len(array),
             support_percent=100*len(array)/len(query) if len(query) else None,
             excluded_insufficient_neighbors=len(query)-len(array),
             mpv_m=float(array[:, 3].mean()) if len(array) else None,
             mpv_p95_m=float(np.quantile(array[:, 3], .95)) if len(array) else None,
             mme_nats=float(array[:, 4].mean()) if len(array) else None,
             neighbor_median=float(np.median(array[:, 5])) if len(array) else None))
    # A fixed lattice is a support diagnostic, not dense-reference completeness.
    # Same corridor wall locations for all maps, unlike map-dependent queries.
    lattice = np.array([[x, y, z] for x in np.arange(.5, 21.01, .5)
                        for y in [0., 2.74] for z in np.arange(.5, 2.01, .3)])
    distance, _ = tree.query(lattice, k=1)
    support = dict(definition='Fraction of fixed corridor-wall lattice sites with any retained map point <=0.20 m; diagnostic, not surveyed completeness',
                   sites=len(lattice), within_20cm=int(np.sum(distance <= .20)),
                   percent_within_20cm=float(100*np.mean(distance <= .20)),
                   median_nearest_m=float(np.median(distance)),
                   p95_nearest_m=float(np.quantile(distance, .95)))
    return dict(roi_point_count=len(roi), metrics=results, fixed_wall_lattice_support=support)


def analyse(row, signature):
    import numpy as np
    import geometry_core as g
    d = O / 'runs' / row['bag_id'] / row['name']
    d.mkdir(parents=True, exist_ok=True)
    resultfile = d / 'measurement.json'
    bag_stat = Path(row['output_bag']).stat()
    identity = dict(source_map=row['source_map'], source_stat=row['source_stat'],
                    output_bag=row['output_bag'], output_bag_size_bytes=bag_stat.st_size,
                    output_bag_mtime_ns=bag_stat.st_mtime_ns)
    if resultfile.exists():
        old = json.loads(resultfile.read_text())
        if old.get('analysis_signature') == signature and old.get('source_identity') == identity:
            return old
    started = time.monotonic()
    result = dict(name=row['name'], bag_id=row['bag_id'], source_duration_s=row['source_duration_s'],
                  analysis_signature=signature, source_identity=identity, created_utc=now(), status='failed')
    try:
        pose = read_pose(row)
        p, extraction = extract(row, d, pose)
        result['extraction'] = extraction
        own = pose[1]
        if row['name'].startswith('LeGO-'):
            own = own[:, [2, 0, 1]]
        np.savez_compressed(d / 'own_poses.npz', stamps=pose[0], positions=own, quaternions=pose[2])
        s, n, good = g.normals(p)
        R0, z0, seed = g.floor_seed(s, n, good, own)
        s = s @ R0.T
        s[:, 2] -= z0
        n = n @ R0.T
        own = own @ R0.T
        own[:, 2] -= z0
        R1, t1, registration = g.wall_pair(s, n, good, own)
        rotation = R1 @ R0
        translation = t1 + R1 @ np.array([0., 0., -z0])
        assert np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-7) and np.linalg.det(rotation) > .99999
        aligned = p @ rotation.T + translation
        route = own @ R1.T + t1
        np.savez_compressed(d / 'aligned.npz', cloud=aligned.astype(np.float32), route=route.astype(np.float32))
        measurement, _ = g.assess(aligned)
        # Pilot combined-check aggregates mixed lengths and one registered
        # coordinate. Keep those only as diagnostics, never as primary MAE.
        for key in ['dimension_mae_m', 'dimension_rmse_m', 'dimension_max_abs_m']:
            measurement['diagnostic_all_checks_' + key] = measurement.pop(key)
        result.update(status='measured', rotation=rotation.tolist(), translation=translation.tolist(),
                      determinant=float(np.linalg.det(rotation)), scale=1., floor_seed=seed,
                      registration=registration, point_count=len(p), **measurement)
        result.update(aggregates(measurement))
        result['local_quality'] = local_quality(aligned, d)
        sensitivity = []
        for key, cm, bm, z in [('margin_minus_10cm', .5, .6, 1.2), ('margin_plus_10cm', .7, .8, 1.2),
                              ('probe_1m', .6, .7, 1.), ('probe_1p4m', .6, .7, 1.4)]:
            mm, _ = g.assess(aligned, corridor_margin=cm, branch_margin=bm, probe_z=z)
            aa = aggregates(mm)
            sensitivity.append(dict(variant=key, corridor_margin_m=cm, branch_margin_m=bm, probe_z_m=z,
                  dimensions=mm['dimensions'], core_dimension_mae_m=aa['core_dimension_mae_m'],
                  core_surface_rmse_m=aa['core_surface_rmse_m'], valid_dimensions=mm['valid_dimensions']))
        result['definition_sensitivity'] = sensitivity
        result['review_status'] = 'pending_semantic_alignment_review'
    except Exception as exc:
        result.update(status='not_measurable_automatically', reason=f'{type(exc).__name__}: {exc}',
                      failure_traceback=traceback.format_exc())
    result['analysis_wall_time_s'] = time.monotonic() - started
    dump(resultfile, result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inventory-only', action='store_true')
    parser.add_argument('--wait-for-timing', action='store_true')
    parser.add_argument('--names', nargs='+')
    parser.add_argument('--bags', nargs='+')
    parser.add_argument('--postprocess', action='store_true', help='After all 180 maps, summarize, render, independently verify and package the local supplement.')
    a = parser.parse_args()
    rows = inventory()
    if a.inventory_only:
        print('SOURCE_INVENTORY_OK: 180 real map outputs, ten common source bags, eighteen configurations')
        return
    if a.names:
        assert set(a.names) <= {r['name'] for r in rows}
        rows = [r for r in rows if r['name'] in a.names]
    if a.bags:
        assert set(a.bags) <= {r['bag_id'] for r in rows}
        rows = [r for r in rows if r['bag_id'] in a.bags]
    # First verify LIO-SAM across ten bags; then extend exactly the same rules.
    rows.sort(key=lambda r: (r['name'] != 'LIO-SAM', r['name'], r['bag_id']))
    lock = gate(a.wait_for_timing)
    # Keep BLAS from spinning all logical cores during map postprocessing.
    for key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ[key] = '1'
    signature = {f: sha(O / f) for f in ['run_geometry.py', 'geometry_core.py', 'protocol.json']}
    signature['../paper_unified_20261120/unified_180.json'] = sha(R / 'paper_unified_20261120/unified_180.json')
    started = now()
    manifest = dict(started_utc=started, analysis_signature=signature, expected_runs=len(rows), results=[])
    ending_status = 'interrupted'
    try:
        for i, row in enumerate(rows):
            state('running', current=row['name'], bag_id=row['bag_id'], completed_maps=i, total_maps=len(rows))
            measured = analyse(row, signature)
            manifest['results'].append(measured)
            dump(O / 'measurements.json', manifest)
            print(i+1, '/', len(rows), row['name'], row['bag_id'], measured['status'],
                  measured.get('valid_dimensions'), measured.get('core_dimension_mae_m'), flush=True)
        manifest.update(status='analysis_complete_review_pending', completed_utc=now())
        dump(O / 'measurements.json', manifest)
        if a.postprocess:
            assert len(rows) == 180, 'Postprocessing requires the complete 180-map campaign.'
            import subprocess
            for script in ['summarize.py', 'render_figures.py', 'verify_measurements.py', 'build_supplement.py']:
                state('running', current_stage=script, completed_maps=180, total_maps=180)
                subprocess.run(['python3', '-u', str(O / script)], check=True, cwd=str(O))
        ending_status = 'analysis_complete_review_pending'
        state('analysis_complete_review_pending', completed_maps=len(rows), total_maps=len(rows))
    except Exception as exc:
        ending_status = 'failed'
        dump(O / 'workflow_error.json', dict(stage='geometry_or_postprocess', error=str(exc), traceback=traceback.format_exc(), utc=now()))
        raise
    finally:
        state(ending_status,
              completed_maps=len(manifest['results']), total_maps=len(rows))
        fcntl.flock(lock, fcntl.LOCK_UN)


if __name__ == '__main__':
    main()
