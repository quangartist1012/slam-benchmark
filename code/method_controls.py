#!/usr/bin/env python3
"""Synthetic software controls ONLY. Never include these in physical-trial counts."""
from pathlib import Path
import numpy as np
from geometry_core import assess
from run_geometry import aggregates, dump, now

O=Path(__file__).resolve().parent


def fixture(width=2.74,length=21.99,branch=5.70,depth=3.10,ty=0.,sparse_left=False,ghost=False):
    rng=np.random.RandomState(20261123)
    parts=[]
    for y in [ty,ty+width]:
        parts.append(np.column_stack([rng.uniform(0,length,7000),np.full(7000,y),rng.uniform(.25,2.3,7000)]))
    for x in [length,length+branch]:
        n=5 if sparse_left and x==length else 1500
        parts.append(np.column_stack([np.full(n,x),rng.uniform(-depth+ty,ty,n),rng.uniform(.25,2.3,n)]))
    parts.append(np.column_stack([rng.uniform(length,length+branch,3000),np.full(3000,ty-depth),rng.uniform(.25,2.3,3000)]))
    if ghost:
        extra=parts[1].copy();extra[:,1]+=.15;parts.append(extra)
    return np.vstack(parts)


def control_checks():
    exact,_=assess(fixture())
    direct={d['key']:d for d in exact['dimensions']}
    assert all(abs(direct[k]['signed_error_m'])<1e-10 for k in ['width_3','width_7','width_20','straight_length','branch_width'])
    deformed,_=assess(fixture(width=3.014,length=22.2099,branch=5.985,ty=-.137))
    dd={d['key']:d for d in deformed['dimensions']}
    assert abs(dd['width_3']['measured_m']-3.014)<1e-9
    assert abs(dd['straight_length']['measured_m']-22.2099)<1e-9
    assert abs(dd['branch_width']['measured_m']-5.985)<1e-9
    assert abs(dd['branch_depth']['measured_m']-3.237)<1e-9
    agg=aggregates(deformed)
    expected=(abs(22.2099-21.99)+abs(5.985-5.7))/2
    assert abs(agg['core_dimension_mae_m']-expected)<1e-9
    assert agg['core_dimension_keys']==['straight_length','branch_width']
    missing,_=assess(fixture(sparse_left=True))
    mm={d['key']:d for d in missing['dimensions']}
    assert mm['straight_length']['measured_m'] is None and mm['branch_width']['measured_m'] is None
    assert aggregates(missing)['core_dimension_mae_m'] is None
    ghost,_=assess(fixture(ghost=True))
    assert ghost['stations'][0]['walls'][1]['rmse_m']>.05
    assert exact['stations'][0]['walls'][1]['rmse_m']<1e-10
    checks=[dict(name='known reference lengths',passed=True),
            dict(name='changed dimensions remain changed; reference not forced by assessment',passed=True,
                 synthetic_expected_mae_m=expected),
            dict(name='reference-centred back coordinate is excluded from primary lengths',passed=True,
                 physical_depth_m=3.10,registered_coordinate_m=3.237),
            dict(name='insufficient wall support yields missing, never zero error',passed=True),
            dict(name='all ROI returns retained; ghost surface increases dispersion',passed=True)]
    report=dict(status='passed',created_utc=now(),data_type='Synthetic deterministic software controls; not robot measurements, not additional trials',checks=checks)
    dump(O/'method_control_checks.json',report)
    return report


if __name__=='__main__':
    print(control_checks())
