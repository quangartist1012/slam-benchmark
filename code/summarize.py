#!/usr/bin/env python3
"""Bag-level and paired statistics. Missing measurements never become zeros."""
import csv
import itertools
import json
from pathlib import Path
import statistics

from run_geometry import CORE_DIMS, CORE_PATCHES, dump, now

O = Path(__file__).resolve().parent


def stats(values):
    values = [float(v) for v in values if v is not None]
    return dict(n=len(values), mean=statistics.mean(values) if values else None,
                sd=statistics.stdev(values) if len(values) > 1 else None,
                median=statistics.median(values) if values else None,
                min=min(values) if values else None, max=max(values) if values else None)


def csv_write(path, rows):
    if not rows:
        return
    with path.open('w', newline='', encoding='utf-8-sig') as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    manifest = json.loads((O / 'measurements.json').read_text())
    inventory = json.loads((O / 'source_inventory.json').read_text())
    rows = manifest['results']
    assert len(rows) == 180 and len({(r['name'], r['bag_id']) for r in rows}) == 180
    assert manifest['status'] == 'analysis_complete_review_pending'
    names, bags = inventory['configurations'], inventory['bag_ids']
    by_key = {(r['name'], r['bag_id']): r for r in rows}
    common_dimension_bags = [b for b in bags if all(by_key[n, b].get('core_dimensions_valid') for n in names)]
    common_surface_bags = [b for b in bags if all(by_key[n, b].get('core_surfaces_valid') for n in names)]
    dimensions, surfaces, bag_rows, local, sensitivities = [], [], [], [], []
    for r in rows:
        name, bag = r['name'], r['bag_id']
        for dd in r.get('dimensions', []):
            dimensions.append(dict(name=name, bag_id=bag,
                 quantity_type='registered_coordinate_diagnostic' if dd['key']=='branch_depth' else 'direct_projected_length',
                 included_in_primary_aggregate=dd['key'] in CORE_DIMS, **dd))
        for key, patch in r.get('surface_details', {}).items():
            surfaces.append(dict(name=name, bag_id=bag, surface=key, valid=patch['valid'],
                                 points=patch.get('n', patch.get('point_count')),
                                 support_bins=patch.get('support_bins', patch.get('longitudinal_20cm_bins')),
                                 rmse_m=patch.get('rmse_m') if patch['valid'] else None))
        for lm in r.get('local_quality', {}).get('metrics', []):
            local.append(dict(name=name, bag_id=bag, **lm))
        sensitivity_values = [r.get('core_dimension_mae_m')]
        surface_sensitivity_values = [r.get('core_surface_rmse_m')]
        for variant in r.get('definition_sensitivity', []):
            sensitivities.append(dict(name=name, bag_id=bag, variant=variant['variant'],
                 valid_dimensions=variant['valid_dimensions'], core_dimension_mae_m=variant['core_dimension_mae_m'],
                 core_surface_rmse_m=variant['core_surface_rmse_m']))
            sensitivity_values.append(variant['core_dimension_mae_m'])
            surface_sensitivity_values.append(variant['core_surface_rmse_m'])
        ss, tt = stats(sensitivity_values), stats(surface_sensitivity_values)
        lm = next((v for v in r.get('local_quality', {}).get('metrics', []) if v['radius_m'] == .6), {})
        lattice = r.get('local_quality', {}).get('fixed_wall_lattice_support', {})
        bag_rows.append(dict(name=name, bag_id=bag, status=r['status'],
             core_dimension_mae_m=r.get('core_dimension_mae_m'),
             core_dimension_rmse_m=r.get('core_dimension_rmse_m'),
             core_dimension_max_m=r.get('core_dimension_max_m'),
             core_surface_rmse_m=r.get('core_surface_rmse_m'),
             valid_dimensions=r.get('valid_dimensions', 0), point_count=r.get('point_count'),
             mpv_m=lm.get('mpv_m'), mme_nats=lm.get('mme_nats'), neighbourhood_support_percent=lm.get('support_percent'),
             fixed_wall_support_percent=lattice.get('percent_within_20cm'),
             sensitivity_dimension_mae_range_m=ss['max']-ss['min'] if ss['n'] == 5 else None,
             sensitivity_surface_rmse_range_m=tt['max']-tt['min'] if tt['n'] == 5 else None,
             sensitivity_valid_variants=ss['n'], failure_reason=r.get('reason'),
             review_status=r.get('review_status', 'not_available')))
    metrics = ['core_dimension_mae_m', 'core_dimension_rmse_m', 'core_dimension_max_m',
               'core_surface_rmse_m', 'mpv_m', 'mme_nats', 'neighbourhood_support_percent',
               'fixed_wall_support_percent', 'sensitivity_dimension_mae_range_m',
               'sensitivity_surface_rmse_range_m']
    summaries = []
    for name in names:
        rr = [r for r in bag_rows if r['name'] == name]
        summaries.append(dict(name=name, attempted_bags=10,
              aligned_maps=sum(r['status'] == 'measured' for r in rr),
              metrics={key: stats(r[key] for r in rr) for key in metrics},
              balanced_dimension_mae_m=stats(by_key[name, b]['core_dimension_mae_m'] for b in common_dimension_bags),
              balanced_surface_rmse_m=stats(by_key[name, b]['core_surface_rmse_m'] for b in common_surface_bags),
              per_dimension={key: stats(r['signed_error_m'] for r in dimensions if r['name'] == name and r['key'] == key)
                             for key in ['width_3', 'width_7', 'width_20', *CORE_DIMS, 'branch_depth']}))
    paired = []
    indexed = {(r['name'], r['bag_id']): r for r in bag_rows}
    for a, b in itertools.combinations(names, 2):
        for metric in ['core_dimension_mae_m', 'core_surface_rmse_m', 'mpv_m']:
            valid = [bag for bag in bags if indexed[a, bag][metric] is not None and indexed[b, bag][metric] is not None]
            delta = [indexed[a, bag][metric] - indexed[b, bag][metric] for bag in valid]
            paired.append(dict(first=a, second=b, metric=metric, paired_bags=';'.join(valid), **stats(delta)))
    for filename, values in [('dimension_measurements.csv', dimensions), ('surface_measurements.csv', surfaces),
                              ('bag_metrics.csv', bag_rows), ('local_metrics.csv', local),
                              ('definition_sensitivity.csv', sensitivities), ('paired_differences.csv', paired)]:
        csv_write(O / filename, values)
    dimension_summary=[]
    for name in names:
        for key in ['width_3','width_7','width_20',*CORE_DIMS,'branch_depth']:
            valid=[r for r in dimensions if r['name']==name and r['key']==key and r['measured_m'] is not None]
            errors=[r['signed_error_m'] for r in valid]
            measured=stats(r['measured_m'] for r in valid)
            absolute=stats(abs(v) for v in errors)
            relative=stats(r['percent_error'] for r in valid)
            dimension_summary.append(dict(name=name,key=key,
                quantity_type='registered_coordinate_diagnostic' if key=='branch_depth' else 'direct_projected_length',
                valid_bags=len(valid),attempted_bags=10,reference_m=valid[0]['reference_m'] if valid else None,
                measured_mean_m=measured['mean'],measured_sd_m=measured['sd'],
                signed_error_mean_m=statistics.mean(errors) if errors else None,
                mae_m=absolute['mean'],rmse_m=(statistics.mean(v*v for v in errors)**.5) if errors else None,
                max_abs_error_m=absolute['max'],mean_abs_error_percent=relative['mean'],max_abs_error_percent=relative['max']))
    csv_write(O / 'dimension_summary.csv',dimension_summary)
    summary_flat = []
    for row in summaries:
        flat = dict(name=row['name'], attempted_bags=10, aligned_maps=row['aligned_maps'])
        for metric, st in row['metrics'].items():
            for key, value in st.items():
                flat[metric + '_' + key] = value
        summary_flat.append(flat)
    csv_write(O / 'configuration_summary.csv', summary_flat)
    final = dict(generated_utc=now(), maps_attempted=180, geometric_checks_attempted=1080,
                 direct_dimension_checks_attempted=900, registered_coordinate_checks_attempted=180,
                 maps_automatically_aligned=sum(r['status'] == 'measured' for r in rows),
                 valid_dimensions=sum(d['measured_m'] is not None for d in dimensions),
                 valid_direct_dimensions=sum(d['measured_m'] is not None and d['key'] != 'branch_depth' for d in dimensions),
                 valid_registered_coordinates=sum(d['measured_m'] is not None and d['key'] == 'branch_depth' for d in dimensions),
                 core_dimensions=CORE_DIMS, core_surfaces=CORE_PATCHES,
                 common_dimension_bag_ids=common_dimension_bags, common_surface_bag_ids=common_surface_bags,
                 bag_ids=bags, configurations=names, summaries=summaries, bag_metrics=bag_rows,dimension_summary=dimension_summary,
                 missing_maps=[{k:r.get(k) for k in ['name', 'bag_id', 'reason']} for r in rows if r['status'] != 'measured'],
                 interpretation='Descriptive results on ten recorded sessions in one environment. Do not rank mixed-denominator means; use paired common bags or balanced subsets. Reference uncertainty is undocumented; local quality metrics do not establish global correctness.',
                 semantic_alignment_review_complete=False)
    dump(O / 'summary.json', final)
    print('SUMMARIZED', final['maps_automatically_aligned'], '/180 maps;', final['valid_dimensions'], '/1080 geometric checks;',
          'all-18 common dimension bags:', len(common_dimension_bags), flush=True)


if __name__ == '__main__':
    main()
