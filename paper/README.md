# Condensed paper: quantitative tables and merged figures

This folder accompanies the FoodMate benchmark: **16 upstream algorithms, 18 configurations,
10 physical recordings**. The paper manuscript itself is not included.

## Rebuild the paper assets

Tested on Ubuntu 20.04 with Python 3.8.10, NumPy 1.24.4 and Matplotlib 3.1.2.
No ROS, robot or hardware connection is required for these commands:

```bash
python3 rebuild_tables.py
python3 render_figures.py
```

Commands work from any current directory and write to `generated/`. They produce
the 18-row geometric and computational tables, bag-paired dimensional differences,
and five merged figures as PNG/PDF. Dimensions use measured lengths, without scaling.

## Evidence levels

- `data/primary_lengths.csv`: the same two directly measured lengths and their
  project-reference values for every configuration/recording; blanks are missing support.
  The rebuild recomputes per-map MAE/RMSE/MAPE and mean/sample SD across valid recordings.
- `data/geometry_summary.json`: retained per-map surface and neighbourhood measurements;
  plots use their original definitions and show individual recording values.
- `data/map_previews.npz`: **display subsamples**, up to 14,000 points per panel,
  from the already rigidly aligned maps of the fixed illustrative recording 164656.
  These are not full clouds and are not used to recalculate geometric accuracy.
- `data/resources_common_recording.json`: common-recording measurements, 1 s sampling.
  100% CPU means one logical core; RSS sums process memory. Both VoxelMap entries
  cover mapping-core resources only.
- `data/latency_summary.json`: three full timing replays per configuration, one
  recording, 10 s warmup excluded. New VoxelMap latency includes live preparation;
  DLIO ends at registered-cloud availability. Latency is not pure CPU time.
- `data/poses_and_stability_per_run.csv`: 180 native-output audits. Simulated ROS
  arrival intervals are not wall-clock latency. Completion screening is not accuracy.
- `data/software_versions.json`: upstream URLs/commits. Exact timing parameters,
  historical snapshots, patches and resource/pose/event observations are in the release below.

The archival [reproducibility release](https://github.com/quangartist1012/slam-benchmark/releases/tag/reproducibility-20261123)
contains the approximately 81 MB evidence ZIP, integrity checks, original measurement
scripts, retained configurations and instructions. Its portable analysis was tested
after extraction in a separate directory. The ten sensor bags and full output clouds
are not uploaded. Historical drivers require the original ROS project; a fresh-machine
rebuild and full raw-to-map replay are not claimed.

## Interpretation

Ten bags are separate hand-pushed recordings in one corridor, not ten environments
or repetitions of identical controlled motion. Different valid bag subsets must be
reported with n; pairwise contrasts use common recordings. Map dimensions have no
retained survey uncertainty budget. Surface compactness and trajectory disagreement
are internal consistency diagnostics, not absolute 3D or trajectory accuracy. No
overall best algorithm is inferred from one metric.

Interactive [180 maps](https://quangartist1012.github.io/slam-benchmark/)
and [trajectory consistency](https://quangartist1012.github.io/slam-benchmark/trajectories/)
remain available. Full parameter and sampling limitations are retained in their protocols
and in the release README. Third-party code retains its upstream notices; this folder
does not grant a new license to third-party material.
