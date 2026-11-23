FoodMate: Supplementary Map Geometry Evidence

Scope: 16 algorithms, 18 configurations, 10 real sensor recordings (13 November 2026), 180 saved map outputs.

Open index.html locally, keeping all sibling directories. For GitHub Pages, publish the main branch root; .nojekyll disables Jekyll processing. All paths are relative, including per-map dynamically loaded JavaScript files. No build step or external JavaScript dependency is required.

Map previews: deterministic samples of at most 10,000 points per map, not the complete clouds used for the calculations. Full sensor bags and full point clouds are not included. The download tables contain measurement results, support counts, missing cases, parameter sensitivity and paired comparisons. See downloads/protocol.json and downloads/provenance.json for definitions and provenance.

Limitations: no independent dense 3D ground truth, and the project reference sketch has no documented survey uncertainty. Dimension agreement, plane RMSE and local MPV/MME are reported with their different meanings. Missing registration/measurement is not labelled as a SLAM failure.

Computational evidence and condensed paper assets

The paper/ directory provides portable scripts, derived inputs and main-table generation for all 18 configurations. See paper/README.md. The versioned evidence ZIP (original resource CSV, native poses, timing events, configuration snapshots and checksums) is available at https://github.com/quangartist1012/slam-benchmark/releases/tag/reproducibility-20261123 . The manuscript draft and full raw sensor bags are not published.
