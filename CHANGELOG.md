# Changelog

All notable changes to MuRAT are documented here.
Dates follow ISO 8601 (YYYY-MM-DD). Versions follow [Semantic Versioning](https://semver.org).

---

## [Unreleased]

### Added
- First-difference (gradient) smoothing for the `'Tikhonov'` inversions of
  Qc and Q. `smoothingValueQc` and `smoothingValueQ` were not used by the
  Tikhonov branch; they now set the weight `s` of a first-difference
  regulariser on the model grid, so the objective is
  `norm(A*m-d)^2 + (dampingValue*norm(A))^2*norm(m)^2 +
  (s*norm(A))^2*norm(D*m)^2` (m >= 0).
  A uniform model is not penalised by `D`, so the solution is no longer
  pulled towards zero. Empty or 0 keeps the previous behaviour; if smoothing
  is set and `dampingValue*` is empty, damping defaults to 0. Other
  inversion methods keep using `smoothingValue*` as the correlation length
  of `Murat_smoothing`. New `bin/Murat_gradientOperator.m` and
  `bin/Murat_regulariser.m`; added `Tests/TestTikhonovSmoothing.m`.
- Optional `Murat.input.dampingSweepValues` and
  `bin/Murat_dampingRuns.m` for damping sweeps using Tikhonov, Particle,
  Annealing or Global inversion. Each value applies to Qc and Q at all
  configured frequencies, preserving frequency and smoothing settings.
  Waveform processing and forward modelling are shared between runs.
- Separate output folders and full model caches for each damping value.
  Before a new sweep, matching existing folders are archived as `_old`,
  `_old_2`, etc., preserving earlier results.
- Damping-sweep settings and usage guidance in the MSH, Romania and Toba
  input files, disabled by default with `dampingSweepValues = []`.
- `Murat.input.diffusionConstant` in the input files: `'Theoretical'` fixes
  the diffusion constant D = v/(3·B0·Le⁻¹) in the coda-normalisation
  inversion; `'Inverted'` (default) estimates D from the data and prints it
  next to the theoretical value.
- `bin/Murat_codaConstant.m` — coda constant per ray from the diffusion
  model, shared by `Murat_inversion`, `Murat_lsqlinQmean` and `Murat_plot`.

### Changed
- **Damping is now relative to the inversion matrix** (Qc and Q, all
  inversion methods). The penalty is `(dampingValue·‖A‖)²·‖m‖²`, with ‖A‖
  the largest singular value of the (weighted) matrix, so a value of 1
  weighs model norm and data misfit equally and the same value acts alike
  for Qc and Q. Before, the value was absolute: the Qc matrix has
  ‖A‖² ≈ 10⁹–10¹¹ for Toba, so damping values up to about 100 had no effect
  on Qc and damping sweeps returned identical models. Empty damping now
  defaults to 0.1 for both Qc and Q. The Tikhonov smoothing weight is
  relative to ‖A‖ in the same way. **Values from earlier input files are
  not interchangeable**; the MSH, Romania and Toba examples now use the
  default (`[]`). Added `Tests/TestRelativeDamping.m`.
- `Murat_tikhonov_nonneg_fista` starts from the regionalised starting
  model (`opts.x0`) instead of zero, stops on a relative step size, and
  warns (`MuRAT:fistaMaxIterations`) when it reaches `MaximumIterations`
  without converging. Relative damping below about 0.1 may need more
  iterations.
- `MuRAT.m` delegates nonempty damping sweeps to the new helper while
  preserving the original inversion, plotting and saving blocks for
  normal runs. An empty or absent sweep uses the existing per-frequency
  damping settings, including automatic selection when those are empty.
- MSH, Romania and Toba examples use the Tikhonov inversion method.
- `Murat_lsqlinQmean` removes the coda constant ray by ray and solves for the
  average Q⁻¹ by least squares. In `'Inverted'` mode, D is the value
  minimising the misfit (log grid over 1e-3–1e3 km²/s, refined locally).
  New signature: `[d0, Q0, mDest, cRay] = Murat_lsqlinQmean(cf_k, l,
  time0_k, rapsp_k, Qc_k, tCm, te, v, tl, Dfix)`.
- Checkerboard and spike tests for Q use the same per-ray coda constants as
  the data inversion.
- `Murat.Q.averageQ` (2×2×nFreq) and `Murat.Q.estimatedDiffusion` (1×nFreq)
  now store every frequency.
- Documentation (section 2.4.1) and wiki updated with the corrected CN
  equations.

### Fixed
- Clear the previous `Murat` configuration before loading a selected
  input file, preventing settings from a previous dataset or damping
  sweep from carrying over into the next run.
- Peak-delay colour limits in `Murat_image3D` are symmetric around zero,
  using the maximum absolute finite value in the volume at each frequency.
- `Murat_inversionQc` returns the exit flag in the Global solver branch
  when inversion plotting is disabled.
- Coda energy modelled as the window integral measured by `Murat_body`
  instead of the difference of diffusive energies at the window ends, which
  is negative once the coda decays and made the coda constant undefined.
- The coda constant now includes the direct-wave normalisation term
  −ln(4π·v·t_l)/(2πf).
- The grid search for D used `(tCm + te)` instead of `te` in the coda model.
- Comments in the MSH and Romania input files gave D = vS/(3·Le⁻¹·(1−B0));
  corrected to B0.
- `bin/Murat_codaCheck.m` — with `lapseTimeMethod = 'Travel'` the origin
  time was added twice, so for SAC files with a non-zero origin time (`o`)
  the coda window (used for Qc and the coda-normalization energy ratios)
  started `o` seconds too late. `tCoda` is now the lapse time from the
  origin for all methods. Added `Tests/TestCodaCheck.m`.
- `bin/Murat_unfoldXYZ.m` — failed with "Dimensions of arrays being
  concatenated are not consistent" when an axis had a single node, because
  `repelem` of a scalar returns a row. This happens in `Murat_kernels` for
  short lapse times relative to the model grid step. Added
  `Tests/TestUnfoldXYZ.m`.
- `bin/Murat_data.m`, `bin/Murat_dataParallelized.m` — `envelopeSmoothTime`
  from the input file is now passed to `Murat_envelope`; previously the 1 s
  default was always used. Inputs without the field keep the 1 s default.
- `Murat_inputToba.m`, `Murat_inputMSH.m`, `Murat_inputRomania.m` — the
  inversion section named the genetic-algorithm option `Genetic`; the code
  expects `'Global'`. Added `Tests/TestInputMethods.m`.
- `bin/Murat_smoothing.m` — coordinate columns were swapped; the function
  expected [lat lon z] but received [lon lat z], affecting non-Tikhonov
  smoothing. Fixed column assignment and updated documentation.
- `bin/Murat_checks.m` — third consecutive run deleted earlier results by
  unconditionally removing the `_old` backup folder. Now uses numbered
  backups (`_old1`, `_old2`, etc.) to preserve run history.
- `bin/Murat_inversion.m` — checkerboard and spike test inversions overwrote
  L-curve figures from main inversions. Each test type now saves to separate
  output directories (e.g., `L-curve_Qc_Checkerboard_*`, `L-curve_Q_Spike_*`).
- `bin/Murat_location.m` — undefined variable `listaSac_i` in error messages;
  now accepts `waveformName` parameter from callers. Also fixed error checking:
  replaced `find()` with `any()` for detecting missing -12345 sentinel values,
  and moved ×1000 depth conversion after validation.
- `bin/Murat_inversionQc.m` — Qc weight matrix applied incorrect double
  inversion and squaring. Correlation coefficients are now stored directly in
  `uncertaintyQc`, and weights are computed as correlation² for proper
  weighting in `Murat_plot.m` line 87.
- `bin/Murat_plot.m` — crash when peak-delay model contains no negative values
  (e.g., all positive anomalies). Added empty-array handling for threshold
  calculations.
- `bin/Murat_kernels.m` — Peak lapse-time method incorrectly averaged lapse
  times across frequencies; each frequency now retains its own lapse time for
  kernel computation.
- `bin/Murat_declustering.m` — index mismatch after `unique()` operation.
  Indices from the decimated array were used to index into the original
  uncertainty array, causing out-of-bounds errors. Fixed by tracking original
  indices.
- `bin/Murat_modv1D.m`, `bin/Murat_modv3D.m`, `bin/Murat_imageRays.m` —
  azimuth computation used `atan()` instead of `atan2()`, giving incorrect
  results in quadrants 2 and 3. Replaced with proper two-argument arctangent.
- `MuRAT.m` — `uigetfile` dialog blocked headless/batch operation. Now checks
  if `inputFile` is pre-defined; interactive mode uses the dialog only if no
  file is provided.
- Example input files (`Murat_inputToba.m`, `Murat_inputMSH.m`,
  `Murat_inputRomania.m`) — added commented example for `dampingSweepValues`
  to help users understand how to use damping sweeps with the default 0.1
  damping.

## [4.0.0] – 2026-10-06

### Added
- `CONTRIBUTING.md` — full contributor guidelines including bug reporting,
  feature requests, code style, test requirements, and third-party code rules.
- `THIRD_PARTY_NOTICES.md` — consolidated licence table for all third-party
  components in `Utilities_Matlab/`, clarifying notices for `inferno.m`,
  `redblue.m`, and `inpaintn.m`, and documenting the GPL v3 resolution.
- `CITATION.cff` — machine-readable citation metadata with 23 named authors
  and ORCIDs, linked to the MuRAT 4.0 Zenodo DOI.
- `AUTHORS.md` — lists principal developer and all contributors.
- GitHub Actions CI (`run_test.yml`) — runs `matlab-actions/run-tests@v2` on
  every push and pull request; produces JUnit XML and Cobertura coverage
  artefacts; back-tested on MATLAB R2021b.
- `Tests/TestRepositoryStructure.m` — class-based tests verifying critical
  files exist, all `.m` sources parse without errors, no editor backup files
  are committed, no duplicate `Murat_*` function names exist, and key
  regression guards on `Murat_testData` and `Murat_testAll` hold.
- `Tests/TestSacIo.m` — tests that `fget_sac` is on the path and that SAC
  header fields required by MuRAT are populated; skips gracefully when sample
  data folders are absent.
- `Tests/TestParallel.m` — smoke-tests `parfor` and verifies computed results;
  skips when Parallel Computing Toolbox is not installed.
- `bin/Murat_testData.m` — new `getFieldValue` helper reads SAC header fields
  by dot-path string without `eval()`; returns the SAC header struct as a
  third output (`sacHeader`).
- `Murat_checkerboard3D.m` (in `bin/`) — MIT-licensed, drop-in replacement
  for the GPL v3 `GIBBON/checkerBoard3D.m`; removes the only GPL dependency.
- `fillNaN3` local function inside `bin/Murat_rescale.m` — replaces
  `GIBBON/inpaintn.m` using MATLAB's built-in `fillmissing`; no external
  dependency.
- Selective `addpath` in `MuRAT.m` and `Murat_checks.m` — filters out `.git`,
  `.svn`, `private`, and hidden folders; only adds folders containing `.m`
  files.
- Automatic save format selection in `MuRAT.m` — uses `-v7.3` only when the
  `Murat` struct exceeds 2 GB; otherwise uses the faster `-v7` format.
- `Tests/TestFunctional.m` — toolbox-free functional smoke test verifying
  the data-reading pipeline (fget_sac → Murat_test → Murat_testData) on
  the MSH sample dataset; skips gracefully when data is absent
- `paper/` — JOSS submission files (paper.md, paper.bib, Fig1–3.png)

### Changed
- `MuRAT.m` — parallel pool now started correctly with `parpool(useParallel)`;
  `warning` in the catch block now passes `ME.identifier` and `ME.message`
  correctly; Live Script appendix metadata restored.
- `bin/Murat_testData.m` — `flag` encoding changed from integer accumulator
  (`flag = 0`, `+1`, `+2`) to boolean sentinels (`origMissing + 2*sMissing`);
  function signature extended to three outputs `[muratHeader, flag, sacHeader]`;
  `createsList` updated to single output; mandatory fields now raise `error`
  instead of `warning`.
- `bin/Murat_checks.m` — `createsList` call updated to single output; `flag`
  warning block updated to cover all four states (0, 1, 2, 3).
- `Tests/` — migrated from ad-hoc scripts (`test.m`, `test_sac_io.m`,
  `test_syntax_check.m`) to proper `matlab.unittest.TestCase` classes.
- `Utilities_Matlab/MyUtilities/Murat_testAll.m` — `createsList` call updated
  to single output.
- `README.md` — updated DOI badge to MuRAT 4.0 Zenodo record; expanded
  contributing and testing sections; added `CHANGELOG` and `CITATION` links.
- `_config.yml` — description field cleaned to plain prose.

### Removed
- `Utilities_Matlab/GIBBON/checkerBoard3D.m` and `iseven.m` — replaced by
  `bin/Murat_checkerboard3D.m` (MIT); eliminates GPL v3 incompatibility.
- `Utilities_Matlab/GIBBON/inpaintn.m` — replaced by `fillNaN3` in
  `bin/Murat_rescale.m`; eliminates unlicensed external dependency.
- `Utilities_Matlab/MyUtilities/Murat_test.m` — duplicate of `bin/Murat_test.m`;
  `bin/` is the canonical location per `CONTRIBUTING.md`.
- `Tests/test.m`, `Tests/test_sac_io.m`, `Tests/test_syntax_check.m` — replaced
  by the class-based test suite.
- `bin/Murat_testData.m.bak` — stale editor backup file.
- `sac_Romania/.DS_Store` — accidentally committed macOS metadata file.
- `Utilities_Matlab/regtu/licenseInverse.txt` — stale GPL v3 text with no
  referent in `regtu/`; removed to avoid licence confusion.

### Fixed
- `bin/Murat_testData.m` — `eval()` removed from all three pick-field guards;
  `getFieldValue` used consistently for both reading and the `-12345` check.
- `MuRAT.m` parallel pool — `pool = useParallel` (no-op) replaced by
  `parpool(useParallel)`; `warning` format string corrected.

### Migrating a 3.0 input file
3.0 inputs (`Murat_input*.mlx`) do not run unchanged in 4.0. Start from a 4.0
example input and copy values across:

| 3.0 field | 4.0 field | Note |
|---|---|---|
| `extinctionLength` (scalar) | `iExtinctionLength` | One value per central frequency. |
| `albedo` (scalar) | `albedo` | One value per central frequency. |
| `lCurve` | `PlotInversion` | `1` shows L-curves and solver progress. |
| `lCurveQc`, `lCurveQ` | `dampingValueQc`, `dampingValueQ` | One value per frequency; `[]` selects a default. Values are not interchangeable with 3.0 because the solvers changed. |
| — | `smoothingValueQc`, `smoothingValueQ` | New; `[]` or `0` disables smoothing. |
| `highCheck`, `lowCheck` | removed | Checkerboard values are now the min and max of the measured inverse Qc (reused for Q). |
| `spectralDecay` | removed | `bin/Murat_Qc.m` uses a fixed geometrical-spreading exponent of 1.5. |
| `inversionMethod = 'Iterative'` or `'Hybrid'` | `'Tikhonov'`, `'Particle'`, `'Annealing'` or `'Global'` | IRtools-based solvers removed. |
| — | `MaximumIterations`, `MaximumStallIterations` | Solver limits. |
| — | `envelopeSmoothTime` | Envelope smoothing window (s); 1 s matches 3.0. |
| — | `PlotRays`, `PlotTests`, `PlotResults`, `PlotCheckers`, `PlotSpikes` | Per-figure plot switches. |

---

## [3.0.0] – 2026-01-20

Legacy release. Source archived at Zenodo:
Luca De Siena et al. (2026). LucaDeSiena/MuRAT: Legacy MuRAT3.0 Code
(v3.26.01.20). <https://doi.org/10.5281/zenodo.18314469>

[4.0.0]: https://github.com/LucaDeSiena/MuRAT/compare/v3.26.01.20...v4.0.0
[3.0.0]: https://github.com/LucaDeSiena/MuRAT/releases/tag/v3.26.01.20
