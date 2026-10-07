# Changelog

All notable changes to MuRAT are documented here.
Dates follow ISO 8601 (YYYY-MM-DD). Versions follow [Semantic Versioning](https://semver.org).

---

## [4.0.0] – 2026-10-06

### Added
- `CONTRIBUTING.md` — full contributor guidelines including bug reporting,
  feature requests, code style, test requirements, and third-party code rules.
- `THIRD_PARTY_NOTICES.md` — consolidated licence table for all third-party
  components in `Utilities_Matlab/`, clarifying notices for `inferno.m`,
  `redblue.m`, and `inpaintn.m`, and documenting the GPL v3 resolution.
- `CITATION.cff` — machine-readable citation metadata with 19 named authors
  and ORCIDs, linked to the MuRAT 4.0 Zenodo DOI.
- `AUTHORS.md` — lists principal developer and all contributors.
- GitHub Actions CI (`run_test.yml`) — runs `matlab-actions/run-tests@v2` on
  every push and pull request; produces JUnit XML and Cobertura coverage
  artefacts; pinned to MATLAB R2025b.
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

---

## [3.0.0] – 2026-01-20

Legacy release. Source archived at Zenodo:
Luca De Siena et al. (2026). LucaDeSiena/MuRAT: Legacy MuRAT3.0 Code
(v3.26.01.20). <https://doi.org/10.5281/zenodo.18314469>

[4.0.0]: https://github.com/LucaDeSiena/MuRAT/compare/v3.26.01.20...v4.0.0
[3.0.0]: https://github.com/LucaDeSiena/MuRAT/releases/tag/v3.26.01.20
