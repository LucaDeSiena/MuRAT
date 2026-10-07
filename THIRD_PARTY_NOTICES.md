# Third-Party Notices

MuRAT itself is released under the [MIT License](LICENSE.md). The folder
`Utilities_Matlab/` also contains third-party code that is **not** covered by
that licence and remains under the terms chosen by its authors. Each component
keeps its original licence file and copyright notices; please do not remove
them.

## Components in `Utilities_Matlab/`

| Folder | Contents | Authors / copyright | Licence |
|---|---|---|---|
| `MatSAC/` | `fget_sac.m` (replaces all previous MatSAC `.m` files), sample files `MYJH.*` and `N.MYJH.Z.sac` | MuRAT Authors, Copyright (c) 2026. The function replicates the SAC binary reader originally written by Zhigang Peng and Xianglei Huang; the reimplementation is independent and carries no code from those authors. | MIT (`LICENSE.md`) |
| `F_SAC/` | `fread_sac.m`, `fwrite_sac.m` | Whyjay Zheng, Copyright (c) 2015 | BSD 2-Clause (`F_SAC/license.txt`) |
| `regtu/` | `corner.m`, `fil_fac.m`, `l_corner.m`, `l_curve.m`, `l_curve_tikh_svd.m`, `lcfun.m`, `picard.m`, `plot_lc.m`, `tikhonov.m` — a subset of Regularization Tools 4.1 | Per Christian Hansen, DTU Compute, Copyright (c) 2015. Source: <https://www.mathworks.com/matlabcentral/fileexchange/52-regtools> | BSD 3-Clause (`regtu/license.txt`) |
| `COLORMAP/` | `colMapGen.m`, `inferno.m`, `redblue.m` | `colMapGen.m`: Timothy Olsen, Copyright (c) 2019 (UCSF). `inferno.m`: colormap data from the matplotlib/viscm project — Nathaniel J. Smith, Stefan van der Walt, Eric Firing (CC0/public domain); MATLAB wrapper author not recorded. `redblue.m`: Adam Auton, 9 October 2009. | `colMapGen.m`: BSD 3-Clause (`COLORMAP/LICENSE`). `inferno.m`: colormap data is CC0; wrapper carries no licence — `COLORMAP/LICENSE_inferno.txt` added as a clarifying notice. `redblue.m`: originally posted to MATLAB File Exchange under the BSD licence; `COLORMAP/LICENSE_redblue.txt` added. |
| `MyUtilities/` | `Murat_changeHdr.m`, `Murat_plotMore.m`, `Murat_testAll.m`, `freq_analysis.m`, `hitmap.m` | MuRAT authors (see [AUTHORS.md](AUTHORS.md)) | MIT, as MuRAT |

## Notes

### Regularization Tools (`regtu/`)

If you use any of these routines, please cite:

- P. C. Hansen, *Regularization Tools: A Matlab package for analysis and solution of discrete ill-posed problems*, Numerical Algorithms **6** (1994),  pp. 1–35.

- P. C. Hansen, *Regularization Tools Version 4.0 for Matlab 7.3*, Numerical
  Algorithms **46** (2007).

### GIBBON (resolved)

All three files previously taken from GIBBON have been replaced:

- `checkerBoard3D.m` and `iseven.m` → replaced by `bin/Murat_checkerboard3D.m`
  (MIT-licensed, drop-in compatible, ~30 lines).
- `inpaintn.m` → replaced by a `fillNaN3` local function inside
  `bin/Murat_rescale.m` using MATLAB's built-in `fillmissing` (R2016b+,
  no toolbox required). The `GIBBON/` folder and its licence files
  (`LICENSE`, `licenseBoilerPlate.txt`, `LICENSE_inpaintn.txt`) should now
  be deleted from the repository.

The GPL v3 compatibility issue is therefore fully resolved.

### `MatSAC/`

The original eight `.m` files from Zhigang Peng's SAC tutorial (five of which
had unknown provenance) have been replaced by a single MIT-licensed
reimplementation: `fget_sac.m` (MuRAT Authors, 2026). The new file reads the
same SAC binary format and returns an identical `SAChdr` struct, so all calling
code is unchanged. The sample SAC files (`MYJH.*`, `N.MYJH.Z.sac`) are kept
for testing and are not executable code.

### `inferno.m`

The 256-entry RGB table is the matplotlib `inferno` colormap, released under
CC0 (public domain) by Nathaniel J. Smith, Stefan van der Walt, and Eric
Firing. The MATLAB wrapper function has no stated author or licence; a
clarifying notice is provided in `COLORMAP/LICENSE_inferno.txt`.

### `redblue.m`

Written by Adam Auton (2009) and originally posted to MATLAB File Exchange.
At the time of submission MATLAB File Exchange applied a BSD licence to all
submissions by default. A clarifying notice is provided in
`COLORMAP/LICENSE_redblue.txt`.

## Dependencies that are not distributed

- **MATLAB** and the MathWorks toolboxes listed in the [README](README.md)
  (proprietary software, licensed separately by MathWorks).
- **GeophysicalModelGenerator.jl** and **ParaView**, suggested for visualising
  output models; they are used outside MuRAT and are not redistributed.

## Updating this file

When you add, update, or remove third-party code: update the table above, keep
the original licence file next to the code, and verify that the licence permits
redistribution alongside MuRAT. See [CONTRIBUTING.md](CONTRIBUTING.md).
