# Third-Party Notices

MuRAT itself is released under the [MIT License](LICENSE.md). The folder `Utilities_Matlab/` also contains third-party code that is **not** covered by that licence and remains under the terms chosen by its authors. Each component keeps its original licence file and copyright notices; please do not remove them.

## Components in `Utilities_Matlab/`

| Folder | Contents | Authors / copyright | Licence |
|---|---|---|---|
| `MatSAC/` | `fget_sac.m`, `sachdr.m`, `sac.m`, `rdSac.m`, `rdSacHead.m`, `wtSac.m`, `newSacHeader.m`, `sacfft.m`, sample files `MYJH.*` and `N.MYJH.Z.sac` | `fget_sac.m` and `sachdr.m`: Zhigang Peng. `sac.m`: Xianglei Huang (03/2000). The author of the remaining routines is not recorded (the folder's `readme.txt` suggests they may come from Jeff McGuire, WHOI). Source: <http://geophysics.eas.gatech.edu/classes/SAC/> | **No licence stated** in the folder or file headers. Redistribution terms need to be confirmed with the authors. |
| `F_SAC/` | `fread_sac.m`, `fwrite_sac.m` | Whyjay Zheng, Copyright (c) 2015 | BSD 2-Clause (`F_SAC/license.txt`) |
| `regtu/` | `corner.m`, `fil_fac.m`, `l_corner.m`, `l_curve.m`, `l_curve_tikh_svd.m`, `lcfun.m`, `picard.m`, `plot_lc.m`, `tikhonov.m`: a subset of Regularization Tools 4.1 | Per Christian Hansen, DTU Compute, Copyright (c) 2015. Source: <https://www.mathworks.com/matlabcentral/fileexchange/52-regtools> | BSD 3-Clause (`regtu/license.txt`) |
| `COLORMAP/` | `colMapGen.m`, `inferno.m`, `redblue.m` | `colMapGen.m`: Timothy Olsen, Copyright (c) 2019. `inferno.m`: the colormap data originates from the matplotlib/viscm project (Nathaniel J. Smith, Stefan van der Walt, Eric Firing); the MATLAB wrapper carries no explicit copyright or licence header — redistribution terms need confirming with the uploader. `redblue.m`: Adam Auton, 9 October 2009; no licence stated in the file. | `colMapGen.m`: BSD 3-Clause (`COLORMAP/LICENSE`). `inferno.m` and `redblue.m`: **no licence stated** in the files — redistribution terms need confirming. |
| `GIBBON/` | `checkerBoard3D.m` and `iseven.m` (Kevin Mattheus Moerman), `inpaintn.m` (Damien Garcia, 2010–2017) | GIBBON: Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors. Source: <https://github.com/gibbonCode/GIBBON>. `inpaintn.m`: Damien Garcia, 2010/06, last update 2017/08; originally published on MATLAB File Exchange — **no licence text** in the file header or in this folder. | **GNU GPL v3** (`GIBBON/LICENSE` and `GIBBON/licenseBoilerPlate.txt`). `inpaintn.m`: licence not stated in the distributed copy — redistribution terms need confirming with the author. |
| `MyUtilities/` | `Murat_changeHdr.m`, `Murat_plotMore.m`, `Murat_test.m`, `Murat_testAll.m`, `freq_analysis.m`, `hitmap.m` | MuRAT authors (see [AUTHORS.md](AUTHORS.md)) | MIT, as MuRAT |

## Notes

- **Regularization Tools** (`regtu/`): if you use these routines, please also cite P. C. Hansen, *Regularization Tools: A Matlab package for analysis and solution of discrete ill-posed problems*, Numerical Algorithms 6 (1994), pp. 1-35, and P. C. Hansen, *Regularization Tools Version 4.0 for Matlab 7.3*, Numerical Algorithms 46 (2007).
- **`regtu/licenseInverse.txt`** contains the full text of the GNU GPL v3. No file in `regtu/` refers to it; its presence is unexplained. If none of the Regularization Tools files are actually GPL-licensed, this file should be removed to avoid confusion.
- **`inpaintn.m`**: please cite the references listed in the header of the file — Garcia, *Computational Statistics & Data Analysis*, 2010, and Wang et al.
- **`GIBBON/` and GPL compatibility**: the GIBBON files (`checkerBoard3D.m`, `iseven.m`) are GPL v3. GPL v3 code cannot be redistributed as part of an MIT-licensed project without the combined work also being GPL v3. If MuRAT is intended to remain MIT-licensed, these files should either be replaced with MIT/BSD-compatible alternatives or kept strictly as an optional add-on that is not distributed alongside the main package.
- **`MatSAC/`**: no licence is stated for any file in this folder. The `readme.txt` (last updated 2006 by Zhigang Peng) acknowledges that the provenance of `rdSac.m`, `rdSacHead.m`, `wtSac.m`, `newSacHeader.m`, and `sacfft.m` is uncertain ("maybe Jeff McGuire at WHOI"). Redistribution terms for all files in this folder need to be confirmed before the next release.
- **`inferno.m`**: the colormap data matches the matplotlib `inferno` colormap (Nathaniel J. Smith, Stefan van der Walt, Eric Firing), which is CC0 / public domain. However, the MATLAB wrapper file itself carries no copyright or licence statement. The uploader should be identified and a licence confirmed.
- **`redblue.m`**: written by Adam Auton (2009), originally posted to MATLAB File Exchange under the BSD licence at the time of submission. No licence file accompanies the copy in this repository; a licence file should be added.
- **`UTMDEG/`**: referenced in the previous version of this file but the folder is **not present** in the repository. If it was removed, this notice is now accurate. If it should be present, it needs to be re-added together with its licence.

## Dependencies that are not distributed

- **MATLAB** and the MathWorks toolboxes listed in the [README](README.md) (proprietary software, licensed by MathWorks).
- **GeophysicalModelGenerator.jl** and **ParaView**, suggested for visualising the output models; they are used outside MuRAT and are not redistributed.

## Updating this file

When you add, update or remove third-party code, update the table above, keep the original licence file next to the code, and check that the licence permits redistribution with MuRAT. See [CONTRIBUTING.md](CONTRIBUTING.md).
