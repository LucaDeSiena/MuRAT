# Third-Party Notices

MuRAT itself is released under the [MIT License](LICENSE.md). The folder `Utilities_Matlab/` also contains third-party code that is **not** covered by that licence and remains under the terms chosen by its authors. Each component keeps its original licence file and copyright notices; please do not remove them.

## Components in `Utilities_Matlab/`

| Folder | Contents | Authors / copyright | Licence |
|---|---|---|---|
| `MatSAC/` | `fget_sac.m`, `sachdr.m`, `sac.m`, `rdSac.m`, `rdSacHead.m`, `wtSac.m`, `newSacHeader.m`, `sacfft.m`, sample files `MYJH.*` and `N.MYJH.Z.sac` | `fget_sac.m` and `sachdr.m`: Zhigang Peng. `sac.m`: Xianglei Huang (03/2000). The author of the remaining routines is not recorded (the folder's `readme.txt` suggests they may come from Jeff McGuire, WHOI). Source: <http://geophysics.eas.gatech.edu/classes/SAC/> | **No licence stated** in the folder or file headers. Redistribution terms need to be confirmed with the authors. |
| `F_SAC/` | `fread_sac.m`, `fwrite_sac.m` | Whyjay Zheng, Copyright (c) 2015 | BSD 2-Clause (`F_SAC/license.txt`) |
| `regtu/` | `corner.m`, `fil_fac.m`, `l_corner.m`, `l_curve.m`, `l_curve_tikh_svd.m`, `lcfun.m`, `picard.m`, `plot_lc.m`, `tikhonov.m`: a subset of Regularization Tools 4.1 | Per Christian Hansen, DTU Compute, Copyright (c) 2015. Source: <https://www.mathworks.com/matlabcentral/fileexchange/52-regtools> | BSD 3-Clause (`regtu/license.txt`) |
| `COLORMAP/` | `colMapGen.m`, `inferno.m`, `redblue.m` | The folder's licence names Timothy Olsen, Copyright (c) 2019. The authors of `inferno.m` and `redblue.m` are not recorded in the files. | BSD 3-Clause (`COLORMAP/LICENSE`); the licence for `inferno.m` and `redblue.m` needs confirming |
| `GIBBON/` | `checkerBoard3D.m` and `iseven.m` (Kevin Mattheus Moerman), `inpaintn.m` (Damien Garcia, 2010-2017) | GIBBON: Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors. Source: <https://github.com/gibbonCode/GIBBON> | **GNU GPL v3** (`GIBBON/LICENSE`). `inpaintn.m` carries no licence text in its header. |
| `MyUtilities/` | `Murat_changeHdr.m`, `Murat_plotMore.m`, `Murat_testAll.m`, `freq_analysis.m`, `hitmap.m` | MuRAT authors (see [AUTHORS.md](AUTHORS.md)) | MIT, as MuRAT |

The folder `Utilities_Matlab/UTMDEG/` (coordinate conversion) keeps its own `license.txt`; add its author and licence here once confirmed.

## Notes

- **Regularization Tools** (`regtu/`): if you use these routines, please also cite P. C. Hansen, *Regularization Tools: A Matlab package for analysis and solution of discrete ill-posed problems*, Numerical Algorithms 6 (1994), pp. 1-35, and P. C. Hansen, *Regularization Tools Version 4.0 for Matlab 7.3*, Numerical Algorithms 46 (2007).
- **`regtu/licenseInverse.txt`** contains the full text of the GNU GPL v3. No file in that folder refers to it, so its purpose is unclear.
- **`inpaintn.m`**: please cite the references listed in the header of the file (Garcia, *Computational Statistics & Data Analysis*, 2010, and Wang et al.).

## Dependencies that are not distributed

- **MATLAB** and the MathWorks toolboxes listed in the [README](README.md) (proprietary software, licensed by MathWorks).
- **GeophysicalModelGenerator.jl** and **ParaView**, suggested for visualising the output models; they are used outside MuRAT and are not redistributed.

## Updating this file

When you add, update or remove third-party code, update the table above, keep the original licence file next to the code, and check that the licence permits redistribution with MuRAT. See [CONTRIBUTING.md](CONTRIBUTING.md).
