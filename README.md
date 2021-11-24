![MuRAT is a code for attenuation, scattering and absorption tomography.](https://github.com/LucaDeSiena/MuRAT/blob/master/MuRAT.jpg)

[![test](https://github.com/deconvolution/MuRAT/actions/workflows/run_test.yml/badge.svg?branch=with_test)](https://github.com/deconvolution/MuRAT/actions/workflows/run_test.yml)

MuRAT - Multi-Resolution seismic Attenuation Tomography
=======

MuRAT is a Matlab Package for seismic Attenuation, Scattering and Absorption Tomography using Body and Coda Waves at multiple frequencies.

MuRAT measures seismic attenuation, scattering, and absorption from passive and active data, and models 3D variations of these parameters in space.

The group of active users (providing questions, feedback, snippets of code) is the [Volcano Earth Imaging group](https://www.lucadesiena.com).

*History*
-------

* *2006-2010*: MuRAT was built first by Luca De Siena during his PhD at the INGV-Osservatorio Vesuviano (Italy) using Matlab, c++, csh and Fortran codes.

* *2010-2013*: MuRAT1.0 was developed and published in 2014 while De Siena was research assistant at the Westfälisches Wilhelms Universität, Münster (Germany). Murat1.0 allowed 3D total attenuation imaging with the coda-normalization method. Important contributions were given by Christine Thomas (WWU Münster) and Richard Aster (Colorado State University).

* *2014*: MuRAT1.0 is published in [De Siena et al. 2014, JVGR](https://www.sciencedirect.com/science/article/abs/pii/S0377027314000961) with two sample datasets (Mount St. Helens and Vesuvius).

* *2018*: [MuRAT2D](https://github.com/LucaDeSiena/MuRAT2D) is the result of the activity of the [Volcano Earth Imaging group](https://www.lucadesiena.com), led by De Siena during his stint as Lecturer at the University of Aberdeen (UK). It images 2D seismic scattering (peak delay) and absorption (Qc at late lapse time - kernel-based). It is ideal for small datasets and in case no previous velocity information is available.

* *2021*: [MuRAT3D](https://github.com/LucaDeSiena/MuRAT) is released as a multi-frequency parallelized code for full 3D attenuation, scattering and absorption imaging using peak delays, coda attenuation and coda-normalized energies.

*Documentation*
-------------

The file Documentation.pdf in this folder serves as full documentation for MuRAT3.0. This README file and the *Input_.mlx* files in this folder act as additional documentation.

*System*
------------

The program works on Mac, Linux and Windows systems equipped with Matlab version R2017a or higher.

Necessary Toolboxes: Signal Processing, Curve Fitting, Image Processing and Mapping. The Parallel Computing Toolbox is recommended for speed.

Custom toolboxes not included in standard Matlab installations are also provided with the package. These are:

1. Routines to read SAC files created by Zhigang Peng and available from [his SAC tutorial page](http://geophysics.eas.gatech.edu/classes/SAC/).
2. The [Regularization Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/52-regtools?s_tid=prof_contriblnk) created by Per Christian Hansen and available from Matlab File Exchange.
3. The [IRTools](https://github.com/jnagy1/IRtools/tree/ebd70d4036c3cd8c82fc1e17033351491fddf11f), included in MuRAT as a zipped folder.
4. Functions from the [Geometry and Image-Based Bioengineering add-On for MATLAB](https://github.com/gibbonCode/GIBBON).

Three sample datasets (Mount St. Helens, Romania, and Toba) are included and allow the user to obtain sample models. The  datasets work with the three corresponding *input.mlx* files that show examples of what the user can obtain with the code.

*Instructions in a nutshell*
------------

The current version works following these steps:

1. Download or clone the package at <https://github.com/LucaDeSiena/MuRAT>.

2. Work in the downloaded folder after moving it to an appropriate location on your system.

3. Check that the IRTools have been downloaded as a zipped folder in the corresponding folder in the working directory. Otherwise download them from <https://github.com/jnagy1/IRtools/tree/ebd70d4036c3cd8c82fc1e17033351491fddf11f>.

4. Open one of the three input .mlx files, providing a step-by-step explanation of all inputs (*Murat_inputMSH.mlx*, *Murat_inputRomania.mlx*, or *Murat_inputToba.mlx*) and create your own.

5. Use a velocity model, storing it in the corresponding folder. The format is [Latitude, Longitde, Altitude (meters)]

6. MuRAT works with [SAC files](https://ds.iris.edu/files/sac-manual/) that must be stored into a single folder and corrected for the instrument function. The files must have populated headers, although the code can work using just the following header fields:

              a. The P-wave picking in the reference time of the waveform.
              b. The coordinates of the event.
              c. The coordinates of the station.
              d. The origin time of the event (optional).

Test your SAC headers with the functions Murat\_test and Murat\_testAll in the folder **Utilities**.

7. Run MuRAT3 and select the name of the input file desired.

*What the code does*
--------

To understand what MuRAT3D does:

  1. ***Start from the Murat_input..mlx files***

The input files are self-explanatory and provide detailed descriptions of every input and references to papers you can use to set them. If you have a 3D velocity model use *MuRAT_InputMSH.mlx* otherwise start from either *MuRAT_InputRomania.mlx* or *MuRAT_InputToba.mlx*, the examples for 2- and 3-component data.

  2. ***Read the Documentation***

The Documentation includes a summary of the theory underlying attenuation imaging: read it to understand the approximations used to process data, forward model kernels, and invert observations.

  3. ***Understand the output text files***

All the output files (.mat, .txt and xlsx), figures and .vtk files (for visualisation in Paraview) are stored in the **TXT** and **VTK** sub-directories in the **Label** folder, created in the working directory. In the following, a list of the output files and what they contain is provided.

  4. **Understand the output figure files**

Beware, *.fig* figures are created with the invisible option in Matlab. Use the function *openfig(..,'visible','on')* to open them from the command window. All the figures are stored in subdirectories in the **Label** folder, created in the working directory:

*Structure of the Label Folder, where results are stored*
--------

  1. **TXT directory** and **VTK directory**

*peakdelay__.txt*, *Qc__.txt* and *Q__.txt*:  The 3D models of the parameters at different frequencies. The first three columns of all text files correspond to WE, SN, and altitude. The fourth column is the mapped parameter. They contain a minimum of five columns (for *Peak Delay*) that can be imported to show the locations of the anomalies in a simple (x,y,z) reference system. The fifth columns shows blocks hit by at least one ray. Qc and Q are solved with an inversion and thus have: (1) sixth and seventh columns that corresponds to the input and output of the checkerboard test; (2) eight and ninth columns that corresponds to the input and output of the spike test. All the .vtk files are stored in omonimous folder.

*Murat.mat*: A Matlab structure containing all inputs and data produced by the code.

*DataHeaders.xls*: A file containing all headers variables of the SAC files used for the mapping.

  ------------

  2. **Checkerboard directory**

      Qc subdirectory

*Qc-Checkerboard__.tif* and *Qc-Checkerboard__.fig*: These figures show input and output of the Qc checkerboard test in the 3D space (*.fig*) and across sections (*.tif*).

      Q subdirectory

*Q-Checkerboard__.tif* and *Q-Checkerboard__.fig*: These figures show input and output of the Q checkerboard test in the 3D space (*.fig*) and across sections (*.tif*).

  ------------

  3. **RaysKernels directory**

*Rays__.tif*: These figures show how rays develop in 3D for the Peak Delay and Q measurements. It plots them on three slices (WE, SN, Z). The fourth panel shows the location of the area on the Earth.

*Kernel__.tif* and *Kernel__.fig*: Each *.fig* figure has two panels showing the sensitivity kernels in the entire 3D space (left) and the normalised kernels in the chosen inversion grid (right). This reduction implies several hypotheses: among these the most important is that most of the energy is still comprised in the grid (the difference is general < 1% if all source and stations are in the inversion grid. The *.tif* figures are sections in the WE, SN, and Z directions. Figures are produced for all frequencies.

  ------------

  4. **Results directory**

      Parameter subdirectory

*Parameter__.tif* and *Parameter__.fig*: Parameter maps in 3D (*.fig*) and across sections (*.tif*).

      PeakDelay subdirectory

*Peak-Delay__.tif* and *Peak-Delay__.fig*: Peak delay maps in 3D (*.fig*) and across sections (*.tif*).

      Q subdirectory

*Q__.tif* and *Q__.fig*: Total attenuation maps in 3D (*.fig*) and across sections (*.tif*).

      Qc subdirectory

*Qc__.tif* and *Qc__.fig*: Coda attenuation maps in 3D (*.fig*) and across sections (*.tif*).

  ------------

  5. **Spike directory**

      Qc subdirectory

*Qc-Spike__.tif* and *Qc-Spike__.fig*: These figures show input and output of the Qc spike test in the 3D space (*.fig*) and across sections (*.tif*).

      Q subdirectory

*Q-Spike__.tif* and *Q-Spike__.fig*: These figures show input and output of the Q spike test in the 3D space (*.fig*) and across sections (*.tif*).

  ------------

  6. **Tests directory**

*Qc_Analysis__.tif*, *PD_Analysis__.tif*, and *CN_Analysis__.tif*

Three figures to evaluate the appropriate peak-delay and coda inputs. Read the documentation for further clarifications.

*L_curve__.fig*: L-curves and cost functions (depending on inversion method) for the Qc and Q inversions necessary to set the damping parameters. The user can ask for a prompt or set the damping parameters from start.

*Qc_vs_frequency*: Relationship between coda attenuation and frequency.

*Velocity_model.fig*: The 3D velocity model is also available as a figure in Matlab format. They can be loaded in Matlab and will show the vertical and horizontal slices defined in *Figures Sections*.

*Citing MuRAT*
------------

If you use MuRAT for your research and publications, please consider mentioning the GitHub internet site and citing the following papers, depending on the techniques you are going to use

**Q (Total attenuation)**:

1. De Siena, L., C. Thomas, and R. Aster. "Multi-scale reasonable attenuation tomography analysis (MuRAT): An imaging algorithm designed for volcanic regions." Journal of Volcanology and Geothermal Research 277 (2014): 22-35. - *Older release that discusses the code for coda-normalisation, also used in the early works of Prudencio et al. 2015,a,b, GJI*

2. De Siena, L., G. Chiodini, G. Vilardo, E. Del Pezzo, M. Castellano, S. Colombelli, N. Tisato, and G. Ventura, 2017. Source and dynamics of a volcanic caldera unrest: Campi Flegrei, 1983–84. Scientific reports: Nature Journals 7, 8099. - *Recent implementation of the Coda Normalization method with correction for coda attenuation variations*

3. Sketsiou P., L. De Siena, S. Gabrielli, F. Napolitano, 2021. 3-D attenuation image of fluid storage and tectonic interactions across the Pollino fault network. Geophysical Journal International, 226(1), 536-547. - *Most recent application of Q imaging with MuRAT*

**Qc and Peak Delay (Absorption and scattering)**:

1. De Siena L., Calvet, M., Watson, K.J., Jonkers, A.R.T. and Thomas, C., 2016. Seismic
scattering and absorption mapping of debris flows, feeding paths, and tectonic units at Mount St. Helens volcano. Earth and Planetary Science Letters, 442, pp.21-31. - *Implementation of the older peak delay and Qc technique, both with regionalisation*

2. De Siena L., A. Amoruso, E. Del Pezzo, Z. Wakeford, M. Castellano, L. Crescentini, 2017.
Space-weighted seismic attenuation mapping of the aseismic source of Campi Flegrei 1983–84
unrest. Geophysical Research Letters, 44.4 pp. 1740-1748. - *First implementation with kernels for Qc*

3. Del Pezzo, E., De La Torre, A., Bianco, F., Ibanez, J., Gabrielli, S., and De Siena, L. (2018). Numerically Calculated 3D Space-Weighting Functions to Image Crustal Volcanic Structures Using Diffuse Coda Waves. - *Numerical implementation of kernel functions*

4. Sketsiou P., F. Napolitano, A. Zenonos, L. De Siena, (2020). New insights into seismic absorption imaging. Physics of the Earth and Planetary Interiors, 298, 106337. - *Comprehensive review of the method and future outlooks*

*Disclaimer*
------------

Although we have cross-checked the whole code, we cannot warranty it is exempt of bugs. The package is provided as-is, we will neither be held responsible for any use you make of it nor for the results and conclusions you may derive using MuRAT.

*Licence*
------------

MuRAT is released under EUPL v1.1
