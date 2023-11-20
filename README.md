MuRAT - Multi-Resolution seismic Attenuation Tomography
=======

![MuRAT is a code for attenuation, scattering and absorption tomography.](./img/muratlogo.jpg)

[![test](https://github.com/deconvolution/MuRAT/actions/workflows/run_test.yml/badge.svg?branch=with_test)](https://github.com/deconvolution/MuRAT/actions/workflows/run_test.yml)

MuRAT is a Matlab Package for seismic Attenuation, Scattering and Absorption Tomography using Body and Coda Waves at multiple frequencies.

MuRAT measures seismic attenuation, scattering, and absorption from passive and active data and models 3D variations of these parameters in space.

The group of active users (providing questions, feedback, and snippets of code) is the [Volcano Earth Imaging group](https://www.lucadesiena.com).

*Documentation*
-------------

The file Documentation.pdf in this folder serves as complete documentation for MuRAT3.0. This README file and the *Input_.mlx* files in this folder act as additional documentation.

The Wiki for MuRAT is under construction, but you can already check a bit of the history of the code.

*System*
------------

The program works on Mac, Linux and Windows systems with Matlab version R2017a or higher.

Necessary Toolboxes: Signal Processing, Curve Fitting, Image Processing and Mapping. The Parallel Computing Toolbox is recommended for speed.

Custom toolboxes not included in standard Matlab installations are also provided with the package. These are:

1. Routines to read SAC files created by Zhigang Peng and available from [his SAC tutorial page](http://geophysics.eas.gatech.edu/classes/SAC/).
2. The [Regularization Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/52-regtools?s_tid=prof_contriblnk) was created by Per Christian Hansen and available from Matlab File Exchange.
3. The [IRTools](https://github.com/jnagy1/IRtools/tree/ebd70d4036c3cd8c82fc1e17033351491fddf11f), included in MuRAT as a zipped folder.
4. Functions from the [Geometry and Image-Based Bioengineering add-On for MATLAB](https://github.com/gibbonCode/GIBBON).

Three sample datasets (Mount St. Helens, Romania, and Toba) are included and allow the user to obtain sample models. The datasets work with the three corresponding *input.mlx* files that show examples of what the user can get with the code.

*Instructions in a nutshell*
------------

The current version works following these steps:

1. Download or clone the package at <https://github.com/LucaDeSiena/MuRAT>.

2. Work in the downloaded folder after moving it to an appropriate location on your system.

3. Check that the IRTools have been downloaded as a zipped folder in the corresponding folder in the working directory. Otherwise, download them from <https://github.com/jnagy1/IRtools/tree/ebd70d4036c3cd8c82fc1e17033351491fddf11f>.

4. Open one of the three input .mlx files, providing a step-by-step explanation of all inputs (*Murat_inputMSH.mlx*, *Murat_inputRomania.mlx*, or *Murat_inputToba.mlx*) and create your own.

5. Use a velocity model, storing it in the corresponding folder. The format is [Latitude, Longitude, Altitude (meters)]

6. MuRAT works with [SAC files](https://ds.iris.edu/files/sac-manual/) that must be stored in a single folder and corrected for the instrument function. The files must have populated headers. Your SAC headers get tested anyway; the result is shown in an Excel file. The code takes from the header the following fields:
***a)*** The P-wave picking in the reference time of the waveform (in seconds);
***b)*** The coordinates of the event in degrees - beware, *the earthquake depth must be in kilometers*;
***c)*** The coordinates of the station - beware, *the station elevation must be in meters*;
***d)*** The origin time of the event (optional) in seconds.

7. Run MuRAT3 and select the name of the input file desired.

*Workflow*
--------

A. ***Start from the Murat_input..mlx files***

The input files are self-explanatory and provide detailed descriptions of every input and references to papers you can use to set them. If you have a 3D velocity model, use *MuRAT_InputMSH.mlx* otherwise start from either *MuRAT_InputRomania.mlx* or *MuRAT_InputToba.mlx*, the examples for 3-component data.

B. ***Read the Documentation***

The Documentation includes a summary of the theory underlying attenuation imaging: read it to understand the approximations used to process data, forward model kernels, and invert observations.

C. ***Understand the output text files***

All the output files (.mat, .txt and .xlsx) and figures are stored in the **TXT** sub-directories in the **Label** folder, created in the working directory. A list of the output files and what they contain is provided in the following. Ascii files contain the models in degrees and UTM. We strongly suggest imaging the TXT files using the [GeophysicalModelGenerator](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl).

D. **Understand the output figure files**

Beware, *.fig* figures are created with the invisible option in Matlab. There are two ways to open them:

(1) Use the function *openfig(..,'visible')* to open them from the command window.

(2) Click twice on the figure file and write *shg* in the command window.

All the figures are stored in subdirectories in the **Label** folder, created in the working directory:

*Structure of the Label Folder*
--------

------------

* **TXT directory**

------------

*peakdelay__.txt*, *Qc__.txt* and *Q__.txt*:  The 3D models of the parameters at different frequencies. The first three columns of all text files correspond to Latitude, Longitude, and altitude. The fourth column is the mapped parameter. They contain a minimum of five columns (for *Peak Delay*) that can be imported to show the locations of the anomalies in a simple (x,y,z) reference system. The fifth column shows blocks hit by at least one ray. Qc and Q are solved with an inversion and thus have (1) the sixth and seventh columns that correspond to the input and output of the checkerboard test; (2) the eighth and ninth columns that correspond to the input and output of the spike test.

*Murat.mat*: A Matlab structure containing all inputs and data the code produces.

*DataHeaders.xls*: A file containing all header variables of the SAC files used for the mapping, useful for data selection.

------------

* ***Checkerboard directory***

------------

      Qc subdirectory

*Qc-Checkerboard__.tif* and *Qc-Checkerboard__.fig*: These figures show the input and output of the Qc checkerboard test in the 3D space (*.fig*) and through cross-sections (*.tif*).

      Q subdirectory

*Q-Checkerboard__.tif* and *Q-Checkerboard__.fig*: These figures show the input and output of the Q checkerboard test in the 3D space (*.fig*) and through cross-sections (*.tif*).

------------

* ***RaysKernels directory***

------------

*Clustering.tif*: This figure shows all rays used on the map (black, discarded) with those after declustering (red). 

*Rays__.tif*: These figures show how rays develop in 3D for the Peak Delay and Q measurements. It plots them on three slices (WE, SN, Z). The fourth panel shows the location of the area on the Earth.

*Kernel__.tif* and *Kernel__.fig*: Each *.fig* figure has two panels showing the sensitivity kernels in the entire 3D space (left) and the normalised kernels in the chosen inversion grid (right). This reduction implies several hypotheses: the most important is that most energy is still in the grid (the difference is generally < 1% if all sources and stations are in the inversion grid). The *.tif* figures are sections in the WE, SN, and Z directions. Figures are produced for all frequencies.

------------

* ***Results directory***

------------

      Parameter subdirectory

*Parameter__.tif* and *Parameter__.fig*: Parameter maps in 3D (*.fig*) and across sections (*.tif*).

      PeakDelay subdirectory

*Peak-Delay__.tif* and *Peak-Delay__.fig*: Peak delay maps in 3D (*.fig*) and across sections (*.tif*).

      Q subdirectory

*Q__.tif* and *Q__.fig*: Total attenuation maps in 3D (*.fig*) and across sections (*.tif*).

      Qc subdirectory

*Qc__.tif* and *Qc__.fig*: Coda attenuation maps in 3D (*.fig*) and across sections (*.tif*).

  ------------

* ***Spike directory***

------------
      Qc subdirectory

*Qc-Spike__.tif* and *Qc-Spike__.fig*: These figures show the input and output of the Qc spike test in the 3D space (*.fig*) and through cross sections (*.tif*).

      Q subdirectory

*Q-Spike__.tif* and *Q-Spike__.fig*: These figures show the input and output of the Q spike test in the 3D space (*.fig*) and through cross sections (*.tif*).

------------

* ***Tests directory***

------------

*Qc_Analysis__.tif*, *PD_Analysis__.tif*, and *CN_Analysis__.tif*

Three figures to evaluate the appropriate peak-delay and coda inputs. Read the documentation for further clarifications.

*L_curve__.fig*: L-curves and cost functions (depending on the inversion method) for the Qc and Q inversions necessary to set the damping parameters. The user can ask for a prompt or set the damping parameters.

*Qc_analysis__*: Relationship between coda attenuation and frequency.

*Velocity_model.fig*: The 3D velocity model is also available as a figure in Matlab format. They can be loaded in Matlab and show the vertical and horizontal slices defined in *Figures Sections*.

------------

*Citing MuRAT*
------------

If you use MuRAT for your research and publications, please consider mentioning the GitHub internet site and citing the following papers, depending on the techniques you are going to use

**Q (Total attenuation)**:

1. De Siena, L., C. Thomas, and R. Aster. "Multi-scale reasonable attenuation tomography analysis (MuRAT): An imaging algorithm designed for volcanic regions." Journal of Volcanology and Geothermal Research 277 (2014): 22-35. - *Older release that discusses the code for coda-normalisation, also used in the early works of Prudencio et al. 2015, a,b, GJI*

2. De Siena, L., G. Chiodini, G. Vilardo, E. Del Pezzo, M. Castellano, S. Colombelli, N. Tisato, and G. Ventura, 2017. Source and dynamics of a volcanic caldera unrest: Campi Flegrei, 1983–84. Scientific reports: Nature Journals 7, 8099. - *Recent implementation of the Coda Normalization method with correction for coda attenuation variations*

3. Sketsiou P., L. De Siena, S. Gabrielli, F. Napolitano, 2021. 3-D attenuation image of fluid storage and tectonic interactions across the Pollino fault network. Geophysical Journal International, 226(1), 536-547. - *Most recent application of Q imaging with MuRAT*

**Qc and Peak Delay (Absorption and scattering)**:

1. De Siena L., Calvet, M., Watson, K.J., Jonkers, A.R.T. and Thomas, C., 2016. Seismic scattering and absorption mapping of debris flows, feeding paths, and tectonic units at Mount St. Helens volcano. Earth and Planetary Science Letters, 442, pp.21-31. - *Implementation of the older peak delay and Qc technique, both with regionalisation*

2. De Siena L., A. Amoruso, E. Del Pezzo, Z. Wakeford, M. Castellano, L. Crescentini, 2017.
Space-weighted seismic attenuation mapping of the aseismic source of Campi Flegrei 1983–84
unrest. Geophysical Research Letters, 44.4 pp. 1740-1748. - *First implementation with kernels for Qc*

3. Del Pezzo, E., De La Torre, A., Bianco, F., Ibanez, J., Gabrielli, S., and De Siena, L. (2018). Numerically Calculated 3D Space-Weighting Functions to Image Crustal Volcanic Structures Using Diffuse Coda Waves. - *Numerical implementation of kernel functions*

4. Sketsiou P., F. Napolitano, A. Zenonos, L. De Siena, (2020). New insights into seismic absorption imaging. Physics of the Earth and Planetary Interiors, 298, 106337. - *Comprehensive review of the method and future outlooks*

*Disclaimer*
------------

Although we have cross-checked the whole code, we cannot warranty it is exempt from bugs. The package is provided as-is; we will neither be held responsible for any use you make of it nor for the results and conclusions you may derive using MuRAT.

*Licence*
------------

MuRAT is released under EUPL v1.1

*Funding*
------------

Some developments of this software package were funded by the Deutsche Forshungsgemeinshaft under grant number SI1748/4-1.

