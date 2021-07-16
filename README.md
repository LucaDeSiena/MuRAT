![Campi Flegrei direct-wave attenuation model.](https://github.com/LucaDeSiena/MuRAT/blob/Murat_ID/MuRAT.jpeg)

MuRAT - Multi-Resolution Seismic Attenuation Tomography Analysis
=======

MuRAT is a Matlab Package for Seismic Attenuation, Scattering and Absorption Tomography using Body and Coda Waves. 

MuRAT measures seismic attenuation, scattering, and absorption from passive and active data, and models 3D variations of these parameters in space.

The group of active users (providing questions, feedback, snippets of code) is the [Volcano Earth Imaging group](https://www.lucadesiena.com).

*History*
-------

* *2006-2010*: MuRAT is built first by Luca De Siena during his PhD at the INGV-Osservatorio Vesuviano (Italy) using Matlab, c++, csh and Fortran codes.

* *2010-2013*: MuRAT1.0 was developed and published in 2014 while De Siena was research assistant at the Westfälisches Wilhelms Universität, Münster (Germany). Murat1.0 allowed 3D total attenuation imaging with the coda-normalization method. Important contributions were given by Christine Thomas (WWU Münster) and Richard Aster (Colorado State University).

* *2014*: MuRAT1.0 is published in [De Siena et al. 2014, JVGR](https://www.sciencedirect.com/science/article/abs/pii/S0377027314000961) with two sample datasets (Mount St. Helens and Vesuvius).

* *2018*: [MuRAT2D](https://github.com/LucaDeSiena/MuRAT2D) is the result of the activity of the [Volcano Earth Imaging group](https://www.lucadesiena.com), led by De Siena during his stint as Lecturer at the University of Aberdeen (UK). It images 2D seismic scattering (peak delay) and absorption (Qc at late lapse time - kernel-based). It is ideal for small datasets and in case no previous velocity information is available.

* *2021*: [MuRAT3D](https://github.com/LucaDeSiena/MuRAT) is released as a parallelized code for full 3D attenuation, scattering and absorption imaging using peak delays, coda attenuation and coda-normalized energies.

*Documentation*
-------------

The full documentation for MuRAT3.0 can be found in html files associated with the software. This README file and the linked internet sites are to be used as a reference. 

*System*
------------
The program works on Mac, Linux and Windows systems equipped with Matlab R2019a.

Necessary Toolboxes: Signal Processing, Curve Fitting, Image Processing and Mapping. The Parallel Computing Toolbox is recommended for speed.

Two sample datasets (Mount St. Helens and Romania) are included and allow the user to obtain sample models. The  datasets work with the two input .mlx files that are provided and show examples of what the user can obtain with the code. 

*Instructions in a nutshell*
------------

The current version works following these steps:

1. Download or clone the package at https://github.com/LucaDeSiena/MuRAT.

2. Work in the downloaded folder after moving it to an appropriate location on your system.

3. Open one of the two input .mlx files, providing a step-by-step explanation of all inputs (Murat_input_MSH.mlx or Murat_input_Romania.mlx) and create your own. 

4. MuRAT works with [SAC files](https://ds.iris.edu/files/sac-manual/) that must be stored into a single folder and corrected for the instrument function. The files must have populated headers, although the code can work using just the following header fields:

              a. The P-wave picking in the reference time of the waveform.
              b. The coordinates of the event.
              c. The coordinates of the station.
              d. The origin time of the event (optional).

5. Run MuRAT3 and select the name of the input file desired.

*What the code does*
--------

To understand what MuRAT3D does:

1. **Start from the Murat_input..mlx files**

The input files are self-explanatory and provide detailed descriptions of every input and references to papers you can use to set them. If you have a 3D velocity model use *MuRAT_InputMSH.mlx* otherwise start from *MuRAT_InputRomania.mlx*.

2. **Read the html**

The package has a html folder where each file explains what the code does: read ethem to understand the approximations used to process data, forward model kernels, and invert observations. Particularly important are:

* MuRAT3.html: is the *make* code and can be run in sections. Here you find the names of the primary functions.
* Murat_plot.htlm: the plot function shows at which lines of the code figures are produced, providing additional information on the output.

3. **Understand the output text files**

All the output files (.txt), figures and .vtk files (for visualisation in Paraview) are stored in sub-directories in the **Label** folder, created in the working directory. Use the html Murat_plot html file to have information about what each plot means and how it is created. In the following, a list of the output files and what they contain is provided.

Inside the *TXT* subfolder, the first three columns of each output file correspond to WE, SN, and depth. The fourth column is the mapped parameter. In ascii format, they contain a minimum of five columns (for *Peak Delay*) that can be imported to show the locations of the anomalies in a simple (x,y,z) reference system. The fifth columns shows blocks hit by at least one ray.

*Qc.txt* and *Q3D.txt* are solved with an inversion and thus have:

      A. a sixth and seventh columns that corresponds to the input and output of the checkerboard test;
      B. eight and nineth columns that corresponds to the input and output of the spike test;
      C. a 10th column for the model resolution matrix.

All the .vtk files are stored into the *VTK* subfolder. A Matlab structure ('Murat.mat') containing all inputs and data is stored in the working directory.

4. **Understand the output figure files**

All the figures (in the **Figures Format** defined by the user) are stored in subdirectories in the **Label** folder, created in the working directory. 

*Rays.Figures format*

A figure to show how rays develop in 3D for the Peak Delay and Q measurements. It plots them on three slices (WE, SN, Z). The fourth panel shows the location of the area on the Earth.

------------
*Kernel-Qc.Figures format* 

Two panels showing the sensitivity kernels in the entire 3D space (left) and the normalised kernels in the chosen inversion grid (right). This reduction implies several hypotheses: among these the most important is that most of the energy is still comprised in the grid (the difference is general < 1% if all source and stations are in the inversion grid.

------------
*Qc_Peak_Delay.Figures format*

A figure to evaluate the appropriate peak-delay and coda inputs. The upper panel should show a constant inverse Qc with travel time. The lower panel should show a linearly-increasing peak delay with travel time.

------------
*Lc_Qc.Figures format and Lc_CN.Figures format*

L-curves corresponding to the coda-attenuation and total-attenuation inversion. After they appear, a prompt asks which damping parameter the user wants to pick. 

------------
*Picard_Qc.Figures format and Picard_CN.Figures format*

These plots show the result of the Picard criterium, necessary to evaluate how many of the inversion parameters are correctively solved in the coda-attenuation and total-attenuation inversions, respectively. The two figures do not appear during computation.

------------
*Peak-delay-3D.fig*, *Qc-3D.fig* and *Q-3D.fig*

These plots show the result of the peak-delay, Qc and Q 3D tomography in the grid's reference system. All in Matlab .fig format, use the .vtk and Paraview for publication-quality figures.

------------
*Qc-checkerboard.fig*, *Qc-spike.fig*

These plots show input and output of the checkerboard and spike tests for the Qc and Q mapping in the grid's reference system.

------------
*Parameter_space_variations.Figures format*

The plot shows the separation of the scattering and absorption parameters in their parameter space. Grey dots correspond to parameters too near to the average to be interpreted as scattering or absorption variations - the threshold is pre-defined at 5% of the maximum variation of each parameter. Red = High scattering and absorption; Cyan = High scattering and low absorption; Orange = Low scattering and high absorption; Green = Low scattering and absorption.

------------
*Parameter-Map.fig*

The parameter space separation as it is apparent in the 3D space. Each block is characterized by the color corresponding to its scattering and absorption characteristics

------------
*V_model.fig*

The 3D velocity model is also available as 3D figures in Matlab format. They can be loaded in Matlab and will show the vertical and horizontal slices defined in **Figures Sections**.

*Citing MuRAT*
------------

If you use MuRAT for your research and publications, please consider mentioning the GitHub internet site and citing the following papers, depending on the techniques you are going to use


*Q (Total attenuation)*:

1. De Siena, L., C. Thomas, and R. Aster. "Multi-scale reasonable attenuation tomography analysis (MuRAT): An imaging algorithm designed for volcanic regions." Journal of Volcanology and Geothermal Research 277 (2014): 22-35. - **Older release that discusses the code for coda-normalisation, also used in the early works of Prudencio et al. 2015,a,b, GJI**

2. De Siena, L., G. Chiodini, G. Vilardo, E. Del Pezzo, M. Castellano, S. Colombelli, N. Tisato, and G. Ventura, 2017. Source and dynamics of a volcanic caldera unrest: Campi Flegrei, 1983–84. Scientific reports: Nature Journals 7, 8099. - **Recent implementation of the Coda Normalization method with correction for coda attenuation variations**

3. Sketsiou P., L. De Siena, S. Gabrielli, F. Napolitano, 2021. 3-D attenuation image of fluid storage and tectonic interactions across the Pollino fault network. Geophysical Journal International, 226(1), 536-547. - **Most recent application of MuRAT**

*Qc and Peak Delay (Absorption and scattering)*:

1. De Siena L., Calvet, M., Watson, K.J., Jonkers, A.R.T. and Thomas, C., 2016. Seismic
scattering and absorption mapping of debris flows, feeding paths, and tectonic units at Mount St. Helens volcano. Earth and Planetary Science Letters, 442, pp.21-31. - **Implementation of the older peak delay and Qc technique, both with regionalisation**

2. De Siena L., A. Amoruso, E. Del Pezzo, Z. Wakeford, M. Castellano, L. Crescentini, 2017.
Space-weighted seismic attenuation mapping of the aseismic source of Campi Flegrei 1983–84
unrest. Geophysical Research Letters, 44.4 pp. 1740-1748. - **First implementation with kernels for Qc**

3. Del Pezzo, E., De La Torre, A., Bianco, F., Ibanez, J., Gabrielli, S., and De Siena, L. (2018). Numerically Calculated 3D Space-Weighting Functions to Image Crustal Volcanic Structures Using Diffuse Coda Waves. - **Numerical implementation of kernel functions**

4. Sketsiou P., F. Napolitano, A. Zenonos, L. De Siena, (2020). New insights into seismic absorption imaging. Physics of the Earth and Planetary Interiors, 298, 106337. - **Comprehensive review of the method and future outlooks**

*Disclaimer*
------------

Although we have cross-checked the whole code, we cannot warranty it is exempt of bugs. The package is provided as-is, we will neither be held responsible for any use you make of it nor for the results and conclusions you may derive using MuRAT.

*Licence*
------------

MuRAT is released under EUPL v1.1
