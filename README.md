MuRAT - Multi-Resolution Seismic Attenuation Tomography
=======

Getting Started
-------

MuRAT is a Matlab Package for Seismic Attenuation Tomography at multiple Earth scales using Body and Coda Waves 


*Features*
--------

* MuRAT is a complete software package for measuring seismic attenuation, scattering, and absorption from passive and active data, and map 2D and 3D variations of these parameters in space.

* MuRAT also comes in Python and C++ versions (under development) to provide a a fully-integrated solution for seismic attenuation imaging. 

* MuRAT1.0 was first developed by Luca De Siena (Johannes Gutenberg University, Mainz, Germany) during its PhD at the INGV-Osservatorio Vesuviano, Italy, and published in 2014 while he was research assistant at the Westfälisches Wihelms Universität, Münster.

* MuRAT2.0 is the the result of the activity of the Volcano Imaging group, led by De Siea during his stint as permanent Lecturer at the University of Aberdeen, UK.

* The group of active users (providing questions, feedback, snippets of code) comprises De Siena's PhD students at the University of Aberdeen. It includes PhD students he co-supervises internationally. 

*History*
-------

* 2006-2010: The core of MuRAT (at the time named "Multi-scale reasonable attenuation tomography analysis") is based on Matlab, c++, csh and fortran codes developped at INGV-Osservatorio.

* 2010-2013: MuRAT1.0 is developed as a 3D direct-wave attenuation imaging Matlab-only code with the contribution of Christine Thomas and Richard Aster.

* 2014: MuRAT1.0 is published in De Siena et al. 2014, JVGR, with two sample datasets (Mount St. Helens and Vesuvius) [JVGR article](https://www.sciencedirect.com/science/article/abs/pii/S0377027314000961)

* 2019: MuRAT2.0 is released including 2D scattering/absorption mapping and kernel-based inversion and re-branded Multi-Resolution Seismic Attenuation Tomography. [GitHub Repository] (https://github.com/LucaDeSiena/MuRAT)

*Documentation*
-------------

The full documentation for MuRAT2.0 is under construction. This README file and the linked internet sites are to be used as reference. 

*Installation and running*
------------

SYSTEM: The program works on a Macbook pro with High Sierra, Matlab R2017a.

Necessary Toolboxes: Signal Processing (compulsory) and Mapping (optional, for geolocalisation - Romania example).

Two sample datasets (Mount St. Helens and Romania) can be downloaded at https://doi.pangaea.de/10.1594/PANGAEA.893893 to test the code. A third sample input (Pollino) is provided and the related dataset may be requested to Luca De Siena.

The current version works following these steps:

1. Download the package at https://github.com/LucaDeSiena/MuRAT.

2. Download the two sample datasets at https://doi.pangaea.de/10.1594/PANGAEA.893893. Unzip the MSH (Mount St. Helens) and Romania datasets and put the folders in Murat-master folder. Delete the zipped files.

3. MuRAT works with a single input file in .xsl format (for now). Build the Input file starting from the samples (either Input_MSH.xsl or Input_Romania.xsl). Read the instructions presented in the next sections.

4. Run MuRAT2.m.

5. After the L curves are produced, pick the smoothing parameter.

6. When building your example, use one of the Input files as template. Read attentively the comments and edit only the required parameters. Start with a “pa=1” or “pa=2” analysis. “pa=3” only works with a suitable velocity model.

*Instructions - the input file*
------------

**Analysis**

*Available indexes: 1, 2, or 3 - Default -> 1*

Analysis = 1 performs a pick-delay and Qc analysis without kernels, as *De Siena et al. 2016 (EPSL)*

Analysis = 2 performs a pick-delay and Qc analysis with kernels, as *De Siena et al. 2017 (GRL)*

Analysis = 3 performs a pick delay, Qc with kernels, and P/S wave attenuation with the coda normalisation method, as *De Siena et al. 2017 (GRL) and 2017 (Scientific Reports)* combined.

------------
**Label**

Name of the folder containing output files and figures, will appear in the working directory

------------
**Working Directory**

Directory where MuRAT.m is contained - *Default = ./*

------------
**Data Folder**

Directory containing the SAC files - *Default = ./sac_Label*

------------
**Figures Format**

Figures' output format as per Matlab functions - *Default -> jpeg*

------------
**Figures Visibility**

Set visibility during computation, set on or off - *Default -> on*

------------
**Figures Sections**

Sections in the 3D V and Q models. Only available for Analysis=3

------------
**P or S**

*Available indexes: 2 or 3 - Default -> 2*

Index one in the time files contains the origin time of events. Index two in the time files contains the P-wave phase. Index three in the time files contains the S-wave phase.

------------
**Registration Components**

*Available indexes: 1, 2 or 3 - Default -> 1*

Works with one vertical (1) or two horizontal (2) recordings, or with the three components(3) of motion.

------------
**Spectrogram**

*Available indexes: 0 or any number n in the data range - Default -> 0*

*n* corresponds to the SAC recording of the spectrogram you want to visualize.

------------
**Central frequency**

This is the frequency where the analyses will be carried on, depending on data recorded.

------------
**Maximum Peak Delay Time**

This is the maximum time from picking where the maximum of the envelope will be searched for in the peak delay analysis

------------
**Minimum Peak Delay Time**

This is the minimum time from picking where the maximum of the envelope will be searched for in the peak delay analysis

------------
**Starting Lapse Time**

This is the time where the coda window for the coda attenuation analysis starts.

------------
**Coda Window Length**

This is the length of the coda window for the coda attenuation analysis.

------------
**Average Velocity**

This is the average velocity of the medium, in the case that the origin time is absent from recordings or a velocity model is unavailable - *Default -> 0*

------------
**SAC Origin Time**

The origin time of the event. As best practice, this should be saved in the SAC file as *o* variable -> *SAChdr.time.o*. If this is unavailable set it as *Default -> \[]*

------------
**SAC P Time**

The P-wave time of the event - this parameter is compulsory. As best practice, this should be saved in the SAC file as *a* variable -> *SAChdr.time.a*

------------
**SAC S Time**

The S-wave time of the event. As best practice, this should be saved in the SAC file as *t0* variable -> *SAChdr.time.t0*. If this is unavailable set it as *Default -> \[]*

------------
**Length Window for Body Wave and Noise**

Window length set to compute energy for P or S waves and noise. Only available for Analysis=3

------------
**Start Window Noise**

Seconds before P-wave arrival, where the window to compute noise starts. Only available for Analysis=3

------------
**Import Option for Events and Stations**

*Available indexes: 1 or 2 - Default -> 2*

It is necessary to define source and station locations for mapping: you can import event origin time and coords of event and station from and external .txt file (1) or from the SAC files directly (2). Index one is the original format of MuRAT 1.0 and requires even.txt and staz.txt files as per sample files. Index two is the ideal format, where event and station info are stored in the SAC header, in lat/long format. The field that must be filled are:

*Event*: Name of event - *SAChdr.event.kevnm*; Station latitude - *SAChdr.event.evla*; Station longitude- *SAChdr.event.evlo*; Station depth - *SAChdr.event.evdp*

*Station*: Name of station - *SAChdr.station.kstnm*; Station latitude - *SAChdr.station.stla*; Station longitude- *SAChdr.station.stlo*; Station elevation - *SAChdr.station.stel*

-------------
**Origin of the Grid**

The origin of the 2D or 3D grid for imaging is either in UTM or lat/long with depths in meters or km

-------------
**Name Events File**

Name of the file containing events' coordinates, for Import Option = 1. The format is:

*even.txt*: Column (1) = twelve numbers for the origin time of the event (date+time in seconds); Column (2) = UTM (WE) or latitude; Column (3) = UTM (SN) or longitude; Column (4) = Depth in meters or km

-------------
**Name Stations File**

Name of the file containing stations' coordinates, for Import Option = 1. The format is:

*staz.txt*: Column (1) = Max four characters; Column (2) = UTM (WE) or latitude; Column (3) = UTM (SN) or longitude; Column (4) = Depth in meters or km

-------------
**Number of X Nodes**, **Number of Y Nodes** and **Grid Step**

Number of X and Y layers and step of the mapping grid. Must be defined by the user to include all events' and stations' locations

-------------
**Meters or Degrees**

Can be set to either 1 (meters) or 111 (meters in degrees)

-------------
**Set 1 to Create Rays**

Can be set to either 0 (does not create) or 1 (creates) to build a unique variable containing all rays in 3D. This will be stored in the Murat structure in the Murat.geometry.createrays field. NOT RECOMMENDED as the variable will become BIG for tomographic datasets, but can be useful as input for other programs.

-------------
**Depth Synthetic Model Layer**

Set the depth (negative) to create a two-layer synthetic test. Only available for Analysis=3.

-------------
**Dimension V**

*Available indexes: 1 or 3*

Sets the dimension of the velocity model. In the case no 3D model is available, *iasp91* is used. Only available for Analysis=3.

-------------
**Name Velocity Model**

Name of the file containing the 3D velocity model, if avialable. Only available for Analysis=3.

-------------
**Size Checkerboard Anomalies**

*Available indexes: 2 or 4*

Sets the dimension of the checkerboard anomalies to either twice or four times that of the imaging grid.

-------------
**Highs Checkerboard Anomalies** and **Lows Checkerboard Anomalies**

RULE: Sets the high-attenuation and low-attenuation values (as inverse Q and inverse Qc) for the checkerboard tests.

-------------
**Linear or Nonlinear**

*Available indexes: 0 or 1 - Default -> 0*

Inverts coda attenuation with either a linearised (0) or non-linear (1) approach.

-------------
**Treshold Linear Fit**

The maximum inverse-Qc uncertainty used to weight the coda attenuation data. To be set to an appropriate parameter for the linearised inversion

-------------
**Length Nonlinear Fit Window** and **Number Nonlinear Fit Window**

The length and number of windows used to compute coda energies for the inversion. E.g., for a total coda window of 15 seconds, the user can set either 3 windows of 5 seconds each or 5 windows of 3 seconds each. To be set to appropriate parameters for the non-linear inversion

-------------
**Minimum Nonlinear Inverse Qc**, **Maximum Nonlinear Inverse Qc**, and **Total Nonlinear Inverse Qc**

The search for the Qc minimising the inversion is done starting from a minimum inverse Q (e.g. 0) to a maximum inverse Qc (e.g. 0.01). The total number of Qc we search are equally-spaced and defined between minimum and maximum. To be set to appropriate parameters for the non-linear inversion

*Instructions - the output files*
------------

All the output files (.txt) and figures (in the format defined by the user) are stored in the **Label** folder, created in the **Working Directory**. The first three columns of the output files correspond to WE, SN, and depth. The fourth column is the mapped parameter. 

------------
*peakdelay.txt*

Variations of Log10 Peak delay with respect to the average - the parameter used to map scattering attenuation. The depth is set to -1000 m. Only four columns as the parameters are just regionalised, not inverted.

------------
*Qc.txt*

Variations of inverse Qc, alias coda attenuation - the parameter used to map absorption. The depth is set to -1000 m.
    
**Analysis = 1**: Only four columns as the parameters are just regionalised, not inverted.

**Analyses = 2 or 3**: The fifth and sixth columns correspond to input and output of the checkerboard test, respectively. The seventh and eigth columns are the first two 2D kernels.

------------
*Q3D.txt* - Only for **Analysis = 3**

Variations of inverse Q, alias total direct-wave attenuation. The fifth and sixth columns correspond to input and output of the 3D checkerboard test, respectively. The seventh column provides the values of the resolution matrix. The eigth column is the result of a synthetic test using the actual result as input - this is improper, but used in some studies. The ninth and tenth columns are the input and output of the two-layer synthetic test, where the interface is set at **Depth Synthetic Model Layer**. 

*Instructions - the output figures*
------------

All the figures (in the **Figures Format** defined by the user) are stored in the **Label** folder, created in the **Working Directory**. If the Mapping Toolbox is available and the coordinates are in latitude and longitude, the figures will show coast lines.

------------
*Rays.Figures format*

A figure to plot rays either in 2D (Analyses + 1 or 2) or 3D (Analysis = 3) - it will show rays in the reference system of the pre-defined grid with effective sources and stations used.

------------
*Qc_sensitivity.Figures format* - Only available for **Analyses = 2 or 3**

A figure to plot the source-station normalised kernel for the first source-station pair. It will show the sensitivity in the reference system of the pre-defined grid.

------------
*Spectrogram.Figures format - Optional*

An optional figure to plot the first spectrogram in the dataset. The plot can be used to set appropriate frequency, peak delay, and coda inputs.

------------
*Qc_Peak_Delay.Figures format*

A figure to evaluate how well the peak-delay and coda inputs have been set. The upper panel should show a constant inverse Qc with travel time. The lower panel should show an increasing peak delay with travel time. Red dots correspond to outliers. For Q, the maximum uncertainty is defined in **Treshold Linear Fit**, for the non-linear it is pre-defined. For the peak-delay, we define as outliers all the values over twice the standard deviation with respect to the average.

------------
*CN_plot.Figures format - Optional*

A figure to evaluate if the direct-energy decrease with travel time follows the assumption of the coda-normalisation method (upper panel). Differently from the standard coda-normalisation method, we input the Qc measurements obtained in previous steps in the equation. This removes the effect of coda attenuation variations, improving the fit for the direct wave attenuation. The lower panel shows the effect of resigual geometrical spreading on coda attenuation, with the evental presence of outliers.

------------
*Lc_Qc.Figures format (for **Analyses = 2 or 3**) and Lc_CN.Figures format (for **Analyses = 2 or 3**)*

These are the only plots appearing if **Figures Visibility = 0** during computation. They show the L-curves corresponding to the coda-attenuation and total-attenuation inversion. After they appear, a prompt asks which damping parameter the user wants to pick.

------------
*Picard_Qc..Figures format (for **Analyses = 2 or 3**) and Picard_CN.Figures format (for **Analysis = 3**)*

These plots show the result of the Picard analysis, necessary to evaluate how many of the inversion parameters are correctively solved in the coda-attenuation and total-attenuation inversions, respectively. The two figures are not shown during computation.

------------
*Peak_delay_map..Figures format* and *Qc_map.Figures format*

These plots show the result of the peak-delay and Qc mapping in the grid reference system.

------------
*Qc_checkerboard_input.Figures format* and *Qc_checkerboard_output.Figures format - Only for **Analyses = 2 or 3***

These plots show the result of the checkerboard test for the Qc mapping in the grid reference system.

------------
*Parameter_space_variations.Figures format*

Separation of scattering and absorption parameters in their parameter space. Gray dots correspond to parameters too near to the average to be interpreted as scattering or absorption variations - the treshold is pre-defined at 5% of the maximum variation of each parameter. Red = High scattering and absorption; Cyan = High scattering and low absorption; Orange = Low scattering and high absorption; Green = Low scattering and absorption.

------------
*Q_model.fig; V_model.fig; 3DQ_check_input.fig; 3DQ_check_output.fig - Only for **Analysis = 3***

The 3D velocityy model, the total attenuation model, and the checkerboard test input and output are available as 3D figures in Matlab format. They can be loaded in matlab, and will show the vertical and horizontal slices defined in **Figures Sections**.

*Citing MuRAT*
------------

If you use MuRAT for your research and publications, please consider mentioning the GitHub internet site and citing the first release of the code, published as:

*De Siena, L., C. Thomas, and R. Aster. "Multi-scale reasonable attenuation tomography analysis (MuRAT): An imaging algorithm designed for volcanic regions." Journal of volcanology and geothermal research 277 (2014): 22-35.*


*Disclaimer*
------------

Although we have cross-checked the whole code, we cannot warranty it is exempt of bugs. The package is provided as-is, we will not be held responsible for any use you make of it, nor for the results and conclusions you may find using MuRAT.

*Licence*
------------

MuRAT is released under EUPL v1.1
