MuRAT - Multi-Resolution Seismic Attenuation Tomography
=======

A Matlab Package for Seismic Attenuation Tomography at multiple Earth scales using Body and Coda Waves 

MuRAT is a complete software package for measuring seismic attenuation, scattering, and absorption from passive and active data, and map 2D and 3D variations of these parameters in space.

MuRAT also comes in Python and C++ versions (under development) to provide a a fully-integrated solution for seismic attenuation imaging. 

MuRAT1.0 was first developed by Luca De Siena (Johannes Gutenberg University, Mainz, Germany) during its PhD at the INGV-Osservatorio Vesuviano, Italy, and published in 2014 while he was research assistant at the Westfälisches Wihelms Universität, Münster.

MuRAT2.0 is the the result of the activity of the Volcano Imaging group, led by De Siea during his stint as permanent Lecturer at the University of Aberdeen, UK.

The group of active users (providing questions, feedback, snippets of code) comprises De Siena's PhD students at the University of Aberdeen. It includes PhD students he co-supervises internationally. 

History
-------

* 2006-2010: The core of MuRAT (at the time named "Multi-scale reasonable attenuation tomography analysis") is based on Matlab, c++, csh and fortran codes developped at INGV-Osservatorio.
* 2010-2013: MuRAT1.0 is developed as a 3D direct-wave attenuation imaging Matlab-only code with the contribution of Christine Thomas and Richard Aster.
* 2014: MuRAT1.0 is published in De Siena et al. 2014, JVGR, with two sample datasets (Mount St. Helens and Vesuvius) [JVGR article](https://www.sciencedirect.com/science/article/abs/pii/S0377027314000961)
* 2019: MuRAT2.0 is released including 2D scattering/absorption mapping and kernel-based inversion and re-branded Multi-Resolution Seismic Attenuation Tomography. [GitHub Repository] (https://github.com/LucaDeSiena/MuRAT)


Documentation
-------------
The full documentation can be found on: ???.


Installation
------------

SYSTEM: The program works on a Macbook pro with High Sierra, Matlab R2017a.
Necessary Toolboxes: Signal Processing and Mapping (for latitute-longitude coordinates - Romania example).

Two sample datasets (Mount St. Helens and Romania) can be downloaded at https://doi.pangaea.de/10.1594/PANGAEA.893893 to test the code. A third sample input (Pollino) is provided and the related dataset may be requested to Luca De Siena.

The current version works following these steps:

1. Download the package at https://github.com/LucaDeSiena/MuRAT and unzip folder Utilities_Matlab.

2. Download the two sample datasets at https://doi.pangaea.de/10.1594/PANGAEA.893893. Unzip the MSH (Mount St. Helens) and Romania datasets and put the folders in Murat-master folder.

3. Run MuRAT5_3_6.m.


INSTRUCTIONS:

1. When asked insert the name of the input file desired (either Input_MSH.m or Input_Romania.m).

2. After the L curves are produced, pick the smoothing parameter.

3. When building your example, use one of the Input files as template. Read attentively the comments and edit only the required parameters. Start with a “pa=1” or “pa=2” analysis. “pa=3” only works with a suitable velocity model.


Citing MuRAT
--------------

If you use MuRAT for your research and publications, please consider citing it:


Disclaimer
----------

Although we have cross-checked the whole code, we cannot warranty it is exempt of bugs. The package is provided as-is, we will not be held responsible for any use you make of it, nor for the results and conclusions you may find using MuRAT.

Licence
-------

MuRAT is released under EUPL v1.1


Input file fields
-------

**Analysis = 1, 2, or 3**

-pa=1 - Pick delay and Qc without kernels - **De Siena et al. 2016 (EPSL)**

-pa=2 - Pick delay and Qc with kernels - **De Siena et al. 2017 (GRL)**

-pa=3 - Pick delay, Qc with kernels, and P/S wave attenuation with the coda normalisation method - **De Siena et al. 2017 (GRL) and 2017 (Scientific Reports)

------------
**Label**
RULE: Name of the folder containing output files and figures** 

------------
**Working Directory**
RULE:  Directory where MuRAT.m is contained - Default = ./

------------
**Data Folder**
RULE: Directory containing the SAC files - Default = ./sac_Label

------------
**Figures Format**
RULE: Figures' output format as per Matlab functions - Default -> jpeg
RULE: Available formats: tiff,...

------------
**Figures Visibility**
RULE: Set visibility during computation, set on or off - Default -> on

------------
**Figures Sections**
RULE: Sections in the 3D V and Q models
RULE: Only for Analysis=3

------------
**P or S**
RULE: Available indexes: 2 or 3 - Default -> 2
RULE: Index one in time files contains the origin time of events
RULE: Index two in time files contains the P-wave phase -> 2
RULE: Index three in time files contains the S-wave phase -> 3

------------
**Registration Components**
RULE: Available indexes: 1, 2 or 3 - Default -> 1
RULE: Works with one vertical (1) or two horizontal (2) recordings, or with the three components(3) of motion

------------
**Spectrogram**
RULE: Available indexes: 0 or any number n between the data range - Default -> 0
RULE: n corresponds to the SAC recording of the spectrogram you want to visualize

------------
**Central frequency**
RULE: This is the frequency where the analyses will be carried on, depending on data recorded

------------
**Maximum Peak Delay Time**
RULE: This is the maximum time from picking where the maximum of the envelope will be searched for in the peak delay analysis

------------
**Minimum Peak Delay Time**
RULE: This is the minimum time from picking where the maximum of the envelope will be searched for in the peak delay analysis

------------
**Starting Lapse Time**
RULE: This is the time where the coda window for the coda attenuation analysis starts

------------
**Coda Window Length**
RULE: This is the length of the coda window for the coda attenuation analysis

------------
**Average Velocity**
RULE: This is the average velocity of the medium, in case the origin time is absent from recordings or a velocity model is unavailable
RULE: Default -> 0

------------
**SAC Origin Time**
RULE: The origin time of the event
RULE: As best practice, this should be saved in the SAC file as "o" variable -> "SAChdr.time.o"
RULE: In case this is unavailable set it as "[]" (Default).

------------
**SAC P Time**
RULE: The P-wave time of the event - this parameter is compulsory.
RULE: As best practice, this should be saved in the SAC file as "a" variable -> "SAChdr.time.a"

------------
**SAC S Time**
RULE: The S-wave time of the event
RULE: As best practice, this should be saved in the SAC file as "t0" variable -> "SAChdr.time.t0"
RULE: In case this is unavailable set it as "[]" (Default).

------------
**Length Window for Body Wave and Noise**
RULE: Window length set to compute energy for P/S waves and noise
RULE: Only for Analysis=3

------------
**Start Window Noise**
RULE: Seconds before P-wave arrival, where the window to compute noise starts
RULE: Only for Analysis=3

------------
**Import Option for Events and Stations**
RULE: Available indexes: 1 or 2 - Default -> 2
RULE: It is necessary to define source and station locations for mapping: you can import event origin time and coords of event and station from and external .txt file (1) or from the SAC files directly (2).
RULE: Index one is the original format of MuRAT 1.0 and requires even.txt and staz.txt files in the format of the sample files:

Event:
% column (1) = twelve numbers for the origin time of the event (date+time in seconds)
% column (2) = UTM (WE) or latitude
% column (3) = UTM (SN) or longitude
% column (4) = Altitude above sea level in meters

Station:
% column (1) = Name of station (3 characters)
% column (2) = UTM (WE) or latitude
% column (3) = UTM (SN) or longitude
% column (4) = Altitude above sea level in meters

RULE: Index two is the ideal format, where event and station info are stored in the SAC header, in lat/long format.