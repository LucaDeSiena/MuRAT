# MuRAT
Attenuation tomography code. 

PROGRAM MuRAT for Multi-resolution attenuation tomography analysis.

PROGRAM FOR: 3D direct- and 2D peak delay and coda-wave (with and without kernels) attenuation tomography

Two sample datasets (Mount St. Helens and Romania) have to be downloaded before running the code. They can be found as zipped folders at https://doi.pangaea.de/10.1594/PANGAEA.893893

Two other sample inputs are provided for the Tirreno Sea and Pollino area. The dataset may be requested to Luca De Siena.

SYSTEM: The program works on a Macbook pro with High Sierra, Matlab R2017a, Necessary toolboxes are: Signal Processing and Mapping Toolboxes (for latitute-longitude coordinates - Romania example).

INSTRUCTIONS:
The current version works with the following steps:

Download the package at https://github.com/LucaDeSiena/MuRAT and unzip folder Utilities_Matlab.

Download the two sample datasets at https://doi.pangaea.de/10.1594/PANGAEA.893893. Unzip the MSH (Mount St. Helens) and Romania datasets and put the folders in Murat-master folder.

Run MuRAT5_3_6.m.

When asked insert the name of the input file desired (either Input_MSH.m or Input_Romania.m).

After the L curves are produced, pick the smoothing parameter.

When building your example, use one of the Input files as template. Read attentively the comments and edit only the required parameters. Start with a “pa=1” or “pa=2” analysis. “pa=3” only works with a suitable velocity model.
