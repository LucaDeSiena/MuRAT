%% MuRAT   Multi-Resolution (seismic) Attenuation Tomography
%%
%   SCOPE: 3D total attenuation, scattering, and absorption tomography
%
%%   
% SYSTEM: The program works on all Mac and Linux computers where it has
% been tried.
%%
% MATLAB Version: The code keeps the us of toolboxes to a minimum. However,
% there are four necessary toolboxes: Signal Processing, Curve Fitting,
% Image Processing and Mapping Toolboxes. The parallel computing toolbox is
% highly recommended for speed.
%%
% Two sample datasets (Mount St. Helens and Romania) can be downloaded
% at www.lucadesiena.com to test the code.
%%
% INSTRUCTIONS:
% The current release (3.0) works following these steps:
% 
% # Download the package at https://github.com/LucaDeSiena/MuRAT.
% # Download the two sample datasets at https://www.lucadesiena.com/murat.
% # Unzip the MSH (Mount St. Helens) and Romania datasets and put
% the folders in the Murat-master folder.
% # Build your own input file (.m) - each field is described in the
% attached README file and in the example INPUT code. Build your 
% Input files  fom those.
% # Run MuRAT3.m and select the name of the input file you created.
%
% Author:  L. De Siena, November 2020
%
%% INPUTS AND CHECKS
% The code asks for an input file, build one from the sample files provided
% and select it after prompt
addpath('./Utilities_Matlab')

clear; close all; clc

% Call input file
[file,path]                         =   uigetfile('Murat_input*.mlx');

if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end

run(fullfile(path, file));

Murat                               =   Murat_checks(Murat);
save Murat_checks.mat   Murat

%% SEISMIC DATA PROCESSING AND FORWARD MODELLING
load Murat_checks.mat
tic
disp('Data Section')
if Murat.input.parallelized == 1
    Murat                           =   Murat_dataParallelized(Murat);
else
    Murat                           =   Murat_data(Murat);
end
toc
save Murat_forward.mat Murat

%%  2D PEAK-DELAY AND Qc TOMOGRAPHIC INVERSIONS
load Murat_forward.mat

disp('Inversion Section')

Murat                               =   Murat_inversion(Murat);
save Murat_inverse.mat Murat

%% CREATING PLOTS
load Murat_inverse.mat
close all
disp('Plot Section')

Murat                               =   Murat_plot(Murat);

save('Murat.mat','Murat');
