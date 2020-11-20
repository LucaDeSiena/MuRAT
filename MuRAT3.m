%MURAT   Multi-Resolution (seismic) Attenuation Tomography
%   SCOPE: 3D total attenuation, scattering, and absorption tomography
%
%   Author:  L. De Siena, November 2020
% 
% This is v3.0 of the program.
% 
% Installation
% ------------
% SYSTEM: The program works on all Mac and Linux computers where it has
% been tried.
%
% MATLAB Version: The code keeps the us of toolboxes to a minimum. However,
% there are four necessary toolboxes: Signal Processing, Curve Fitting,
% Image Processing and Mapping Toolboxes. The parallel computing toolbox is
% highly recommended for speed.
%
% Two sample datasets (Mount St. Helens and Romania) can be downloaded
% at www.lucadesiena.com to test the code.

% INSTRUCTIONS:
% The current release (3.0) works following these steps:
%
% 1. Download the package at https://github.com/LucaDeSiena/MuRAT.
%
% 2. Download the two sample datasets at
% https://doi.pangaea.de/10.1594/PANGAEA.893893.
% Unzip the MSH (Mount St. Helens) and Romania datasets and put
% the folders in the Murat-master folder.
%
% 3. Build your own input file (.m) - each field is described in the
% attached README file and in the INPUT sections of the code.
% When building your example, use one of the Input files as template
% (MSH and Romania.
%
% 4. Run MuRAT3.m and select the name of the input file desired.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS AND CHECKS

addpath('./Utilities_Matlab')

clear; close all; clc

disp('Checks and Loops')

[file,path]                         =   uigetfile('Murat_input*.m');

if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end

run(fullfile(path, file)) 

Murat                               =   Murat_checks(Murat);
save Murat_checks.mat   Murat

%% SEISMIC DATA PROCESSING AND FORWARD MODELLING
load Murat_checks.mat

disp('Data Section')
tic
Murat                               =   Murat_data(Murat);
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
