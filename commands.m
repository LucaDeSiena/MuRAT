%COMMANDS for MuRAT

% Choose between P- (2) and S-direct (3) waves
PorS = 2;

% Work with 1 vertical (1) or 2 horizontal (2) regcordings, or with the
% three components(3)
compon = 1;

% Step of the inversion - set to 2 in case of double inversion, 1 otherwise
nstep= 1;

% Name of the file containing the 1st velocity model - according to PorS
modv = load('modv.txt','%f');
v0= mean(modv(:,4));
resol2 = abs(modv(2,3)-modv(1,3))/2;

% Name of the file containing the picking of origin time, P, and S waves
tempi = load('tempi.txt','%f');

% Length of the window used to measure P- or S- wave energy in
% seconds
fin = 1;

% Central frequency
cf = 18;

% Lapse time (from the origin time of the event) must be >> 2*tS
tC = 15;

% Total coda window length
tW = 10;

%UTM coordinates of the origin of the velocity model
originWE=538311;
originSN=5092338;
originz=3350;

%Rays are measured in meters. If the ray length is already in km set um =1;
um = 1000;

%In order to check the images, one can change the smoothing parameter
%(set smoot=1)
smoot=0;

%The minimum coda-to-noise energy ratio for the weighted inversion
soilnoise = 1.4;

if nstep == 2
%     Choose to select just part of the data affected by anomalous
%     behavior of the energy ratios with travel times
%    tA1=0;
%    tA2=2;
%     Otherwise impose a second grid in the regions showing this behavior 
     modv2 = load('modv2.txt','%f');
     resol22 = abs(modv2(2,3)-modv2(1,3))/2;
     smoot2=0;
end
