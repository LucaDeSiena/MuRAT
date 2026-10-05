function CovM   =   Murat_smoothing(coordPrior,smoothValue,sigmaM)
% function CovM   =   Murat_smoothing(coordPrior,smoothValue,sigmaM)
%
% CREATES smoothing matrix
%
% Input parameters:
%    coordPrior:        coords for the starting solution
%    smoothValue:       smoothing value
%    sigmaM:            variance for models
% Output parameters:
%    CovM:              smoothing covariance matrix


% Prior
% coordPrior is N-by-3: [lat(deg), lon(deg), z(km)]
R           =   6371;    % Earth mean radius in km

lat         =   deg2rad(coordPrior(:,1));   % N-by-1
lon         =   deg2rad(coordPrior(:,2));   % N-by-1
z           =   coordPrior(:,3)/1000;       % N-by-1 (km)

% pairwise angular differences (N-by-N)
dlat        =   lat - lat.';
dlon        =   lon - lon.';

% haversine for great-circle distance (horizontal)
aa          =   sin(dlat/2).^2 + (cos(lat) .* cos(lat.')).* sin(dlon/2).^2;
cc          =   2 * atan2(sqrt(aa), sqrt(1 - aa));
dh          =   R * cc;  % horizontal surface distance in km (N-by-N)

% vertical (radial) pairwise differences
dz          =   z - z.'; % N-by-N in km

% 3-D straight-line distances
Dkm         =   sqrt(dh.^2 + dz.^2);

% Setting smoothing parameter
if ~isempty(smoothValue)
    smoothValue     =   abs(smoothValue);
else
    v               =   coordPrior(:,3);
    ref             =   v(1);
    d               =   abs(v - ref);
    d(v == ref)     =   Inf;

    idx             =   find(d == min(d));
    nearestVals     =   v(idx);
    nearVal1        =   nearestVals(idx(1));
    smoothValue     =   abs(3*(nearVal1-ref))/1000;
end

% CovM: if N is large this is expensive; build sparse approximation if needed
CovM_full   =   sigmaM * exp(-0.5*(Dkm.^2) / (2*smoothValue^2));
N           =   size(coordPrior,1);

% If N large, consider thresholding to sparse:
if N > 2000
    % keep only neighbors within 3*smoothValue
    mask    =   Dkm <= 3*smoothValue;
    CovM    =   sparse(double(mask) .* CovM_full);
else
    CovM    =   CovM_full;
end