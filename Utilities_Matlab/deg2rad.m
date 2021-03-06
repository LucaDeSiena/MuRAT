function rad=deg2rad(deg)
% DEG2RAD  Converts decimal degrees to radians. Vectorized.
% Version: 18 Jan 96
% Usage:   rad=deg2rad(deg)
% Input:   deg - vector of angles in decimal degrees
% Output:  rad - vector of angles in radians

% Copyright (c) 2011, Michael R. Craymer
% All rights reserved.
% Email: mike@craymer.com

rad=deg.*pi./180;
