function M = checkerBoard3D(siz, blockSize)
% CHECKERBOARD3D  Create a 3-D logical checkerboard array.
%
%   M = checkerBoard3D(siz)
%   M = checkerBoard3D(siz, blockSize)
%
%   Returns a logical array of size siz (1-, 2-, or 3-element vector) in
%   which adjacent voxels always have opposite values.  The element at
%   index (1,1,1) is true (white).  With blockSize > 1 (positive integer),
%   each "square" of the checkerboard spans blockSize voxels in every
%   direction.
%
%   This function is a drop-in replacement for the identically named
%   routine from the GIBBON toolbox (Kevin Mattheus Moerman, GPL v3).
%   The algorithm is independent and is released under the MIT licence.
%
%   Copyright (c) 2026 MuRAT Authors
%   SPDX-License-Identifier: MIT
%
%   Permission is hereby granted, free of charge, to any person obtaining
%   a copy of this software and associated documentation files (the
%   "Software"), to deal in the Software without restriction, including
%   without limitation the rights to use, copy, modify, merge, publish,
%   distribute, sublicense, and/or sell copies of the Software, and to
%   permit persons to whom the Software is furnished to do so, subject to
%   the following conditions: The above copyright notice and this
%   permission notice shall be included in all copies or substantial
%   portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT
%   WARRANTY OF ANY KIND.

if nargin < 2, blockSize = 1; end

% Pad siz to length 3 (handle 1-D and 2-D callers)
siz = siz(:).';
siz(end+1 : 3) = 1;

% Map each voxel coordinate to its block index (1-based)
%   ceil(coord / blockSize) gives the block number for that voxel.
[I, J, K] = ndgrid( ceil((1:siz(1)) / blockSize), ...
                     ceil((1:siz(2)) / blockSize), ...
                     ceil((1:siz(3)) / blockSize) );

% A cell is "white" when the parity sum of block indices is even,
% matching the convention that (1,1,1) is true:
%   mod(1+1+1,2)=1 → odd → flip once → true  ✓
M = mod(I + J + K, 2) == 1;

end
