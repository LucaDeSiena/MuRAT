function [t, data, SAChdr] = fget_sac(filename)
% FGET_SAC  Read a SAC binary file into MATLAB.
%
%   [t, data, SAChdr] = fget_sac(filename)
%
%   Reads the SAC (Seismic Analysis Code) binary format and returns the
%   time vector, waveform data, and a header struct whose nested field
%   layout is identical to the original fget_sac / sachdr routines
%   by Zhigang Peng and Xianglei Huang, so all calling code is unchanged.
%
%   Output:
%       t       double vector  (npts x 1)  time in seconds from b
%       data    double vector  (npts x 1)  waveform samples
%       SAChdr  struct         header (see field layout below)
%
%   SAChdr field layout (subsets used by MuRAT):
%       .times.delta  .times.b  .times.e  .times.o  .times.a
%       .times.t0 … .times.t9
%       .station.stla  .stlo  .stel  .stdp  .cmpaz  .cmpinc
%       .station.kstnm
%       .event.evla  .evlo  .evel  .evdp  .mag
%       .event.nzyear .nzjday .nzhour .nzmin .nzsec .nzmsec .kevnm
%       .event.imagtyp  .event.imagsrc
%       .user.data  .user.label
%       .descrip.iftype .idep .iztype .iinst .istreg .ievreg
%       .descrip.ievtyp .iqual .isynth
%       .evsta.dist  .az  .baz  .gcarc
%       .llnl.xminimum .xmaximum .yminimum .ymaximum
%       .llnl.norid .nevid .nxsize .nysize
%       .response  (1x10 double)
%       .data.trcLen  .data.scale
%
%   SAC binary layout (little-endian IEEE 754, 632-byte header):
%       bytes   1 – 280  : 70 float32  (header floats,  arranged [5 x 14])
%       bytes 281 – 440  : 40 int32    (header integers, arranged [5 x  8])
%       bytes 441 – 632  : 192 uint8   (header chars,   arranged [24 x  8])
%       bytes 633 –  end : npts float32 (waveform data)
%
%   Copyright (c) 2026 MuRAT Authors
%   SPDX-License-Identifier: MIT

if nargin < 1
    error('fget_sac:noInput', 'Provide a filename.');
end

fid = fopen(filename, 'rb');
if fid == -1
    error('fget_sac:cannotOpen', 'Cannot open file: %s', filename);
end

% ── Read raw header blocks ────────────────────────────────────────────────
h1  = reshape(fread(fid, 70,  'float32=>double'), 5, 14)';   % [14 x 5]
h2  = reshape(fread(fid, 40,  'int32=>double'),   5,  8)';   % [ 8 x 5]
h3  = reshape(fread(fid, 192, 'uint8'),          24,  8)';   % [ 8 x 24]

npts = h2(2, 5);
data = fread(fid, npts, 'float32=>double');
fclose(fid);

% ── Build output struct ───────────────────────────────────────────────────
SAChdr = buildHeader(h1, h2, h3);

dt = SAChdr.times.delta;
b  = SAChdr.times.b;
t  = (b : dt : b + (npts - 1) * dt)';

end

% =========================================================================
function S = buildHeader(h1, h2, h3)
% Translate raw header blocks into the nested SAChdr struct.

c = @(row, cols) strtrim(char(h3(row, cols)));   % char helper

% .times
S.times.delta = h1(1,1);
S.times.b     = h1(2,1);
S.times.e     = h1(2,2);
S.times.o     = h1(2,3);
S.times.a     = h1(2,4);
S.times.t0    = h1(3,1);
S.times.t1    = h1(3,2);
S.times.t2    = h1(3,3);
S.times.t3    = h1(3,4);
S.times.t4    = h1(3,5);
S.times.t5    = h1(4,1);
S.times.t6    = h1(4,2);
S.times.t7    = h1(4,3);
S.times.t8    = h1(4,4);
S.times.t9    = h1(4,5);
S.times.k0    = c(2,  9:16);
S.times.ka    = c(2, 17:24);
S.times.kt0   = c(3,  1: 8);
S.times.kt1   = c(3,  9:16);
S.times.kt2   = c(3, 17:24);
S.times.kt3   = c(4,  1: 8);
S.times.kt4   = c(4,  9:16);
S.times.kt5   = c(4, 17:24);
S.times.kt6   = c(5,  1: 8);
S.times.kt7   = c(5,  9:16);
S.times.kt8   = c(5, 17:24);
S.times.kt9   = c(6,  1: 8);
S.times.kf    = c(6,  9:16);

% .station
S.station.stla   = h1(7,2);
S.station.stlo   = h1(7,3);
S.station.stel   = h1(7,4);
S.station.stdp   = h1(7,5);
S.station.cmpaz  = h1(12,3);
S.station.cmpinc = h1(12,4);
S.station.kstnm  = c(1, 1:8);
S.stations.kcmpnm  = c(7, 17:24);
S.stations.knetwk  = c(8,  1: 8);

% .event
S.event.evla   = h1(8,1);
S.event.evlo   = h1(8,2);
S.event.evel   = h1(8,3);
S.event.evdp   = h1(8,4);
S.event.mag    = h1(8,5);
S.event.nzyear = h2(1,1);
S.event.nzjday = h2(1,2);
S.event.nzhour = h2(1,3);
S.event.nzmin  = h2(1,4);
S.event.nzsec  = h2(1,5);
S.event.nzmsec = h2(2,1);
S.event.kevnm  = c(1, 9:24);
S.event.imagtyp = [];
S.event.imagsrc = [];

% .user
S.user.data  = [h1(9,1:5), h1(10,1:5)];
S.user.label = [c(6,17:24), c(7,1:8), c(7,9:16)];

% .descrip
S.descrip.iftype  = h2(4,1);
S.descrip.idep    = h2(4,2);
S.descrip.iztype  = h2(4,3);
S.descrip.iinst   = h2(4,5);
S.descrip.istreg  = h2(5,1);
S.descrip.ievreg  = h2(5,2);
S.descrip.ievtyp  = h2(5,3);
S.descrip.iqual   = h2(5,4);
S.descrip.isynth  = h2(5,5);

% .evsta
S.evsta.dist  = h1(11,1);
S.evsta.az    = h1(11,2);
S.evsta.baz   = h1(11,3);
S.evsta.gcarc = h1(11,4);

% .llnl
S.llnl.xminimum = h1(12,5);
S.llnl.xmaximum = h1(13,1);
S.llnl.yminimum = h1(13,2);
S.llnl.ymaximum = h1(13,3);
S.llnl.norid  = [];
S.llnl.nevid  = [];
S.llnl.nxsize = [];
S.llnl.nysize = [];

% .response
S.response = [h1(5,2:5), h1(6,1:5), h1(7,1)];

% .data
S.data.trcLen = h2(2,5);
S.data.scale  = h1(1,4);

end
