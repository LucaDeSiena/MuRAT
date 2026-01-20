function modvOriginal     =   Murat_readLithos1(filename,queryLat,queryLon)
% modvOriginal     =   Murat_readLithos1(filename,queryLat,queryLon)
%
% READS Lithos1.0 data and builds the 1D model in iaspei91 style, given
% average latitude and longitude
%
% Input parameters:
%    filename:          name of the velocity model (LITHOS1.0.nc)
%    queryLat:          average latitude of the model
%    queryLon:          average longitude of the model
%
% Output parameters:
%    modvOriginal:      original velocity model in iaspei91 format

info        =   ncinfo(filename);
varnames    =   {info.Variables.Name};

lats        =   ncread(filename, "latitude");
lons        =   ncread(filename, "longitude");

% Find nearest indices (assumes lats and lons are 1D)
[~, iLat]   =   min(abs(lats - queryLat));
[~, iLon]   =   min(abs(lons - queryLon));

% Identify depth variables (names that end with '_depth')
isDepth     =   endsWith(varnames, "_depth");
depthVars   =   varnames(isDepth);

% Prepare storage
depthValues =   zeros(numel(depthVars),1);
baseNames   =   strings(numel(depthVars),1);

for k = 1:numel(depthVars)
    dname   =   depthVars{k};
    % read depth matrix and extract scalar at grid cell
    try
        dmat=   ncread(filename, dname);
    catch
        error("Failed to read variable %s", dname);
    end
    % dmat may be [lon x lat] or [lat x lon]; choose by matching sizes
    % Try both index orders safely:
    val = tryIndex(dmat, iLon, iLat, iLat, iLon);
    depthValues(k) = val;
    base = extractBefore(dname, "_depth");
    baseNames(k) = base;
end

% Collect properties for each base (other vars that start with base + "_")
propsMap = containers.Map; % key: property name (e.g., 'vp'), value: vector across layers
% also collect depths associated per layer
for k = 1:numel(baseNames)
    base = baseNames(k);
    depth_k = depthValues(k);
    % find variables that start with base + "_" but are not the depth var
    prefix = base + "_";
    matches = startsWith(varnames, prefix);
    matches = matches & ~isDepth; % exclude the depth var itself
    matchedNames = varnames(matches);
    for m = 1:numel(matchedNames)
        vname = matchedNames{m};
        % property name = suffix after base + '_'
        prop = extractAfter(vname, char(prefix));
        % read variable and extract scalar at grid cell
        vmat = ncread(filename, vname);
        vval = tryIndex(vmat, iLon, iLat, iLat, iLon);
        % append (depth, value) pair into a per-property cell array
        if ~isKey(propsMap, prop)
            propsMap(prop) = [depth_k, vval]; % first row: [depth, value]
        else
            propsMap(prop) = [propsMap(prop); depth_k, vval];
        end
    end
end

% Now for each property sort by depth and create column vectors
propNames = keys(propsMap);
result = struct();
for p = 1:numel(propNames)
    prop = propNames{p};
    pairs = propsMap(prop); % N x 2: [depth, value]
    % remove rows with NaN depth or value
    pairs = pairs(~isnan(pairs(:,1)), :);
    % sort by depth ascending (shallow -> deep)
    [~, idx] = sort(pairs(:,1), 'ascend');
    pairs = pairs(idx, :);
    result.(prop).depth = pairs(:,1);
    result.(prop).value = pairs(:,2);
end

modvp   =   [result.vp.depth 6371-result.vp.depth result.vp.value result.vs.value];
modvp   =   sortrows(modvp,[1,3]);
[g, groupVals] = findgroups(modvp(:,1));
meanCol3 = splitapply(@mean, modvp(:,3), g);
meanCol4 = splitapply(@mean, modvp(:,4), g);

idx = find(meanCol3 == 0 & (1:numel(meanCol3))' < numel(meanCol3));
meanCol3(idx) = meanCol3(idx + 1)/2;

idx = find(meanCol4 == 0 & (1:numel(meanCol4))' < numel(meanCol4));
meanCol4(idx) = meanCol4(idx + 1)/2;

modvp1 = [-5 6376 meanCol3(1) meanCol4(1); groupVals, 6371-groupVals, meanCol3, meanCol4];

xq = min(modvp1(:,1)):1:max(modvp1(:,1));

y3q = interp1(modvp1(:,1), modvp1(:,3), xq, 'pchip', 'extrap');
y4q = interp1(modvp1(:,1), modvp1(:,4), xq, 'pchip', 'extrap');

modvOriginal = [xq' 6371-xq' y3q' y4q'];
