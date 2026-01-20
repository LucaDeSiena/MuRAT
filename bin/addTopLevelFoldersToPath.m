function addTopLevelFoldersToPath(folderName)
% Adds folderName and its first-level subfolders to MATLAB path (no deep recursion)
if ~exist(folderName,'dir')
    return;
end
addpath(folderName);
d = dir(folderName);
isub = [d(:).isdir];
names = {d(isub).name};
names = names(~ismember(names,{'.','..'}));
for k = 1:numel(names)
    addpath(fullfile(folderName,names{k}));
end
end
