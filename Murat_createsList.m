function [listWithFolder,listNoFolder]  =   Murat_createsList(directory)
%CREATES a list of visible files in a folder, outputs both with and
% without folder

% Get general paths/figure options
list                                    =   dir(directory);
list                                    =...
    list(~startsWith({list.name}, '.'));

listWithFolder                          =...
    fullfile({list.folder},{list.name})';

listNoFolder                            =   {list.name}';
