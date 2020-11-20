function M=checkerBoard3D(varargin)

% function M=checkerBoard3D(siz)
% ------------------------------------------------------------------------
% This function creates a checkboard image of the size siz whereby elements
% are either black (0) or white (1). The first element is white.
% example: siz=[12 12 6]; blockSize=2; 
%Block size in pixel units
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2008/08/15
% 2018/12/10 Added checkboard block size input
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        siz=varargin{1};
        blockSize=1;
    case 2        
        siz=varargin{1};        
        blockSize=varargin{2};
end

% if ~isrounded(blockSize) || blockSize<1
%     error('Block size should be positive integer');
% end

%Coping with 1D or 2D input
if numel(siz)==2
    siz(3)=1; 
elseif numel(siz)==1
    siz(2:3)=[1 1];
end

%%

if blockSize==1
    [I,J,K]=ndgrid(1:1:siz(1),1:1:siz(2),1:1:siz(3));
else    
    [I,J,K]=ndgrid(indexRange(blockSize,siz(1)),indexRange(blockSize,siz(2)),indexRange(blockSize,siz(3)));
end
logic_ij=((iseven(I)| iseven(J)) & ((iseven(I)~=iseven(J))));
M=false(siz);
M(iseven(K))=logic_ij(iseven(K));
M(~iseven(K))=~logic_ij(~iseven(K));

end

%%
function i=indexRange(blockSize,s)

i=repmat(1:blockSize:s,blockSize,1); 
i=i(:);
if numel(i)<s
    ii=i;
    i=siz(1)*ones(1,s);
    i(1:num_i)=ii;
else
    i=i(1:s);
end
[~,~,i]=unique(i);
end