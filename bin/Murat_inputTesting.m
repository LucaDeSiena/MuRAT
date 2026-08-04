function [checkInput,spikeInput]	=...
    Murat_inputTesting(I,spike_o,spike_e,x,y,z)
% function [checkInput,spikeInput]	=...
%     Murat_inputTesting(I,spike_o,spike_e,x,y,z)
%
% CREATES checkerboard and spike inputs.
%   
% Input parameters:
%    I:                 Input matrix
%    spike_o:           Marker to see if you define the spike
%    spike_e:           Location of the spike
%    x:                 x vector
%    y:                 y vector
%    z:                 z vector
%
% Output parameters:
%    checkInput:        input matrix for the checkerboard
%    spikeInput:        input matrix for the spike

checkInput      = reshape(permute(I, [3 2 1]), [], 1);

if ~isempty(spike_o)
    r           =   Murat_unfoldXYZ(x,y,z);
    spikeInput  =   r(:,1)>spike_o(2) & r(:,1)<spike_e(2) &...
                    r(:,2)>spike_o(1) & r(:,2)<spike_e(1) &...
                    r(:,3)<spike_o(3) & r(:,3)>spike_e(3);
else
    spikeInput  =   [];
end

end