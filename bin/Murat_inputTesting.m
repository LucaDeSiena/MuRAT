function [checkInput,spikeInput]	=...
    Murat_inputTesting(I,spike_o,spike_e,x,y,z)

% We fold in the checkerboard velocity model and assign alternating values
% to the checkerboard.
    
[nxc,nyc,nzc]                       =   size(I);
index                               =   0;
checkInput                          =   zeros(nxc*nyc*nzc,1);

for i=1:nxc
    for j=1:nyc
        for l=1:nzc
            index                   = index+1;
            checkInput(index)       = I(i,j,l);
        end
    end
end

if ~isempty(spike_o)
        r                           =   Murat_unfold(x',y',z');
        spikeInput                  =...
            r(:,1)>spike_o(2) & r(:,1)<spike_e(2) &...
            r(:,2)>spike_o(1) & r(:,2)<spike_e(1) &...
            r(:,3)<spike_o(3) & r(:,3)>spike_e(3);
end

end