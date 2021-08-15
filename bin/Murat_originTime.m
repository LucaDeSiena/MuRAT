function [theoreticalTime_i, originTime_i] =...
    Murat_originTime(pktime_i,originTime,v,sst_i,SAChdr) %#ok<INUSD>
% function [theoreticalTime_i, originTime_i] =...
%     Murat_originTime(pktime_i,originTime,v,sst_i,SAChdr_i)
%
% CHECKS in case the zero time is missing in the header and CREATES it
%
% Input parameters:
%    pktime_i:          peak time in seconds
%    originTime:        input origin time
%    v:                 average velocity
%    sst_i:             locations of earthquakes and stations
%
% Output parameters:
%    theoreticalTime_i: Theoretical time to include in computations
%    originTime_i:      output origin time

if isequal(originTime,[]) || isequal(eval(originTime),-12345)
    
    Distance                    =...
    sqrt((sst_i(4)-sst_i(1))^2+(sst_i(5)-sst_i(2))^2+...
    (sst_i(6)-sst_i(3))^2);

    theoreticalTime_i           =   Distance/v/1000;
    originTime_i                =   pktime_i-theoreticalTime_i;
    
else
    
    originTime_i                =   eval(originTime);
    theoreticalTime_i           =   pktime_i - originTime_i;
    if theoreticalTime_i < 0
       error(['The picking is set before the origin time for recording '...
            num2str(i_label)]);
    end
    
end
end

