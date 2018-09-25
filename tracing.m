% FUNCTION TRACING: Traces the minimum path between source and receiver
% using pseudo-bending. The initial path is a straight segment. The segment
% mid-point is used for bending. The procedure repeats for a maximum
% number of points, or until the changes are below a certain treshold.

function ma=tracing(ray,gridD,pvel)
    
    maxit=100; %Maximum number of iterations allowed
    maxpoints=2000; %Maximum number of points to bend in one direction
    max3=3*maxpoints; %Three dimensions
    
    dconv=sqrt((gridD(1,2)-gridD(1,1))^2 + (gridD(2,2)-gridD(2,1))^2 +...
        (gridD(3,2)-gridD(3,1))^2);
    tlim=dconv/mean(mean(mean(pvel))); %Time limit
    
    
    %Start bending procedure
    xtemp=[ray(1,1) (ray(1,1)+ray(1,2))/2 ray(1,2)];
    ytemp=[ray(2,1) (ray(2,1)+ray(2,2))/2 ray(2,2)];
    ztemp=[ray(3,1) (ray(3,1)+ray(3,2))/2 ray(3,2)];
    
    n=3;
    v=zeros(n,1);
    
    for i=1:n
        xp=xtemp(i);
        yp=ytemp(i);
        zp=ztemp(i);
        vp=velocity(xp,yp,zp,gridD,pvel);
        v(i,1)=vp;
    end
    
    [ta,tra]=traveling(xtemp,ytemp,ztemp,gridD,v,pvel);
    
    ttsave =0;
    for j=1:maxit
        taa=tra;
        [xtemp,ytemp,ztemp,v]=bending(xtemp,ytemp,ztemp,gridD,v,pvel);
        [ta,tra]=traveling(xtemp,ytemp,ztemp,gridD,v,pvel);
        deltat=taa-tra;
        % check for convergence
        if deltat < tlim
            dtemp=sqrt((xtemp(2)-xtemp(1))^2+(ytemp(2)-ytemp(1))^2 ...
                +(ztemp(2)-ztemp(1))^2);
            if abs(ttsave-tra) < tlim && dtemp < dconv
                break
            else
                ttsave=tra;
                % As the travel time decrease is smaller than travel-time
                % treshold (tlim) and/or the travel time increased, the
                % code doubles the number of segments:
                ntemp=2*n-1;
                if(ntemp > max3)
                    error('Too many points, reduce them or increase ');
                end
                ta(1)=0;
                for i=1:n
                    l=n-i+1;
                    k=2*l-1;
                    xtemp(k)=xtemp(l);
                    ytemp(k)=ytemp(l);
                    ztemp(k)=ztemp(l);
                    ta(k)=ta(l);
                end
                for i=2:2:(2*n-2)
                    xtemp(i)=(xtemp(i-1)+xtemp(i+1))/2;
                    ytemp(i)=(ytemp(i-1)+ytemp(i+1))/2;
                    ztemp(i)=(ztemp(i-1)+ztemp(i+1))/2;
                end
                n=2*n-1;
                for i=1:n
                    xp=xtemp(i);
                    yp=ytemp(i);
                    zp=ztemp(i);
                    vp=velocity(xp,yp,zp,gridD,pvel);
                    v(i)=vp;
                end
            end
        end
    end
    ma(1:n,1:5) = [(1:n)',xtemp',ytemp',ztemp',ta']; %Final ray
end