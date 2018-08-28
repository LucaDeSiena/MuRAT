function [ip,jp,kp]=corner(xx,yy,zz,ngrid,grid)

% Finds the indices of the "minimum" coordinates of the 
% cube of grid points surrounding the point (xx,yy,zz).
% ip = index of minimum x grid coord (min. x coord. = grid(1,ip)).
% jp = index of minimum y grid coord (min. y coord. = grid(2,jp)).
% kp = index of minimum z grid coord (min. z coord. = grid(3,kp)).
x(1)=xx;
x(2)=yy;
x(3)=zz;

%corner
if (xx < grid(1,1)) || (xx > grid(1,ngrid(1))) || ...
    (yy < grid(2,1)) || (yy > grid(2,ngrid(2))) || ...
    (zz < grid(3,1)) || (zz > grid(3,ngrid(3)))
    no = [xx,yy,zz];
    display(no);
    error('RAY BEND ERR: point is outside grid, called by corner');
end

% Loop over dimensions
indx=zeros(3,1);
for j=1:3
    indx(j)=ngrid(j)-1;
    while x(j)<grid(j,indx(j))
        indx(j)=indx(j)-1;
    end
end
ip=indx(1);
jp=indx(2);
kp=indx(3);
