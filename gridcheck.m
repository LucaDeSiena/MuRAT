function [grid,ngrid,pvel]=gridcheck(passox,passoy,passoz,modvPS,PorS)

% Builds the grid for ray-tracing and checks for zeros

%Necessary parameters!
%step of the grid
dx=1000;
dy=1000;
dz=1000;

ix=passox/dx;%numer of x layers,given the step of the grid
iy=passoy/dy;%numer of y layers,given the step of the grid
iz=passoz/dz;%numer of depths layers,given the step of the grid

blkparx=ix+11;%do not change
blkpary=iy+11;%do not change
blkparz=iz+11; %do not change

% grid dimensions
pvel = zeros(blkparx,blkpary,blkparz);

if PorS==2
    mPorS=4;
elseif PorS==3
    mPorS=5;
end

% Interpolation avoids points outside grid

%NUMBER OF X, Y, AND Z LAYERS

for k=1:iz
    la=find(modvPS(:,3)==modvPS(k,3));
    file=[modvPS(:,1:3) modvPS(:,mPorS)];
    index=0;
    for i=1:ix
        for j=1:iy
            index=index+1;
            iii1=file(index,1);
            iii2=file(index,2);
            pvel(i+5,j+5,k+5) = file(index,4);
        end
    end
%  pad array in x
    for j=1:iy
        for i = 1:5
            pvel(i,j+5,k+5) = pvel(6,j+5,k+5);
        end
        for i = ix+6:ix+10
            pvel(i,j+5,k+5) = pvel(ix+5,j+5,k+5);
        end
    end
end
%  set new ix
ix = ix + 10;
for i = 1:ix
    for j = 1:iy
        for k = 1:iz
            if pvel(i,j+5,k+5)==0
                disp('after x interp pvel = 0 for ', i,j,k);
            end
        end
    end
end
disp('completed check for zeros after x interpolation');
%  pad ends of array in z
for k = 1:5
    for j = 1:iy
        for i = 1:ix
            pvel(i,j+5,k) = pvel(i,j+5,6);
        end
    end
end
izend = iz + 10;
izst = iz + 6;
for k = izst:izend
    for j = 1:iy
        for i = 1:ix
           pvel(i,j+5,k) = pvel(i,j+5,iz);
        end
    end
end
iz = izend;
for i = 1:ix
    for j = 1:iy
        for k = 1:iz
            if(pvel(i,j+5,k)==0)
                disp('after y interp pvel = 0 for i,j,k = ', i,j,k);
            end
        end
    end
end
disp('completed check for zeros after z interpolation');
%  pad in y
for j = 1:5
    for k = 1:iz
        for i = 1:ix
           pvel(i,j,k) = pvel(i,6,k);
        end
    end
end
iyend = iy + 10;
iyst = iy + 6;
for j = iyst:iyend
    for k = 1:iz
        for i = 1:ix
           pvel(i,j,k) = pvel(i,iy+5,k);
        end
    end
end
iy = iyend;
%  check for zeros in array
for i = 1:ix
    for j = 1:iy
        for k = 1:iz
            if(pvel(i,j,k)==0)
                disp('after y interp pvel = 0 for i,j,k = ', i,j,k)
            end
        end
    end
end
disp('completed check for zeros after y interpolation');

%   to check if the program reads the velocity model well
% prova=zeros(ix*iy,5);
% index = 0;
% for i = 1:ix
%     for j=1:iy
%         index=index+1;
%         prova(index,1:5)= [i,j,pvel(i,j,6),iii1,iii2];
%     end
% end

% Set parameters for subroutines
ngrid(1) = ix;
ngrid(2) = iy;
ngrid(3) = iz;
% use a costant-step grid
grid = zeros(3,ix+5);
for i = 1:ix
    grid(1,i+5) = (i-1)*dx;
end
for j = 1:iy
    grid(2,j+5) = (j-1)*dy;
end
for k = 1:iz
    grid(3,k+5) = (k-1)*dz;
end

%  set grid locations for regions around model
for i = 1:5
    grid(1,i) = grid(1,6) - 200*(6-i)*dx;
    grid(2,i) = grid(2,6) - 200*(6-i)*dy;
    grid(3,i) = grid(3,6) - 200*(6-i)*dz;
    grid(1,ix-i+1) = grid(1,ix-6) +  200*(6-i)*dx;
    grid(2,iy-i+1) = grid(2,iy-6) +  200*(6-i)*dy;
    grid(3,iz+5-i+1) = grid(3,iz-6) +  200*(6-i)*dz;
end
