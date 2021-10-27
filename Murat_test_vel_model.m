function Murat_test_vel_model(titelname,model,sections,origin, X,Y,Z)

% check the velocity model
% mode          = velocity model
% sections      = as defined in input
% origin        = as defined in input

if nargin == 3
    % the original velociy model is in lat, lon, depth
    lat_res         = unique(model(:,1));
    lon_res         = unique(model(:,2));
    depth_res       = unique(model(:,3));
    
    %3D
    [X,Y,Z]         = meshgrid(lon_res,lat_res,depth_res);
    
    mVgrid          = griddata(model(:,2),model(:,1),model(:,3),...
        model(:,4),X,Y,Z);
    clear model
    model           = mVgrid;
    
    x_sec = sections(2);
    y_sec = sections(1);
    
elseif nargin == 4
    % the centered velociy model is in lon, lat, depth
    lat_res         = unique(model(:,2));
    lon_res         = unique(model(:,1));
    depth_res       = unique(model(:,3));
    
    %3D
    [X,Y,Z]         = meshgrid(lon_res,lat_res,depth_res);
    
    mVgrid          = griddata(model(:,1),model(:,2),model(:,3),...
        model(:,4),X,Y,Z);
    clear model
    model           = mVgrid;
    
end

if nargin > 3
    x_sec = deg2km(sections(2)-origin(2))*1000;
    y_sec = deg2km(sections(1)-origin(1))*1000;
end


[color]=inferno(100);

figure
slice(X,Y,Z, model, x_sec, y_sec,sections(3))
colormap(color)
view(10,20)
colorbar
title(titelname)


end