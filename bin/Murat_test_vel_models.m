function Murat_test_vel_models(orig_model,m_model,interp_model,X,Y,Z,prop_model,Xq,Yq,Zq)

% first look at original model
% use min and maximum lat / lon of model
lat_res = unique(orig_model(:,1));
lon_res = unique(orig_model(:,2));
depth_res = unique(orig_model(:,3));

%3D
[Xorig_model,Yorig_model,Zorig_model] =...
    meshgrid(lon_res,lat_res,depth_res);
mVorig_model   =...
    griddata(orig_model(:,2),orig_model(:,1),orig_model(:,3),...
    orig_model(:,4),Xorig_model,Yorig_model,Zorig_model);

[color]=inferno(100);

figure
slice(Xorig_model,Yorig_model,Zorig_model, mVorig_model, 14.13, 40.83,-1800)
colormap(color)
view(0,90)
colorbar
title('original velocity model')


%% check out centered model
% use min and maximum lon / lat of model
lon_res2 = unique(m_model(:,1));
lat_res2 = unique(m_model(:,2));
depth_res2 = unique(m_model(:,3));

%3D
[Xm_model,Ym_model,Zm_model] = meshgrid(lon_res2,lat_res2,depth_res2);
mVm_model   =...
    griddata(m_model(:,1),m_model(:,2),m_model(:,3),...
    m_model(:,4),Xm_model,Ym_model,Zm_model);


x_slice = 6000;
y_slice = 8000;
z_slice = -1800;

figure
slice(Xm_model,Ym_model,Zm_model, mVm_model,x_slice,y_slice,z_slice)
colormap(color)
view(0,90)
colorbar
title('cartesian velocity model')

%% now show interpolated model

figure
slice(X,Y,Z,interp_model,x_slice,y_slice,z_slice)
colormap(color)
view(0,90)
colorbar
title('inversion velocity model')

%% propagation model

figure
set(gca,'zdir','rev')
slice(Xq,Yq,Zq,prop_model,y_slice,x_slice,z_slice)
colormap(color)
view(0,90)
colorbar
title('propagation velocity model')

end