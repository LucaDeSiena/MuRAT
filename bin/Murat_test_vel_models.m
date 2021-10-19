function Murat_test_vel_models(orig_model,m_model,interp_model,X,Y,Z,prop_model,Xq,Yq,Zq,origin)

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


% figure
% slice(Xorig_model,Yorig_model,Zorig_model, mVorig_model, 14.13, 40.83,-1800)
% colormap(color)
% view(30,15)
% colorbar
% xlim([14.06 14.20])
% ylim([40.7745 40.8700])
% zlim([-5500 500])
% title('cut original velocity model')

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


% figure
% slice(Xm_model,Ym_model,Zm_model, mVm_model,x_slice,y_slice,z_slice)
% colormap(color)
% view(30,15)
% colorbar
% xlim([0 deg2km(14.13-origin(2))*1000])
% ylim([0 deg2km(40.83-origin(1))*1000])
% zlim([-5000 1000])
% title('cut centered velocity model')

%% now show interpolated model

figure
slice(X,Y,Z,interp_model,x_slice,y_slice,z_slice)
colormap(color)
view(0,90)
colorbar
% xlim([0 deg2km(14.13-origin(2))*1000])
% ylim([0 deg2km(40.83-origin(1))*1000])
% zlim([-5000 1000])
title('inversion velocity model')

%% propagation model

figure
set(gca,'zdir','rev')
slice(Xq,Yq,Zq,prop_model,y_slice,x_slice,z_slice)
colormap(color)
view(0,90)
colorbar
% xlim([0 deg2km(14.13-origin(2))*1000])
% ylim([0 deg2km(40.83-origin(1))*1000])
% zlim([-5000 1000])
title('propagation velocity model')

end