function Murat_test_vel_models(orig_model,m_model,interp_model,X,Y,Z,prop_model,Xq,Yq,Zq,origin)

% first look at original model
% use min and maximum lat / lon of model
lat_res = unique(orig_model(:,1));
lon_res = unique(orig_model(:,2));
depth_res = unique(orig_model(:,3));

%3D
[Xorig_model,Yorig_model,Zorig_model] = meshgrid(lon_res,lat_res,depth_res);
mVorig_model   =   griddata(orig_model(:,2),orig_model(:,1),orig_model(:,3),orig_model(:,5),Xorig_model,Yorig_model,Zorig_model);

[color]=inferno(100);

figure
slice(Xorig_model,Yorig_model,Zorig_model, mVorig_model, 36.0, -2.6,-5)
colormap(color)
view(30,15)
colorbar
title('entire original velocity model')


figure
slice(Xorig_model,Yorig_model,Zorig_model, mVorig_model, 36.0, -2.6,-5)
colormap(color)
view(30,15)
colorbar
xlim([35.7 36.3])
ylim([-3.0 -2.4])
zlim([-30 10])
title('cut original velocity model')

%% check out centered model
% use min and maximum lat / lon of model
lat_res2 = unique(m_model(:,1));
lon_res2 = unique(m_model(:,2));
depth_res2 = unique(m_model(:,3));

%3D
[Xm_model,Ym_model,Zm_model] = meshgrid(lon_res2,lat_res2,depth_res2);
mVm_model   =   griddata(m_model(:,2),m_model(:,1),m_model(:,3),m_model(:,4),Xm_model,Ym_model,Zm_model);


x_slice = deg2km(36.0-origin(2))*1000;
y_slice = deg2km(-2.6-origin(1))*1000;
z_slice = -5*1000;

figure
slice(Xm_model,Ym_model,Zm_model, mVm_model,x_slice,y_slice,z_slice)
colormap(color)
view(30,15)
colorbar
title('centered velocity model')


figure
slice(Xm_model,Ym_model,Zm_model, mVm_model,x_slice,y_slice,z_slice)
colormap(color)
view(30,15)
colorbar
xlim([0 deg2km(36.3-origin(2))*1000])
ylim([0 deg2km(-2.4-origin(1))*1000])
zlim([-30000 10000])
title('cut centered velocity model')

%% now show interpolated model

figure
slice(X,Y,Z, interp_model, x_slice,y_slice,z_slice)
colormap(color)
view(30,15)
colorbar
xlim([0 deg2km(36.3-origin(2))*1000])
ylim([0 deg2km(-2.4-origin(1))*1000])
zlim([-30000 10000])
title('reduced velocity model')

%% propagation model

figure
slice(Xq,Yq,Zq, prop_model, x_slice,y_slice,z_slice)
colormap(color)
view(30,15)
colorbar
xlim([0 deg2km(36.3-origin(2))*1000])
ylim([0 deg2km(-2.4-origin(1))*1000])
zlim([-30000 10000])
title('reduced interpreted velocity model')

end