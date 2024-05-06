% make an ensemble plot of pv and w for the vortices in the r=0.01
% simulation

close all; clear;

for ii = 1:3
    
    if ii ==1

% high rossby
path_root = '/nobackup1c/users/mkohl/linear_gradient_formulation/primitive_equations_constant_strat_r001_eps04_final_high_temp_output_with_advection_more_relaxation';
eps = 0.4;
    elseif ii ==2
% mid rossby
path_root = '/nobackup1c/users/mkohl/linear_gradient_formulation/primitive_equations_constant_strat_r001_eps01_final_high_temp_output_with_advection_more_relaxation';
eps = 0.1;   
    else
% low rossby
path_root = '/nobackup1c/users/mkohl/linear_gradient_formulation/primitive_equations_constant_strat_r001_eps001_final_high_temp_output_with_advection';
eps = 0.01;
    end
    
path = [path_root,'/snapshots/snapshots_s1.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nz = h5_dims(4); nt = h5_dims(1);

x = linspace(0,6*pi,nx); [X,Y] = meshgrid(x,x);
y = x;
t = h5read(path,'/scales/sim_time');
dx = x(2)-x(1);
z = linspace(0,1,nz);
dz = z(2)-z(1);


% define the time-tendencies of PV
Q = h5read(path,'/tasks/Q');

% calculate Qt from finite differencing

%T = permute(repmat(t,1,nz,nx,nx),[2,3,4,1]);

% % % backwards differencing
% Qt = (Q(:,:,:,2:end)-Q(:,:,:,1:end-1))./(T(:,:,:,2:end)-T(:,:,:,1:end-1));
% Qt = cat(4,zeros(nz,nx,nx),Qt);

% output Qt directly

%Qt = h5read(path,'/tasks/Qt');

% define size of ensemble mean

indx = 13; indy = indx;

w_mean = zeros(nz,2*indx+1,2*indx+1);
Q_mean = zeros(nz,2*indx+1,2*indx+1);
Q_full_mean = zeros(nz,2*indx+1,2*indx+1);
Qt_mean = zeros(nz,2*indx+1,2*indx+1);

%Qx_mean = zeros(nz,2*indx+1,2*indx+1);
%Qy_mean = zeros(nz,2*indx+1,2*indx+1);
%Qz_mean = zeros(nz,2*indx+1,2*indx+1);

Q_adv_x_mean = zeros(nz,2*indx+1,2*indx+1);
Q_adv_y_mean = zeros(nz,2*indx+1,2*indx+1);
Q_adv_z_mean = zeros(nz,2*indx+1,2*indx+1);

%Q_diab_x_mean = zeros(nz,2*indx+1,2*indx+1);
%Q_diab_y_mean = zeros(nz,2*indx+1,2*indx+1);
Q_diab_z_mean = zeros(nz,2*indx+1,2*indx+1);
Q_diab_hoskins_mean = zeros(nz,2*indx+1,2*indx+1);

%theta_mean = zeros(nz,2*indx+1,2*indx+1);
theta_prime_ensemble = zeros(nz,2*indx+1,2*indx+1);
theta_dot_mean = zeros(nz,2*indx+1,2*indx+1);
theta_z_mean = zeros(nz,2*indx+1,2*indx+1);
theta_dot_z_mean = zeros(nz,2*indx+1,2*indx+1);

% reduction factor

r = 0.01;

% variable to composite on: '1'='w' and '2'='q2'

nvar = 1;

% number of vortices to average over
n = 10;

% smooth on or off
smooth = 0;
smethod = 'movmean';

% smooth window size
nwindow = 7;

counter = 0;

for kk = 1:1
    
%path = [pwd,'/snapshots/snapshots_s',num2str(kk),'.h5'];

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nt = h5_dims(1);

for tt = 1:1:nt
    
t_start = 70;
    
if t(tt)<t_start
    continue
end
    
if nvar==1
start_track = [5 1 1 tt];
count_track = [1 nx nx 1];    
z = h5read(path,'/tasks/w',start_track,count_track);
z = squeeze(z);
else
start_track = [3 1 1 tt];
count_track = [1 nx nx 1]; 
z = h5read(path,'/tasks/Q',start_track,count_track);    
z = squeeze(z);
z = z-mean(z,2); % PV anomalies
%z = -z;
end

start = [1 1 1 tt];
count = [nz nx nx 1];

w = h5read(path,'/tasks/w',start,count);
Q = h5read(path,'/tasks/Q',start,count);

Qz = zeros(nz,nx,nx);
Qz(1,:,:) = (Q(2,:,:)-Q(1,:,:))/dz;
Qz(end,:,:) = (Q(end,:,:)-Q(end-1,:,:))/dz;
Qz(2:end-1,:,:) = (Q(3:end,:,:)-Q(1:end-2,:,:))/(2*dz);

%Qx = h5read(path,'/tasks/Qx',start,count);
%Qy = h5read(path,'/tasks/Qy',start,count);
%Qz = h5read(path,'/tasks/Qz',start,count);

Q_adv_x = h5read(path,'/tasks/Q_adv_x',start,count);
Q_adv_y = h5read(path,'/tasks/Q_adv_y',start,count);
Q_adv_z = h5read(path,'/tasks/Q_adv_z',start,count);

Q_adv_x = -Q_adv_x;
Q_adv_y = -Q_adv_y;
Q_adv_z = -eps*Q_adv_z;
%Q_adv_z =  -eps*w.*Qz;

%Q_diab_x = h5read(path,'/tasks/Qdot_x',start,count);
%Q_diab_y = h5read(path,'/tasks/Qdot_y',start,count);
Q_diab_z = h5read(path,'/tasks/Q_dot_z',start,count);


theta = h5read(path,'/tasks/theta',start,count);
theta_prime = theta - mean(theta,2);
theta_z = h5read(path,'/tasks/theta_z',start,count);
theta_dot = h5read(path,'/tasks/theta_dot',start,count);
theta_dot_z = h5read(path,'/tasks/theta_dot_z',start,count);


%Q_diab_hoskins = Q_diab_z - theta_dot.*Qz;
Q_diab_hoskins = eps*Q.*theta_dot_z./(1+eps*theta_z) - eps*theta_dot.*Qz./(1+eps*theta_z);
%Q_diab_hoskins = eps*Q.*theta_dot_z - eps*theta_dot.*Qz;

% define the anomaly with respect to the zonal mean

Q_full = Q;
Q = Q-mean(Q,3);

% track maxima

% apply some smoothing

if smooth ==1

z = smoothdata(smoothdata(z,1,smethod,nwindow),2,smethod,nwindow);

end

% Find dimensions to set up loop
xdim = size(z,1);
ydim = size(z,2);

% Loop through x dimension to find peaks of each row
xpeaks = zeros(size(z));
xwidths = NaN(size(z));
for i = 1:xdim
    [~,locs,width] = findpeaks(z(i,:),'MinPeakHeight',0,'WidthReference','halfheight');
    xpeaks(i,locs) = 1;
    xwidths(i,locs) = width;
end

% Loop through y dimension to find peaks of each row
ypeaks = zeros(size(z));
ywidths = NaN(size(z));
for i = 1:ydim
    [~,locs,width] = findpeaks(z(:,i),'MinPeakHeight',0,'WidthReference','halfheight');
    ypeaks(locs,i) = 1;
    ywidths(locs,i) = width;
end

% Find indices that were peaks in both x and y
peak_inds = xpeaks+ypeaks == 2;

% Save data to sruct
peakdata = struct;
peakdata.peakZ = z(peak_inds);
peakdata.peakX = X(peak_inds);
peakdata.peakY = Y(peak_inds);
peakdata.peakXWidth = xwidths(peak_inds);
peakdata.peakYWidth = ywidths(peak_inds);

% sort into highest peaks

[peak_height,I] = sort(peakdata.peakZ,'descend');
[peak_x] = peakdata.peakX(I);
[peak_y] = peakdata.peakY(I);
[peak_xwidth] = peakdata.peakXWidth(I);
[peak_ywidth] = peakdata.peakYWidth(I);


% Average over the nth biggest vortices

peak_x = peak_x(1:n); peak_y = peak_y(1:n);


% test: plot z and mark the peaks identified

% figure
% z = z/max(max(z));
% contourf(X,Y,z); hold on;
% caxis([-max(max(z)) max(max(z))])
% colorbar
% colormap(redblue)
% scatter(peak_x,peak_y,'+g');

% make an ensemble plot 
%counter = 0;
for i=1:length(peak_x)
 [a,indcx] = min(abs(x-peak_x(i)));
 [a,indcy] = min(abs(y-peak_y(i)));
 
 % exclude points that lie on the boundaries
 if indcx+indx>=nx || indcx-indx<=0 || indcy+indy>=nx || indcy-indy<=0
     continue
 end
 
 Q_full_mean = Q_full_mean + Q_full(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 Q_mean = Q_mean + Q(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 %Qx_mean = Qx_mean + Qx(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 %Qy_mean = Qy_mean + Qy(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 %Qz_mean = Qz_mean + Qz(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 
 w_mean = w_mean + w(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 %Qt_mean = Qt_mean + Qt(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx,tt);
 Q_adv_x_mean = Q_adv_x_mean + Q_adv_x(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 Q_adv_y_mean = Q_adv_y_mean + Q_adv_y(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 Q_adv_z_mean = Q_adv_z_mean + Q_adv_z(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 
 
 %Q_diab_x_mean = Q_diab_x_mean + Q_diab_x(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 %Q_diab_y_mean = Q_diab_y_mean + Q_diab_y(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 Q_diab_z_mean = Q_diab_z_mean + Q_diab_z(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 Q_diab_hoskins_mean = Q_diab_hoskins_mean + Q_diab_hoskins(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 

 
 %theta_mean = theta_mean + theta(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 theta_prime_ensemble = theta_prime_ensemble + theta_prime(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 
 theta_z_mean = theta_z_mean + theta_z(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 theta_dot_mean = theta_dot_mean + theta_dot(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 theta_dot_z_mean = theta_dot_z_mean + theta_dot_z(:,indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 
 
 
 counter = counter + 1;
end



end

end

% take average 

Q_mean = Q_mean/counter;
Q_full_mean = Q_full_mean/counter;
%Qx_mean = Qx_mean/counter;
%Qy_mean = Qy_mean/counter;
%Qz_mean = Qz_mean/counter;
w_mean = w_mean/counter;

%Qt_mean = Qt_mean/counter;

Q_adv_x_mean = Q_adv_x_mean/counter;
Q_adv_y_mean = Q_adv_y_mean/counter;
Q_adv_z_mean = Q_adv_z_mean/counter;

%Q_diab_x_mean = Q_diab_x_mean/counter;
%Q_diab_y_mean = Q_diab_y_mean/counter;
Q_diab_z_mean = Q_diab_z_mean/counter;
Q_diab_hoskins_mean = Q_diab_hoskins_mean/counter;

%theta_mean = theta_mean/counter;
theta_prime_ensemble = theta_prime_ensemble/counter;
theta_z_mean = theta_z_mean/counter;
theta_dot_mean = theta_dot_mean/counter;
theta_dot_z_mean = theta_dot_z_mean/counter;

% define residual

%residual = Qt_mean - Q_adv_x_mean - Q_adv_y_mean - Q_adv_z_mean - Q_diab_x_mean ...
%    - Q_diab_y_mean - Q_diab_z_mean;

% make a plot of the pv ensemble
xc = -indx*dx:dx:indx*dx;
yc = -indy*dx:dx:indy*dx;
zc = linspace(0,1,nz);
%zc = h5read(path,'/scales/z/1.0');

%[Xc,Yc] = meshgrid(xc,yc);
[Xc,Zc] = meshgrid(xc,zc);

[a,center] = min(abs(yc));

% figure('Position', [10 10 1400 500]);
% 
% subplot(1,2,1)
% contourf(Xc,Yc,squeeze(Q_mean(2,:,:)),'Edgecolor','none'); hold on;
% [C,h] = contour(Xc,Yc,squeeze(w_mean(5,:,:)),'k--','Linewidth',1.2);
% clabel(C,h);
% colorbar
% caxis([-max(max(squeeze(Q_mean(2,:,:)))) max(max(squeeze(Q_mean(2,:,:))))])
% colormap(redblue(10))
% 
% subplot(1,2,2)
% contourf(Xc,Yc,squeeze(Q_mean(7,:,:)),'Edgecolor','none'); hold on;
% [C,h] = contour(Xc,Yc,squeeze(w_mean(5,:,:)),'k--','Linewidth',1.2);
% clabel(C,h);
% colorbar
% caxis([-max(max(squeeze(Q_mean(2,:,:)))) max(max(squeeze(Q_mean(2,:,:))))])
% colormap(redblue(10))

%save('drv_composite_eps001.mat','Xc','Zc','center','Q_mean','Q_diab_z_mean','Q_adv_z_mean');

figure('Position', [10 10 700 500]);

%tend = Q_diab_z_mean + Q_adv_z_mean;
%tend = Q_diab_z_mean;
tend = Q_diab_hoskins_mean;
%tend = theta_prime_ensemble;

cmax = max(max(abs(Q_mean(:,center,:))));
%cmax = 6;
ncont = 15;
%ncont = 100;
cint = cmax/ncont;

contourf(Xc,Zc,squeeze(Q_mean(:,center,:)),-cmax:cint:cmax,'Edgecolor','none'); hold on;
%contour(Xc,Zc,squeeze(tend(:,center,:)),[cint:cint:cint*ncontours],'r-','Linewidth',1.2);
%contour(Xc,Zc,squeeze(tend(:,center,:)),[-cint*ncontours:cint:-cint],'b-','Linewidth',1.2);

ncontours = 20;
cint = max(max(tend(:,center,:)))/ncontours;
%contour(Xc,Zc,squeeze(tend(:,center,:)),[cint:cint:cint*ncontours],'r','Linewidth',1.2);
%contour(Xc,Zc,squeeze(tend(:,center,:)),[-cint*ncontours:cint:-cint],'b','Linewidth',1.2);
colorbar
%caxis([-max(max(squeeze(Q_mean(:,center,:)))) max(max(squeeze(Q_mean(:,center,:))))])
caxis([-cmax cmax])
%caxis([-cmax cmax])
%colormap(redblue(15))
colormap(redblue(ncont));
xlabel('x')
ylabel('z')
title('\rm PV Anomaly')

figure 
nint = 15;
cmax = max(max(abs(Q_mean(:,center,:))));
%cmax = 6;
cinta = cmax/nint;

contourf(Xc,Zc,squeeze(Q_full_mean(:,center,:)),-cmax:cinta:cmax,'Edgecolor','none'); hold on;
%contour(Xc,Zc,squeeze(tend(:,center,:)),[cint:cint:cint*ncontours],'r-','Linewidth',1.2);
%contour(Xc,Zc,squeeze(tend(:,center,:)),[-cint*ncontours:cint:-cint],'b-','Linewidth',1.2);
colorbar
%caxis([-max(max(squeeze(Q_mean(:,center,:)))) max(max(squeeze(Q_mean(:,center,:))))])
caxis([-cmax cmax])
%colormap(redblue(15))
colormap(redblue(nint));
title('\rm Full PV')


% figure('Position', [10 10 800 400]);
% 
% subplot(3,3,1)
% contourf(Xc,Zc,squeeze(Qt_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(abs(Qt_mean(:,center,:))))) max(max(squeeze(abs(Qt_mean(:,center,:)))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $Q_{t}$','Interpreter','latex')
% 
% 
% subplot(3,3,2)
% contourf(Xc,Zc,squeeze(Q_adv_x_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(abs(Q_adv_x_mean(:,center,:))))) max(max(squeeze(abs(Q_adv_x_mean(:,center,:)))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $-uQ_{x}$','Interpreter','latex')
% 
% 
% subplot(3,3,3)
% contourf(Xc,Zc,squeeze(Q_adv_y_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(abs(Q_adv_y_mean(:,center,:))))) max(max(squeeze(abs(Q_adv_y_mean(:,center,:)))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $-vQ_{y}$','Interpreter','latex')
% 
% subplot(3,3,4)
% contourf(Xc,Zc,squeeze(Q_adv_z_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(abs(Q_adv_z_mean(:,center,:))))) max(max(squeeze(abs(Q_adv_z_mean(:,center,:)))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $-wQ_{z}$','Interpreter','latex')
% 
% subplot(3,3,5)
% contourf(Xc,Zc,squeeze(Q_diab_x_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(abs(Q_diab_x_mean(:,center,:))))) max(max(squeeze(abs(Q_diab_x_mean(:,center,:)))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $\dot{Q}_x$','Interpreter','latex')
% 
% subplot(3,3,6)
% contourf(Xc,Zc,squeeze(Q_diab_y_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(abs(squeeze(Q_diab_y_mean(:,center,:))))) max(max(squeeze(abs(Q_diab_y_mean(:,center,:)))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $\dot{Q}_y$','Interpreter','latex')
% 
% subplot(3,3,7)
% contourf(Xc,Zc,squeeze(Q_diab_z_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(abs(Q_diab_z_mean(:,center,:))))) max(max(squeeze(abs(Q_diab_z_mean(:,center,:)))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $\dot{Q}_z$','Interpreter','latex')
% 
% subplot(3,3,8)
% contourf(Xc,Zc,squeeze(residual(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(abs(residual(:,center,:))))) max(max(squeeze(abs(residual(:,center,:)))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $Residual$','Interpreter','latex')

% figure('Position', [10 10 800 400]);
% 
% subplot(3,3,1)
% contourf(Xc,Zc,squeeze(Qt_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(Qt_mean(:,center,:)))) max(max(squeeze(Qt_mean(:,center,:))))])
% caxis([-100 100])
% colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $Q_{t}$','Interpreter','latex')
% 
% 
% subplot(3,3,2)
% contourf(Xc,Zc,squeeze(Q_adv_x_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(Q_adv_x_mean(:,center,:)))) max(max(squeeze(Q_adv_x_mean(:,center,:))))])
% caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $-uQ_{x}$','Interpreter','latex')
% 
% 
% subplot(3,3,3)
% contourf(Xc,Zc,squeeze(Q_adv_y_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(Q_adv_y_mean(:,center,:)))) max(max(squeeze(Q_adv_y_mean(:,center,:))))])
% caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $-vQ_{y}$','Interpreter','latex')
% 
% subplot(3,3,4)
% contourf(Xc,Zc,squeeze(Q_adv_z_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([min(min(squeeze(Q_adv_z_mean(:,center,:)))) -min(min(squeeze(Q_adv_z_mean(:,center,:))))])
% caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $-wQ_{z}$','Interpreter','latex')
% 
% subplot(3,3,5)
% contourf(Xc,Zc,squeeze(Q_diab_x_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(Q_diab_x_mean(:,center,:)))) max(max(squeeze(Q_diab_x_mean(:,center,:))))])
% caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $\dot{Q}_x$','Interpreter','latex')
% 
% subplot(3,3,6)
% contourf(Xc,Zc,squeeze(Q_diab_y_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(Q_diab_y_mean(:,center,:)))) max(max(squeeze(Q_diab_y_mean(:,center,:))))])
% caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $\dot{Q}_y$','Interpreter','latex')
% 
% subplot(3,3,7)
% contourf(Xc,Zc,squeeze(Q_diab_z_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(Q_diab_z_mean(:,center,:)))) max(max(squeeze(Q_diab_z_mean(:,center,:))))])
% caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $\dot{Q}_z$','Interpreter','latex')
% 
% subplot(3,3,8)
% contourf(Xc,Zc,squeeze(residual(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([min(min(squeeze(residual(:,center,:)))) -min(min(squeeze(residual(:,center,:))))])
% caxis([-100 100])
% colormap(redblue(10))
% set(gca,'FontSize',12)
% title('\rm $Residual$','Interpreter','latex')



% figure('Position', [10 10 800 400]);
% 
% subplot(3,3,1)
% contourf(Xc,Zc,squeeze(theta_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(Qt_mean(:,center,:)))) max(max(squeeze(Qt_mean(:,center,:))))])
% %caxis([-100 100])
% %colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $\theta$','Interpreter','latex')
% 
% subplot(3,3,2)
% contourf(Xc,Zc,squeeze(theta_z_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(theta_z_mean(:,center,:)))) max(max(squeeze(theta_z_mean(:,center,:))))])
% %caxis([-100 100])
% %colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $\theta_z$','Interpreter','latex')
% 
% subplot(3,3,3)
% contourf(Xc,Zc,squeeze(theta_dot_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(Qt_mean(:,center,:)))) max(max(squeeze(Qt_mean(:,center,:))))])
% %caxis([-100 100])
% %colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $\dot{\theta}$','Interpreter','latex')
% 
% subplot(3,3,4)
% contourf(Xc,Zc,squeeze(theta_dot_z_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(theta_dot_z_mean(:,center,:)))) max(max(squeeze(theta_dot_z_mean(:,center,:))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $\dot{\theta}_z$','Interpreter','latex')
% 
% subplot(3,3,5)
% contourf(Xc,Zc,squeeze(Qx_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(Qx_mean(:,center,:)))) max(max(squeeze(Qx_mean(:,center,:))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $Q_x$','Interpreter','latex')
% 
% subplot(3,3,6)
% contourf(Xc,Zc,squeeze(Qy_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(Qy_mean(:,center,:)))) max(max(squeeze(Qy_mean(:,center,:))))])
% %caxis([-100 100])
% colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $Q_y$','Interpreter','latex')
% 
% subplot(3,3,7)
% contourf(Xc,Zc,squeeze(Qz_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% %caxis([-max(max(squeeze(Qt_mean(:,center,:)))) max(max(squeeze(Qt_mean(:,center,:))))])
% %caxis([-100 100])
% %colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $Q_z$','Interpreter','latex')
% 
% subplot(3,3,8)
% contourf(Xc,Zc,squeeze(w_mean(:,center,:)),'Edgecolor','none'); hold on
% colorbar;
% caxis([-max(max(squeeze(w_mean(:,center,:)))) max(max(squeeze(w_mean(:,center,:))))])
% %caxis([-100 100])
% %colormap(redblue(10))
% set(gca,'Fontsize',12)
% title('\rm $w$','Interpreter','latex')

%path = '/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/primitive_equations_constant_strat_r001_eps04_final';

% if ii ==1
% save([pwd,'/DRV_storms_eps04_w_composite.mat'],'center','Xc','Zc','Q_mean','Q_full_mean','Q_diab_hoskins_mean','Q_diab_z_mean','tend');
% else
% save([pwd,'/DRV_storms_eps001_w_composite.mat'],'center','Xc','Zc','Q_mean','Q_full_mean','Q_diab_hoskins_mean','Q_diab_z_mean','tend');
% end

if ii ==1
save([pwd,'/DRV_storms_eps04_w_composite_final.mat'],'center','Xc','Zc','Q_mean','Q_full_mean','Q_diab_z_mean','Q_adv_x_mean','Q_adv_y_mean','Q_adv_z_mean','Q_diab_hoskins_mean','theta_z_mean','theta_dot_mean','tend','w_mean','t_start');
elseif ii==2
save([pwd,'/DRV_storms_eps01_w_composite_final.mat'],'center','Xc','Zc','Q_mean','Q_full_mean','Q_diab_z_mean','Q_adv_x_mean','Q_adv_y_mean','Q_adv_z_mean','Q_diab_hoskins_mean','theta_z_mean','theta_dot_mean','tend','w_mean','t_start');  
else
save([pwd,'/DRV_storms_eps001_w_composite_final.mat'],'center','Xc','Zc','Q_mean','Q_full_mean','Q_diab_z_mean','Q_adv_x_mean','Q_adv_y_mean','Q_adv_z_mean','Q_diab_hoskins_mean','theta_z_mean','theta_dot_mean','tend','w_mean','t_start');
end

% if ii ==1
% save([pwd,'/DRV_storms_eps04_q_composite.mat'],'center','Xc','Zc','Q_mean','Q_full_mean','Q_diab_z_mean','Q_adv_z_mean','Q_diab_hoskins_mean','theta_z_mean','tend','w_mean');
% elseif ii==2
% save([pwd,'/DRV_storms_eps01_q_composite.mat'],'center','Xc','Zc','Q_mean','Q_full_mean','Q_diab_z_mean','Q_adv_z_mean','Q_diab_hoskins_mean','theta_z_mean','tend','w_mean');  
% else
% save([pwd,'/DRV_storms_eps001_q_composite.mat'],'center','Xc','Zc','Q_mean','Q_full_mean','Q_diab_z_mean','Q_adv_z_mean','Q_diab_hoskins_mean','theta_z_mean','tend','w_mean');
% end

end
