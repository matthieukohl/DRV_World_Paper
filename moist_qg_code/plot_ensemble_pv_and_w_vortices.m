% make an ensemble plot of pv and w for the vortices in the r=0.01
% simulation

close all; clear;

%list = {'qg_r01_output_tendencies_with_qt_different_method'};
list = {'qg_r001_output_tendencies_with_qt_different_method_linear_damping_alpha_1p7'};


R = [0.01,0.03,0.05,0.1,0.2,0.4,0.6,0.8,1.0];

path_qg = '/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/full_relax_correct/equations_correct_nd/qg_final/';

for kk = 1:1 %length(list)
kk

path = [path_qg,list{kk},'/snapshots/snapshots_s1.h5'];

%path = [path,'/snapshots/snapshots_s1.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta_bc');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nt = h5_dims(1);

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
dx = x(2)-x(1);
t = h5read(path,'/scales/sim_time');

% define the time-tendencies of PV
zeta_bc = h5read(path,'/tasks/zeta_bc');
zeta_bt = h5read(path,'/tasks/zeta_bt');
tau = h5read(path,'/tasks/tau');

q1 = zeta_bt+zeta_bc-tau;
q2 = zeta_bt-zeta_bc+tau;

clear('zeta_bc','zeta_bt','tau')

q1t = h5read(path,'/tasks/q1t');
q2t = h5read(path,'/tasks/q2t');

% T = permute(repmat(t,1,nx,nx),[2,3,1]);
% 
% % % backwards differencing
% q1t = (q1(:,:,2:end)-q1(:,:,1:end-1))./(T(:,:,2:end)-T(:,:,1:end-1));
% q2t = (q2(:,:,2:end)-q2(:,:,1:end-1))./(T(:,:,2:end)-T(:,:,1:end-1));
% 
% q1t = cat(3,zeros(nx,nx),q1t);
% q2t = cat(3,zeros(nx,nx),q2t);

% central differencing

% q1t = (q1(:,:,3:end)-q1(:,:,1:end-2))./(T(:,:,3:end)-T(:,:,1:end-2));
% q2t = (q2(:,:,3:end)-q2(:,:,1:end-2))./(T(:,:,3:end)-T(:,:,1:end-2));
% 
% q1t = cat(3,zeros(nx,nx),q1t);
% q2t = cat(3,zeros(nx,nx),q2t);
% 
% q1t = cat(3,q1t,zeros(nx,nx));
% q2t = cat(3,q2t,zeros(nx,nx));
% 
% clear('q1','q2')

% define size of ensemble mean

indx = 13; indy = indx;
%indx = 26; indy = indx;

tau_mean = zeros(2*indx+1,2*indx+1);
q1_mean = zeros(2*indx+1,2*indx+1);
q2_mean = zeros(2*indx+1,2*indx+1);
q1t_mean = zeros(2*indx+1,2*indx+1);
q2t_mean = zeros(2*indx+1,2*indx+1);
q1_rad_mean = zeros(2*indx+1,2*indx+1);
q2_rad_mean = zeros(2*indx+1,2*indx+1);
q1x_mean = zeros(2*indx+1,2*indx+1);
q2x_mean = zeros(2*indx+1,2*indx+1);
q1_adv_mean = zeros(2*indx+1,2*indx+1);
q2_adv_mean = zeros(2*indx+1,2*indx+1);


v1_mean = zeros(2*indx+1,2*indx+1);
v2_mean = zeros(2*indx+1,2*indx+1);

w_mean = zeros(2*indx+1,2*indx+1);
diab1_mean = zeros(2*indx+1,2*indx+1);
diab2_mean = zeros(2*indx+1,2*indx+1);

rad1_mean = zeros(2*indx+1,2*indx+1);
rad2_mean = zeros(2*indx+1,2*indx+1);

drag_mean = zeros(2*indx+1,2*indx+1);

% pv gradients

beta = 0.78;
q1y = 1+beta;
q2y = -1+beta;

% mean winds
U1 = 1; U2 = -1;

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
for tt = 10:length(t)
    
tt
    
start = [1 1 tt];
count = [nx nx 1];
 
zeta_bc = h5read(path,'/tasks/zeta_bc',start,count);
zeta_bt = h5read(path,'/tasks/zeta_bt',start,count);
tau = h5read(path,'/tasks/tau',start,count);

q1 = zeta_bt+zeta_bc-tau;
q2 = zeta_bt-zeta_bc+tau;

q1x = h5read(path,'/tasks/q1x',start,count);
q2x = h5read(path,'/tasks/q2x',start,count);

q1_adv = h5read(path,'/tasks/q1_adv',start,count);
q2_adv = h5read(path,'/tasks/q2_adv',start,count);

% streamfunction definition psi_i in J(psi_i,q_i) was off by a factor 2 (now corrected). 

q1_adv = -q1_adv;
q2_adv = -q2_adv;

v1 = h5read(path,'/tasks/v1',start,count);
v2 = h5read(path,'/tasks/v2',start,count);
w = h5read(path,'/tasks/w',start,count);

rr = ones(size(w));
rr(w>0) = r;
% diab1 = -(1-rr).*w;
% diab2 = (1-rr).*w;

diab = h5read(path,'/tasks/q_diab',start,count);
diab1 = -diab;
diab2 = diab;

rad = h5read(path,'/tasks/q_rad',start,count);
q1_rad = rad;
q2_rad = -rad;

drag = h5read(path,'/tasks/q_drag',start,count);
drag = 2*drag;

if nvar==1
z = w;
else
z = q2;
end

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

% Plot w with markers for its peaks
% figure
% contourf(X,Y,z);
% caxis([-max(max(z)) max(max(z))])
% colorbar
% colormap(redblue)
% hold on
% plot(X(peak_inds),Y(peak_inds),'+','MarkerSize',24);
%contourf(X(peak_inds),Y(peak_inds),z(peak_inds),'r*','MarkerSize',24)

% Save data to sruct
peakdata = struct;
peakdata.peakZ = z(peak_inds);
peakdata.peakX = X(peak_inds);
peakdata.peakY = Y(peak_inds);
peakdata.peakXWidth = xwidths(peak_inds);
peakdata.peakYWidth = ywidths(peak_inds);

% diff_peak_x = diff(peakdata.peakX); diff_peak_y = diff(peakdata.peakY);
% diff_peak = 0.5*(diff_peak_x+diff_peak_y);
% cut_off = 10;
% 
% peakdata.peakZ = peakdata.peakZ(diff_peak>cut_off);
% peakdata.peakX = peakdata.peakX(diff_peak>cut_off);
% peakdata.peakY = peakdata.peakY(diff_peak>cut_off);
% peakdata.peakXWidth = peakdata.peakXWidth(diff_peak>cut_off);
% peakdata.peakYWidth = peakdata.peakYWidth(diff_peak>cut_off);

% sort into highest peaks

[peak_height,I] = sort(peakdata.peakZ,'descend');
[peak_x] = peakdata.peakX(I);
[peak_y] = peakdata.peakY(I);
[peak_xwidth] = peakdata.peakXWidth(I);
[peak_ywidth] = peakdata.peakYWidth(I);


% Average over the nth biggest vortices

peak_x = peak_x(1:n); peak_y = peak_y(1:n);
%peak_x = peak_x(1); peak_y = peak_y(1);

% peak_x = peak_x(1:6); peak_y = peak_y(1:6);
% 
% diff_peak_x = abs(diff(peak_x)); diff_peak_y = abs(diff(peak_y));
% diff_peak = 0.5*(diff_peak_x+diff_peak_y);
% cut_off = 1;
% 
% peak_x = peak_x(diff_peak>cut_off);
% peak_y = peak_y(diff_peak>cut_off);
% peak_xwidth = peak_xwidth(diff_peak>cut_off);
% peak_ywidth = peak_ywidth(diff_peak>cut_off);
% 
% if length(peak_x)<=1
%     ratio_vs_time(tt) = NaN;
%     continue
% end

%plot w and mark the peaks identified
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
 [a,indcy] = min(abs(x-peak_y(i)));
 
 % exclude points that lie on the boundaries
 if indcx+indx>=512 || indcx-indx<=0 || indcy+indy>=512 || indcy-indy<=0
     continue
 end
 tau_mean = tau_mean + tau(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 q1_mean = q1_mean + q1(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 q2_mean = q2_mean + q2(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 w_mean = w_mean + w(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 v1_mean = v1_mean + v1(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 v2_mean = v2_mean + v2(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 q1x_mean = q1x_mean + q1x(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 q2x_mean = q2x_mean + q2x(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 q1_adv_mean = q1_adv_mean + q1_adv(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 q2_adv_mean = q2_adv_mean + q2_adv(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 q1t_mean = q1t_mean + q1t(indcy-indy:indcy+indy,indcx-indx:indcx+indx,tt);
 q2t_mean = q2t_mean + q2t(indcy-indy:indcy+indy,indcx-indx:indcx+indx,tt);
 diab1_mean = diab1_mean + diab1(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 diab2_mean = diab2_mean + diab2(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 rad1_mean = rad1_mean + q1_rad(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 rad2_mean = rad2_mean + q2_rad(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 
 drag_mean = drag_mean + drag(indcy-indy:indcy+indy,indcx-indx:indcx+indx);
 
  
counter = counter + 1;
end



end

% take average 

tau_mean = tau_mean/counter;

q1_mean = q1_mean/counter;
q2_mean = q2_mean/counter;

q1x_mean = q1x_mean/counter;
q2x_mean = q2x_mean/counter;

q1_adv_mean = q1_adv_mean/counter;
q2_adv_mean = q2_adv_mean/counter;

q1t_mean = q1t_mean/counter;
q2t_mean = q2t_mean/counter;

diab1_mean = diab1_mean/counter;
diab2_mean = diab2_mean/counter;

rad1_mean = rad1_mean/counter;
rad2_mean = rad2_mean/counter;

drag_mean = drag_mean/counter;

v1_mean = v1_mean/counter;
v2_mean = v2_mean/counter;
w_mean = w_mean/counter;

% define tendencies

q1x_mean = -U1*q1x_mean;
q2x_mean = -U2*q2x_mean;

q1y_mean = -q1y*v1_mean;
q2y_mean = -q2y*v2_mean;

drag1_mean = 0*drag_mean;
drag2_mean = -drag_mean;

% define residual (averaged)
residual1 = q1t_mean-q1x_mean-q1y_mean-q1_adv_mean-diab1_mean-rad1_mean-drag1_mean;
residual = q2t_mean-q2x_mean-q2y_mean-q2_adv_mean-diab2_mean-rad2_mean-drag2_mean;

% make a plot of the pv ensemble
xc = -indx*dx:dx:indx*dx;
yc = -indy*dx:dx:indy*dx;

[Xc,Yc] = meshgrid(xc,yc);

figure('Position', [10 10 1000 300]);

%cmax = max(abs(q2_mean(:)));
cmax = 20;
ncont = 20;
cint = cmax/ncont;


subplot(1,2,1)
contourf(Xc,Yc,q2_mean,[-cint*ncont:cint:cint*ncont],'Edgecolor','none'); hold on;
%contour(Xc,Yc,w_mean,[max(max(w_mean))-80:10:max(max(w_mean))],'k--','Linewidth',1.2);
[C,h] = contour(Xc,Yc,w_mean,[0,20,60,100],'k--','Linewidth',1.2);
%[C,h] = contour(Xc,Yc,tau_mean,[0,0.2,0.4,0.5],'k--','Linewidth',1.2);
clabel(C,h);
xlabel('x')
ylabel('y')
colorbar
%caxis([-max(max(q2_mean)) max(max(q2_mean))])
caxis([-cmax cmax])
colormap(redblue(ncont))
title('\rm Lower Layer')
set(gca,'FontSize',12)

subplot(1,2,2)
contourf(Xc,Yc,q1_mean,[-cint*ncont:cint:cint*ncont],'Edgecolor','none'); hold on
colorbar;
xlabel('x')
ylabel('y')
%caxis([-max(max(q2_mean)) max(max(q2_mean))])
caxis([-cmax cmax])
colormap(redblue(ncont))
%contour(Xc,Yc,w_mean,[max(max(w_mean))-80:10:max(max(w_mean))],'k--','Linewidth',1.2)
[C,h] = contour(Xc,Yc,w_mean,[0,20,60,100],'k--','Linewidth',1.2);
%[C,h] = contour(Xc,Yc,tau_mean,[0,0.2,0.4,0.5],'k--','Linewidth',1.2);
clabel(C,h);
title('\rm Upper Layer')
set(gca,'FontSize',12)

annotation('textbox',[0.09, 0.93, 0.01,0.03],'String',"(a)",'EdgeColor','none');
annotation('textbox',[0.53, 0.93, 0.01,0.03],'String',"(b)",'EdgeColor','none');

path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
   'primitive_equations_constant_strat_r001_eps04_geostrophic/figures_drv_world/qg_drv_pv_and_w_linear_damping_alpha_1p7'];
 
% path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
% 'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
%    'primitive_equations_constant_strat_r001_eps06_geostrophic_with_wdry/figures_committee_meeting_aug26th/qg_drv_pv_and_w'];
%  

saveas(gca,path_save,'epsc');

ind_cross = round(length(peak_x)/2);

figure('Position', [10 10 1000 400]);

cmax = max(abs(q2t_mean(:)));
nint = 9;
cint = cmax/nint;

subplot(2,3,1)
contourf(Xc,Yc,q2t_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
colorbar;
caxis([-cmax cmax])
%xlabel('x')
ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'Fontsize',12)
%title('\rm $q_{2t}$','Interpreter','latex')
title('\rm PV Tend.','Interpreter','latex')

cmax = max(abs(q2x_mean(:)));
cint = cmax/cint;

subplot(2,3,2)
contourf(Xc,Yc,q2x_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
%xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $q_{2x}$','Interpreter','latex')
title('\rm Mean Zon. Adv.','Interpreter','latex')
set(cbr,'YTick',-80:40:80)

cmax = max(abs(q2y_mean(:)));
cint = cmax/nint;

subplot(2,3,3)
contourf(Xc,Yc,q2y_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
%xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-v_2 \bar{q}_{2y}$','Interpreter','latex')
title('\rm Mean Merid. Adv.','Interpreter','latex')
set(cbr,'YTick',-0.5:0.25:0.5)

cmax = max(abs(q2_adv_mean(:)));
cint = cmax/nint;

subplot(2,3,4)
contourf(Xc,Yc,q2_adv_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
colorbar;
caxis([-cmax cmax])
xlabel('x')
ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-J(\psi_2,q_2)$','Interpreter','latex')
title('\rm Nonlin. Adv.','Interpreter','latex')

cmax = max(abs(diab2_mean(:)));
cint = cmax/nint;

subplot(2,3,5)
contourf(Xc,Yc,diab2_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
colorbar;
caxis([-cmax cmax])
xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $(1-r(w))w$','Interpreter','latex')
title('\rm Diabatic','Interpreter','latex')

cmax = max(abs(drag2_mean(:)));
%cmax = 5;
cint = cmax/nint;

subplot(2,3,6)
contourf(Xc,Yc,drag2_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-R \nabla^2 \psi_2$','Interpreter','latex')
title('\rm Drag','Interpreter','latex')
set(cbr,'YTick',-3:1.5:3)

% a,b,c,d

annotation('textbox',[0.08, 0.93, 0.01,0.03],'String',"(a)",'EdgeColor','none');
annotation('textbox',[0.36, 0.93, 0.01,0.03],'String',"(b)",'EdgeColor','none');
annotation('textbox',[0.64, 0.93, 0.01,0.03],'String',"(c)",'EdgeColor','none');
annotation('textbox',[0.08, 0.45, 0.01,0.03],'String',"(d)",'EdgeColor','none');
annotation('textbox',[0.36, 0.45, 0.01,0.03],'String',"(e)",'EdgeColor','none');
annotation('textbox',[0.64, 0.45, 0.01,0.03],'String',"(f)",'EdgeColor','none');


path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
   'primitive_equations_constant_strat_r001_eps04_geostrophic/figures_drv_world/qg_drv_pv_tendencies_alpha_1p7'];
 
saveas(gca,path_save,'epsc');


figure('Position', [10 10 1000 600]);

cmax = max(abs(q2t_mean(:)));
nint = 9;
cint = cmax/nint;

subplot(3,3,1)
contourf(Xc,Yc,q2t_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
colorbar;
caxis([-cmax cmax])
%xlabel('x')
ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'Fontsize',12)
%title('\rm $q_{2t}$','Interpreter','latex')
title('\rm PV Tend.','Interpreter','latex')

cmax = max(abs(q2x_mean(:)));
cint = cmax/cint;

subplot(3,3,2)
contourf(Xc,Yc,q2x_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
%xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $q_{2x}$','Interpreter','latex')
title('\rm Mean Zon. Adv.','Interpreter','latex')
set(cbr,'YTick',-80:40:80)

cmax = max(abs(q2y_mean(:)));
cint = cmax/nint;

subplot(3,3,3)
contourf(Xc,Yc,q2y_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
%xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-v_2 \bar{q}_{2y}$','Interpreter','latex')
title('\rm Mean Merid. Adv.','Interpreter','latex')
set(cbr,'YTick',-0.5:0.25:0.5)

cmax = max(abs(q2_adv_mean(:)));
cint = cmax/nint;

subplot(3,3,4)
contourf(Xc,Yc,q2_adv_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
colorbar;
caxis([-cmax cmax])
%xlabel('x')
ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-J(\psi_2,q_2)$','Interpreter','latex')
title('\rm Nonlin. Adv.','Interpreter','latex')

cmax = max(abs(diab2_mean(:)));
cint = cmax/nint;

subplot(3,3,5)
contourf(Xc,Yc,diab2_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
colorbar;
caxis([-cmax cmax])
%xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $(1-r(w))w$','Interpreter','latex')
title('\rm Diabatic','Interpreter','latex')

cmax = max(abs(drag2_mean(:)));
%cmax = 5;
cint = cmax/nint;

subplot(3,3,6)
contourf(Xc,Yc,drag2_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
%xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-R \nabla^2 \psi_2$','Interpreter','latex')
title('\rm Drag','Interpreter','latex')
set(cbr,'YTick',-3:1.5:3)

cmax = max(abs(q2_mean(:)));
%cmax = 5;
cint = cmax/nint;

subplot(3,3,7)
contourf(Xc,Yc,q2_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
xlabel('x')
ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-R \nabla^2 \psi_2$','Interpreter','latex')
title('\rm PV','Interpreter','latex')
%set(cbr,'YTick',-3:1.5:3)

cmax = max(abs(w_mean(:)));
%cmax = 5;
cint = cmax/nint;

subplot(3,3,8)
contourf(Xc,Yc,w_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-R \nabla^2 \psi_2$','Interpreter','latex')
title('\rm Vertical Velocity','Interpreter','latex')
%set(cbr,'YTick',-3:1.5:3)

cmax = max(abs(v2_mean(:)));
%cmax = 5;
cint = cmax/nint;

subplot(3,3,9)
contourf(Xc,Yc,v2_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-R \nabla^2 \psi_2$','Interpreter','latex')
title('\rm Meridional Velocity','Interpreter','latex')
%set(cbr,'YTick',-3:1.5:3)

% a,b,c,d

annotation('textbox',[0.08, 0.92, 0.01,0.03],'String',"(a)",'EdgeColor','none');
annotation('textbox',[0.36, 0.92, 0.01,0.03],'String',"(b)",'EdgeColor','none');
annotation('textbox',[0.64, 0.92, 0.01,0.03],'String',"(c)",'EdgeColor','none');
annotation('textbox',[0.08, 0.62, 0.01,0.03],'String',"(d)",'EdgeColor','none');
annotation('textbox',[0.36, 0.62, 0.01,0.03],'String',"(e)",'EdgeColor','none');
annotation('textbox',[0.64, 0.62, 0.01,0.03],'String',"(f)",'EdgeColor','none');
annotation('textbox',[0.08, 0.32, 0.01,0.03],'String',"(g)",'EdgeColor','none');
annotation('textbox',[0.36, 0.32, 0.01,0.03],'String',"(h)",'EdgeColor','none');
annotation('textbox',[0.64, 0.32, 0.01,0.03],'String',"(i)",'EdgeColor','none');


path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
 'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
    'primitive_equations_constant_strat_r001_eps04_geostrophic/figures_drv_world/qg_drv_pv_tendencies_with_pv_w_v_linear_damping_alpha_1p7'];
 
saveas(gca,path_save,'epsc');


figure('Position', [10 10 1000 400]);

cmax = max(abs(q1t_mean(:)));
nint = 9;
cint = cmax/nint;

subplot(2,3,1)
contourf(Xc,Yc,q1t_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
colorbar;
caxis([-cmax cmax])
%xlabel('x')
ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'Fontsize',12)
%title('\rm $q_{2t}$','Interpreter','latex')
title('\rm PV Tend.','Interpreter','latex')

cmax = max(abs(q1x_mean(:)));
cint = cmax/cint;

subplot(2,3,2)
contourf(Xc,Yc,q1x_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
%xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $q_{2x}$','Interpreter','latex')
title('\rm Mean Zon. Adv.','Interpreter','latex')
set(cbr,'YTick',-80:40:80)

cmax = max(abs(q1y_mean(:)));
cint = cmax/nint;

subplot(2,3,3)
contourf(Xc,Yc,q1y_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
%xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-v_2 \bar{q}_{2y}$','Interpreter','latex')
title('\rm Mean Merid. Adv.','Interpreter','latex')
%set(cbr,'YTick',-5:0:0.5)

cmax = max(abs(q1_adv_mean(:)));
cint = cmax/nint;

subplot(2,3,4)
contourf(Xc,Yc,q1_adv_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
colorbar;
caxis([-cmax cmax])
xlabel('x')
ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-J(\psi_2,q_2)$','Interpreter','latex')
title('\rm Nonlin. Adv.','Interpreter','latex')

cmax = max(abs(diab1_mean(:)));
cint = cmax/nint;

subplot(2,3,5)
contourf(Xc,Yc,diab1_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
colorbar;
caxis([-cmax cmax])
xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $(1-r(w))w$','Interpreter','latex')
title('\rm Diabatic','Interpreter','latex')

cmax = max(abs(drag1_mean(:)));
%cmax = 5;
cint = cmax/nint;

subplot(2,3,6)
contourf(Xc,Yc,drag1_mean,[-cmax:cint:cmax],'Edgecolor','none'); hold on
cbr = colorbar;
caxis([-cmax cmax])
xlabel('x')
%ylabel('y')
%caxis([-100 100])
colormap(redblue(nint))
set(gca,'FontSize',12)
%title('\rm $-R \nabla^2 \psi_2$','Interpreter','latex')
title('\rm Drag','Interpreter','latex')
%set(cbr,'YTick',-3:1.5:3)

% a,b,c,d

annotation('textbox',[0.08, 0.93, 0.01,0.03],'String',"(a)",'EdgeColor','none');
annotation('textbox',[0.36, 0.93, 0.01,0.03],'String',"(b)",'EdgeColor','none');
annotation('textbox',[0.64, 0.93, 0.01,0.03],'String',"(c)",'EdgeColor','none');
annotation('textbox',[0.08, 0.45, 0.01,0.03],'String',"(d)",'EdgeColor','none');
annotation('textbox',[0.36, 0.45, 0.01,0.03],'String',"(e)",'EdgeColor','none');
annotation('textbox',[0.64, 0.45, 0.01,0.03],'String',"(f)",'EdgeColor','none');



[a,center] = min(abs(yc));
dc = 2;
residual = q2t_mean-q2x_mean-q2y_mean-q2_adv_mean-diab2_mean-rad2_mean-drag_mean;

figure(5)
plot(xc,mean(q2t_mean(center-dc:center+dc,:)),'b','LineWidth',1.2); hold on;
plot(xc,mean(q2x_mean(center-dc:center+dc,:)),'r','LineWidth',1.2); hold on;
plot(xc,mean(q2y_mean(center-dc:center+dc,:)),'g','LineWidth',1.2); hold on;
plot(xc,mean(q2_adv_mean(center-dc:center+dc,:)),'r--','LineWidth',1.2); hold on;
plot(xc,mean(diab2_mean(center-dc:center+dc,:)),'k','LineWidth',1.2); hold on;
%plot(xc,mean(rad2_mean(center-dc:center+dc,:)),'k--','LineWidth',1.2); hold on;
plot(xc,mean(drag2_mean(center-dc:center+dc,:)),'y','LineWidth',1.2); hold on;
%plot(xc,mean(residual(center-dc:center+dc,:)),'m','LineWidth',1.2); hold on;
%legend('$q_{2t}$','$q_{2x}$','$-v\bar{q}_{2y}$','$-J(\psi_2,q_2)$','$(1-r(w))w$','$-R\nabla^2\psi_2$','Residual','Interpreter','latex'); legend boxoff
%legend('PV Tend.','Mean Zon. Adv.','Mean Merid. Adv.','Nonlin. Adv.','Diabatic','Drag','Residua','Interpreter','latex'); legend boxoff
legend('PV Tend.','Mean Zon. Adv.','Mean Merid. Adv.','Nonlin. Adv.','Diabatic','Drag','Interpreter','latex'); legend boxoff
set(gca,'FontSize',12)
xlabel('x')
ylabel('Magnitude')
title('\rm Cross-section PV Tendencies')

set(gca,'box','off')

path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
   'primitive_equations_constant_strat_r001_eps04_geostrophic/figures_drv_world/pv_tendencies_slice_linear_damping_alpha_1p7'];
 
saveas(gca,path_save,'epsc');



% figure(3)
% subplot(2,1,1)
% plot(q1_mean(ind_cross,:),'b','Linewidth',1.4); hold on;
% plot(w_mean(ind_cross,:)/10,'k');
% set(gca,'FontSize',12)
% 
% subplot(2,1,2)
% plot(q2_mean(ind_cross,:),'r','Linewidth',1.4); hold on;
% plot(w_mean(ind_cross,:)/10,'k');
% set(gca,'FontSize',12)

figure('Position', [10 10 800 400]);

subplot(3,3,1)
contourf(Xc,Yc,q2t_mean,'Edgecolor','none'); hold on
colorbar;
caxis([-max(max(q2t_mean)) max(max(q2t_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'Fontsize',12)
title('\rm $q_{2t}$','Interpreter','latex')


subplot(3,3,2)
contourf(Xc,Yc,q2x_mean,'Edgecolor','none'); hold on
colorbar;
caxis([-max(max(q2x_mean)) max(max(q2x_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $q_{2x}$','Interpreter','latex')


subplot(3,3,3)
contourf(Xc,Yc,q2y_mean,'Edgecolor','none'); hold on
colorbar;
caxis([-max(max(q2y_mean)) max(max(q2y_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $q_{2y}$','Interpreter','latex')

subplot(3,3,4)
contourf(Xc,Yc,q2_adv_mean,'Edgecolor','none'); hold on
colorbar;
caxis([min(min(q2_adv_mean)) -min(min(q2_adv_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $J(\psi_2,q_2)$','Interpreter','latex')


subplot(3,3,5)
contourf(Xc,Yc,diab2_mean,'Edgecolor','none'); hold on
colorbar;
caxis([-max(max(diab2_mean)) max(max(diab2_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $\dot{q}_{2,diab}$','Interpreter','latex')

subplot(3,3,6)
contourf(Xc,Yc,rad2_mean,'Edgecolor','none'); hold on
colorbar;
caxis([min(min(rad2_mean)) -min(min(rad2_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $Radiation$','Interpreter','latex')

subplot(3,3,7)
contourf(Xc,Yc,drag2_mean,'Edgecolor','none'); hold on
colorbar;
caxis([min(min(drag2_mean)) -min(min(drag2_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $Drag$','Interpreter','latex')

subplot(3,3,8)
contourf(Xc,Yc,residual,'Edgecolor','none'); hold on
colorbar;
caxis([min(min(residual)) -min(min(residual))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $Residual$','Interpreter','latex')

% path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
% 'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
%    'primitive_equations_constant_strat_r001_eps06_geostrophic_with_wdry/figures_committee_meeting_aug26th/qg_drv_pv_tendencies'];
%  
% saveas(gca,path_save,'epsc');


figure('Position', [10 10 800 400]);

subplot(3,3,1)
contourf(Xc,Yc,q2t_mean,'Edgecolor','none'); hold on
colorbar;
%caxis([-max(max(q2t_mean)) max(max(q2t_mean))])
caxis([-100 100])
colormap(redblue(10))
set(gca,'Fontsize',12)
title('\rm $q_{2t}$','Interpreter','latex')


subplot(3,3,2)
contourf(Xc,Yc,q2x_mean,'Edgecolor','none'); hold on
colorbar;
%caxis([-max(max(q2x_mean)) max(max(q2x_mean))])
caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $q_{2x}$','Interpreter','latex')


subplot(3,3,3)
contourf(Xc,Yc,q2y_mean,'Edgecolor','none'); hold on
colorbar;
%caxis([-max(max(q2y_mean)) max(max(q2y_mean))])
caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $q_{2y}$','Interpreter','latex')

subplot(3,3,4)
contourf(Xc,Yc,q2_adv_mean,'Edgecolor','none'); hold on
colorbar;
%caxis([min(min(q2_adv_mean)) -min(min(q2_adv_mean))])
caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $J(\psi_2,q_2)$','Interpreter','latex')


subplot(3,3,5)
contourf(Xc,Yc,diab2_mean,'Edgecolor','none'); hold on
colorbar;
%caxis([-max(max(diab2_mean)) max(max(diab2_mean))])
caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $\dot{q}_{2,diab}$','Interpreter','latex')

subplot(3,3,6)
contourf(Xc,Yc,rad2_mean,'Edgecolor','none'); hold on
colorbar;
%caxis([min(min(rad2_mean)) -min(min(rad2_mean))])
caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $Radiation$','Interpreter','latex')

subplot(3,3,7)
contourf(Xc,Yc,drag2_mean,'Edgecolor','none'); hold on
colorbar;
%caxis([min(min(drag2_mean)) -min(min(drag2_mean))])
caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $Drag$','Interpreter','latex')

subplot(3,3,8)
contourf(Xc,Yc,residual,'Edgecolor','none'); hold on
colorbar;
%caxis([min(min(residual)) -min(min(residual))])
caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $Residual$','Interpreter','latex')

[a,center] = min(abs(yc));
%residual = q2t_mean-q2x_mean-q2y_mean-q2_adv_mean-diab2_mean-rad2_mean-drag_mean;

figure(4)
plot(xc,q2t_mean(center,:),'b','LineWidth',1.2); hold on;
plot(xc,q2x_mean(center,:),'r','LineWidth',1.2); hold on;
plot(xc,q2y_mean(center,:),'g','LineWidth',1.2); hold on;
plot(xc,q2_adv_mean(center,:),'r--','LineWidth',1.2); hold on;
plot(xc,diab2_mean(center,:),'k','LineWidth',1.2); hold on;
plot(xc,rad2_mean(center,:),'k--','LineWidth',1.2); hold on;
plot(xc,drag2_mean(center,:),'y','LineWidth',1.2); hold on;
plot(xc,residual(center,:),'m','LineWidth',1.2); hold on;
legend('qt','q2x','q2y','q2_{adv}','LH','rad','drag','residual'); legend boxoff


set(gca,'FontSize',12)

[a,center] = min(abs(yc));
dc = 3;
residual = q2t_mean-q2x_mean-q2y_mean-q2_adv_mean-diab2_mean-rad2_mean-drag_mean;

figure(5)
plot(xc,mean(q2t_mean(center-dc:center+dc,:)),'b','LineWidth',1.2); hold on;
plot(xc,mean(q2x_mean(center-dc:center+dc,:)),'r','LineWidth',1.2); hold on;
plot(xc,mean(q2y_mean(center-dc:center+dc,:)),'g','LineWidth',1.2); hold on;
plot(xc,mean(q2_adv_mean(center-dc:center+dc,:)),'r--','LineWidth',1.2); hold on;
plot(xc,mean(diab2_mean(center-dc:center+dc,:)),'k','LineWidth',1.2); hold on;
plot(xc,mean(rad2_mean(center-dc:center+dc,:)),'k--','LineWidth',1.2); hold on;
plot(xc,mean(drag2_mean(center-dc:center+dc,:)),'y','LineWidth',1.2); hold on;
plot(xc,mean(residual(center-dc:center+dc,:)),'m','LineWidth',1.2); hold on;
legend('qt','q2x','q2y','q2_{adv}','LH','rad','drag','residual'); legend boxoff
set(gca,'FontSize',12)


figure('Position', [10 10 800 400]);

subplot(3,3,1)
contourf(Xc,Yc,q1t_mean,'Edgecolor','none'); hold on
colorbar;
caxis([-max(max(q1t_mean)) max(max(q1t_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'Fontsize',12)
title('\rm $q_{1t}$','Interpreter','latex')


subplot(3,3,2)
contourf(Xc,Yc,q1x_mean,'Edgecolor','none'); hold on
colorbar;
caxis([-max(max(q1x_mean)) max(max(q1x_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $q_{2x}$','Interpreter','latex')


subplot(3,3,3)
contourf(Xc,Yc,q1y_mean,'Edgecolor','none'); hold on
colorbar;
caxis([min(min(q1y_mean)) -min(min(q1y_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $q_{2y}$','Interpreter','latex')

subplot(3,3,4)
contourf(Xc,Yc,q1_adv_mean,'Edgecolor','none'); hold on
colorbar;
caxis([min(min(q1_adv_mean)) -min(min(q1_adv_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $J(\psi_1,q_1)$','Interpreter','latex')


subplot(3,3,5)
contourf(Xc,Yc,diab1_mean,'Edgecolor','none'); hold on
colorbar;
caxis([min(min(diab1_mean)) -min(min(diab1_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $\dot{q}_{1,diab}$','Interpreter','latex')

subplot(3,3,6)
contourf(Xc,Yc,rad1_mean,'Edgecolor','none'); hold on
colorbar;
caxis([-max(max(rad1_mean)) max(max(rad1_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $Radiation$','Interpreter','latex')

subplot(3,3,7)
contourf(Xc,Yc,drag1_mean,'Edgecolor','none'); hold on
colorbar;
caxis([min(min(drag1_mean)) -min(min(drag1_mean))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $Drag$','Interpreter','latex')

subplot(3,3,8)
contourf(Xc,Yc,residual1,'Edgecolor','none'); hold on
colorbar;
caxis([-max(max(residual1)) max(max(residual1))])
%caxis([-100 100])
colormap(redblue(10))
set(gca,'FontSize',12)
title('\rm $Residual1$','Interpreter','latex')




end
