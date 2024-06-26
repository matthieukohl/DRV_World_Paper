% plot pauls metric

%close all; clear;

% loop through all QG runs

list = {'r03','r035','r04','r05','r06','r08','r1'};
r_factor = [0.3,0.35,0.4,0.5,0.6,0.8,1];

% list = {'r03','r035','r036','r037','r038','r039','r04','r05','r06','r08','r1'};
% r_factor = [0.3,0.35,0.36,0.37,0.38,0.39,0.4,0.5,0.6,0.8,1];

energy = zeros(size(list));
max_variable = zeros(size(list));

skew_vs_r = zeros(size(list));
lambda_vs_r = zeros(size(list));
metric_drv_norm_vs_r = zeros(size(list));
metric_drv_growth_vs_r = zeros(size(list));
metric_jet_vs_r = zeros(size(list));

skew_q1_vs_r = zeros(size(list));
skew_q2_vs_r = zeros(size(list));

u1_rms_vs_r = zeros(512,length(r_factor));
u1_zonal_vs_r = zeros(512,length(r_factor));

q1_zonal_vs_r = zeros(512,length(r_factor));
q1y_vs_r = zeros(512,length(r_factor));

%nint = 25;


for ii = 1:length(list)
      
part = '/nobackup1c/users/mkohl/qg_';

counter = 0;

ii

for j = [1,2,3,4]

path = [part,list{ii},'_alpha0p15_linear_tau_damping/snapshots/snapshots_s',num2str(j),'.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/w');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nt = h5_dims(1);

% if i ==1
%     nint = nt-5;
% else
%     nint = nt - 1;
% end

%nint = nt-1;

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
y = x;
t = h5read(path,'/scales/sim_time');

for tt = 1:nt
    
%     if i==1 && tt<50
%         continue
%     end

if t(tt)<130 || t(tt)>270
    continue
end
%     
    
% read in the PV and tendency fields

start = [1 1 tt];
count = [nx nx 1];
stride = [1 1 1];

zeta_bc = h5read(path,'/tasks/zeta_bc',start,count,stride);
zeta_bt = h5read(path,'/tasks/zeta_bt',start,count,stride);
phi = h5read(path,'/tasks/tau',start,count,stride);
tau = h5read(path,'/tasks/tau',start,count,stride);
psi1 = phi+tau;
zeta1 = zeta_bt+zeta_bc; 
zeta2 = zeta_bt-zeta_bc;
q1 = zeta1 - tau;
q2 = zeta2 + tau;
beta = 0.78;

energy(ii) = energy(ii) + mean(mean(mean(tau.^2 + zeta_bt.^2 + zeta_bc.^2)));


max_variable(ii) = max_variable(ii) + mean(min(min(zeta1,[],1),[],2),3);


clear('zeta1','zeta2','tau');



% q1 = q1 - mean(q1,2);
% q2 = q2 - mean(q2,2);

% calculate PV gradients

q1_zonal = squeeze(mean(q1,2))+beta*x';
q1_zonal_vs_r(:,ii) = q1_zonal;
q1y = zeros(size(q1_zonal));
dy = 12*pi/nx;
q1y(2:end-1,:) = (q1_zonal(3:end,:)-q1_zonal(1:end-2,:))/(2*dy);
q1y(1,:) = (q1_zonal(2,:)-q1_zonal(1,:))/dy;
q1y(end,:) = (q1_zonal(end,:)-q1_zonal(end-1,:))/dy;

clear('q1_zonal')

%q1y_vs_r(:,ii) = q1y_vs_r(:,ii) + mean(mean(q1y(:,end-nint:end),2),3);

skew_q1_vs_r(ii) = skew_q1_vs_r(ii) + Skew(q1);
skew_q2_vs_r(ii) = skew_q2_vs_r(ii) + Skew(q2);

% calculate u
u1 = zeros(size(psi1));
u1(2:end-1,:,:) = -(psi1(3:end,:,:)-psi1(1:end-2,:,:))/(2*dy);
u1(1,:,:) = -(psi1(2,:,:)-psi1(1,:,:))/dy;
u1(end,:,:) = -(psi1(end,:,:)-psi1(end-1,:,:))/dy;

%u1_rms_vs_r(:,ii) = mean(max(abs(u1(:,:,end-20:end)),[],2),3);
%u1_rms_vs_r(:,ii) = max(max(abs(u1(:,:,end-20:end)),[],2),[],3);

u1_rms_vs_r(:,ii) = u1_rms_vs_r(:,ii) + rms(abs(u1(:,:)),2);


u1_zonal_vs_r(:,ii) = u1_zonal_vs_r(:,ii) + squeeze(mean(u1(:,:),2));


% u1_max_vs_zonal_vs_r(:,ii) = max(max(u1(:,:,end-20:end),[],2),[],3)./u1_vs_r(:,ii);
%metric_jet_vs_r(ii) = rms(u1_rms_vs_r(:,ii)./u1_zonal_vs_r(:,ii));
metric_jet_vs_r(ii) = metric_jet_vs_r(ii) + rms(u1_zonal_vs_r(:,ii)./u1_rms_vs_r(:,ii));
%metric_jet_vs_r(ii) = max(abs(u1_zonal_vs_r(:,ii))./u1_rms_vs_r(:,ii));

w = h5read(path,'/tasks/w',start,count,stride);

lambda_vs_r(ii) = lambda_vs_r(ii) + Lambda(w);

skew_vs_r(ii) = skew_vs_r(ii) + Skew(w(:,:));
 
r = r_factor(ii);
rr = ones(size(w));
rr(w>0) = r;

diab1 = -(1-rr).*w;
diab2 = (1-rr).*w;

metric_drv = q1.*diab1 + q2.*diab2;
metric_eddy = q1.*q1 + q2.*q2;
metric_tend = diab1.^2 + diab2.^2;

% try only bottom layer
% metric_drv = q2.*diab2;
% metric_eddy = q2.*q2;
% metric_tend = diab2.^2;

metric_drv = squeeze(metric_drv);
metric_eddy = squeeze(metric_eddy);
metric_tend = squeeze(metric_tend);

metric_drv(metric_drv<0) = 0;

% average over the 10 biggest maxima in metric_drv;

z = metric_drv;
n = 6;

% smooth on or off
smooth = 0;
smethod = 'movmean';

% smooth window size
nwindow = 10;

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

if isempty(peak_x) ==1
    continue
end

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
 
%  % exclude points that lie on the boundaries
%  if indcx+indx>=nx || indcx-indx<=0 || indcy+indy>=nx || indcy-indy<=0
%      continue
%  end
 
 metric_norm = metric_drv(indcy,indcx)^2./(max(max(metric_eddy,[],1),[],2).*max(max(metric_tend,[],1),[],2));

 metric_growth = metric_drv(indcy,indcx)^2./(max(max(metric_eddy.^2,[],1),[],2));

 metric_drv_norm_vs_r(ii) = metric_drv_norm_vs_r(ii) + metric_norm;

 metric_drv_growth_vs_r(ii) = metric_drv_growth_vs_r(ii) + metric_growth;
 
 
 counter = counter + 1;
end



end
end


energy(ii) = energy(ii)/counter;
max_variable(ii) = max_variable(ii)/counter;
u1_zonal_vs_r(:,ii) = u1_zonal_vs_r(:,ii)/counter;
metric_drv_norm_vs_r(ii) =  metric_drv_norm_vs_r(ii)/counter;
metric_drv_growth_vs_r(ii) = metric_drv_growth_vs_r(ii)/counter;
end




% figure
% semilogx(r_factor,energy,'b-o');
% xlabel('reduction factor r')
% ylabel('energy')
% 
% figure
% semilogx(r_factor,-max_variable,'b-o');
% xlabel('Reduction factor r');
% ylabel('Metric')

figure('Position', [10 10 1200 300]);

y = linspace(0,12*pi,nx);

subplot(1,3,1)

plot(u1_zonal_vs_r(:,1),y,'b','Linewidth',1.2); hold on; 
%plot(u1_zonal_vs_r(:,4),y,'b'); hold on;
%plot(u1_zonal_vs_r(:,7),y,'k--'); hold on;
plot(u1_zonal_vs_r(:,end),y,'r','Linewidth',1.2);
xlabel('u_1')
ylabel('y')
set(gca,'FontSize',12)
set(gca','box','off');
title('\rm Zonal Jets (Upper Layer)')
ylim([y(1) y(end)])
yticks([0 6*pi 12*pi]);
yticklabels({'0','6pi','12pi'})
legend('r=0.3','r=1','Location','NorthWest'); legend boxoff;

subplot(1,3,2)


%metric_drv_norm_vs_r(end) = 0;

semilogx(r_factor,metric_drv_norm_vs_r,'b-o','Linewidth',1.2); hold on;
%semilogx(r_factor,(1-r_factor).^2,'k');
xlabel('Reduction factor r')
ylabel('$\mathcal{M}_1$','Interpreter','latex')
xlim([0 1])
%title('\rm DRV Metric PV Moist QG Simulation')
%title('\rm DRV Metric - PV/PV Dot Normalization')
title('\rm Moist Storm Metric')
set(gca,'FontSize',12);
set(gca','box','off');

% saveas(gcf,['/nfs/pool002/users/mkohl/run_correct/'...
% 'test_runs/adv_correct/full_relax_correct/equations_correct_nd/'...
% 'primitive_equations_final/primitive_equations_constant_strat_r001_eps04_geostrophic/'...
% 'figures_paper/metric_norm'],'epsc')

% saveas(gcf,'drv_metric_qg_simulation','epsc');

subplot(1,3,3)

semilogx(r_factor,metric_drv_growth_vs_r,'b-o','Linewidth',1.2); hold on;
%semilogx(r_factor,(1-r_factor).^2,'k');
xlabel('Reduction factor r')
ylabel('$\mathcal{M}_2$','Interpreter','latex')
title('\rm Growth Rate Metric');
%title('\rm DRV Metric PV Moist QG Simulation')
%title('\rm DRV Metric - PV/PV Dot Normalization')
set(gca,'FontSize',12);
set(gca','box','off');
%saveas(gcf,'drv_metric_qg_simulation','epsc');

annotation('textbox',[0.07, 0.94, 0.01,0.03],'String',"(a)",'EdgeColor','none');
annotation('textbox',[0.35, 0.94, 0.01,0.03],'String',"(b)",'EdgeColor','none');
annotation('textbox',[0.63, 0.94, 0.01,0.03],'String',"(c)",'EdgeColor','none');


saveas(gcf,['/nfs/pool002/users/mkohl/run_correct/'...
'test_runs/adv_correct/full_relax_correct/equations_correct_nd/'...
'primitive_equations_final/primitive_equations_constant_strat_r001_eps04_geostrophic/'...
'figures_paper/transition_metrics_alpha_0p15'],'epsc')

figure
semilogx(r_factor,rms(u1_zonal_vs_r,1),'k-o');
xlabel('Reduction factor r')
ylabel('Metric jet')
set(gca,'FontSize',12);
title('\rm rms of zonal-time mean u velocity')

figure
semilogx(r_factor,metric_jet_vs_r,'k-o');
xlabel('Reduction factor r')
ylabel('Metric jet')
set(gca,'FontSize',12);
title('\rm rms of zonal-time mean u velocity / max u velocity')

y = linspace(0,12*pi,nx);

figure
plot(u1_zonal_vs_r(:,1),y,'k--'); hold on; 
plot(u1_zonal_vs_r(:,3),y,'b'); hold on;
%plot(u1_zonal_vs_r(:,7),y,'k--'); hold on;
plot(u1_zonal_vs_r(:,end),y,'r');
xlabel('u1')
ylabel('y')
set(gca,'FontSize',12)
ylim([y(1) y(end)])
yticks([0 6*pi 12*pi]);
yticklabels({'0','6pi','12pi'})
legend('r=0.01','r=0.3','r=1'); legend boxoff;

figure
y = linspace(0,12*pi,nx);
%plot(q1_zonal_vs_r(:,1),y,'b'); hold on; 
plot(q1_zonal_vs_r(:,1),y,'b'); hold on; 
plot(q1_zonal_vs_r(:,end),y,'r');
xlabel('q1')
ylabel('y')
ylim([y(1) y(end)])
set(gca,'FontSize',12)
legend('r=0.3','r=1','Location','NorthWest'); legend boxoff;
title('\rm Potential Vorticity Top Layer');

figure
y = linspace(0,12*pi,nx);

%q1y_vs_r = movmean(q1y_vs_r,100,1);

plot(q1y_vs_r(:,1),y,'b'); hold on; 
plot(q1y_vs_r(:,end),y,'r');
xlabel('q1y')
ylabel('y')
ylim([y(1) y(end)])
set(gca,'FontSize',12)

figure
semilogx(r_factor,lambda_vs_r,'b-o'); hold on;
semilogx(r_factor,1./(1+sqrt(r_factor)),'r-o');
xlabel('Redution factor')
ylabel('Asymmetry parameter \lambda');
legend('Simulation','Pablo theory'); legend boxoff;

figure
semilogx(r_factor,-skew_q1_vs_r,'b-o'); hold on;
semilogx(r_factor,skew_q2_vs_r,'r-o'); hold on;
xlabel('Redution factor')
ylabel('Skewness');

figure
y = linspace(0,12*pi,nx);
%plot(q1_zonal_vs_r(:,1),y,'b'); hold on; 
plot(q1_zonal_vs_r(:,1),y,'b'); hold on; 
plot(q1_zonal_vs_r(:,end),y,'r');
xlabel('q1')
ylabel('y')
ylim([y(1) y(end)])
set(gca,'FontSize',12)
legend('r=0.01','r=1','Location','NorthWest'); legend boxoff;
title('\rm Potential Vorticity Top Layer');

figure
y = linspace(0,12*pi,nx);

%q1y_vs_r = movmean(q1y_vs_r,100,1);

plot(q1y_vs_r(:,1),y,'b'); hold on; 
plot(q1y_vs_r(:,end),y,'r');
xlabel('q1y')
ylabel('y')
ylim([y(1) y(end)])
set(gca,'FontSize',12)

figure
plot(u1_zonal_vs_r(:,1),y,'b'); hold on; 
plot(u1_zonal_vs_r(:,6),y,'k'); hold on;
plot(u1_zonal_vs_r(:,7),y,'k--'); hold on;
plot(u1_zonal_vs_r(:,end),y,'r');
xlabel('u1')
ylabel('y')
set(gca,'FontSize',12)
ylim([y(1) y(end)])
yticks([0 6*pi 12*pi]);
yticklabels({'0','6pi','12pi'})
legend('r=0.01','r=0.4','r=0.6','r=1'); legend boxoff;


figure
semilogx(r_factor,rms(u1_zonal_vs_r,1),'k-o');
xlabel('Reduction factor r')
ylabel('Metric jet')
set(gca,'FontSize',12);
title('\rm rms of zonal-time mean u velocity')

figure
semilogx(r_factor,rms(q1_zonal_vs_r,1),'k-o');
xlabel('Reduction factor r')
ylabel('Metric jet')
set(gca,'FontSize',12);

figure
semilogx(r_factor,metric_jet_vs_r,'k-o');
xlabel('Reduction factor r')
ylabel('Metric jet')
set(gca,'FontSize',12);
title('\rm rms of zonal-time mean u velocity / max u velocity')


% make spectra of velocities

u1_r001 = u1_vs_r(:,1);
u1_r1 = u1_vs_r(:,end);

spec_r001 = abs(fft(u1_r001)).^2;
spec_r1 = abs(fft(u1_r1)).^2;

loglog(spec_r001(1:nx/2)/max(spec_r001(1:nx/2))); hold on;
loglog(spec_r1(1:nx/2)/max(spec_r1(1:nx/2)));
