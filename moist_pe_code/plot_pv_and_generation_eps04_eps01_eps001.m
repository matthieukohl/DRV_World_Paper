% make a plot of PV composite and generation 

close all; clear;

ncontours = 40; % pv tend contours
nint = 280; % pv contours

figure('Position', [10 10 1300 300]);

subplot(1,3,1)

%load('DRV_storms_eps001_w_composite.mat')
load('DRV_storms_eps001_q_composite.mat')

eps = 0.01;

%tend = Q_diab_z_mean;
tend = Q_diab_z_mean + Q_adv_z_mean;
%tend = Q_diab_hoskins_mean;

%ncontours = 12;
cint = max(max(tend(:,center,:)))/ncontours;
%cint = 0.18/ncontours;
disp(cint)

%nint = 55;
%nint = 25;
%cmax = max(max(abs(Q_mean(:,center,:))));
cmax = 0.12;
%cmax = 19;
cinta = cmax/nint;

%tend = Q_diab_z_mean + Q_adv_z_mean;

contourf(Xc,Zc,squeeze(Q_mean(:,center,:)),-cmax:cinta:cmax,'Edgecolor','none'); hold on;
contour(Xc,Zc,squeeze(tend(:,center,:)),[cint:cint:cint*ncontours],'r-','Linewidth',1.2);
contour(Xc,Zc,squeeze(tend(:,center,:)),[-cint*ncontours:cint:-cint],'b-','Linewidth',1.2);
colorbar
%caxis([-max(max(squeeze(Q_mean(:,center,:)))) max(max(squeeze(Q_mean(:,center,:))))])
caxis([-cmax cmax])
%colormap(redblue(15))
colormap(redblue(nint));
title('\rm \epsilon=0.01');
set(gca,'FontSize',12)
xlabel('x')
ylabel('z')


subplot(1,3,2)

%load('DRV_storms_eps01_w_composite.mat')
load('DRV_storms_eps01_q_composite.mat')

eps = 0.1;

%tend = Q_diab_z_mean;
tend = Q_diab_z_mean + Q_adv_z_mean;
%tend = Q_diab_hoskins_mean;
%tend = theta_prime_ensemble;

%ncontours = 12;
%ncontours = 12;
cint = max(max(tend(:,center,:)))/ncontours;
%cint = 6/ncontours;
disp(cint)

%nint = 55;
%nint = 25;
%cmax = max(max(abs(Q_mean(:,center,:))));
cmax = 1;
%cmax = 6;
%cmax = 19;
cinta = cmax/nint;

contourf(Xc,Zc,squeeze(Q_mean(:,center,:)),-cmax:cinta:cmax,'Edgecolor','none'); hold on;
contour(Xc,Zc,squeeze(tend(:,center,:)),[cint:cint:cint*ncontours],'r-','Linewidth',1.2); hold on;
contour(Xc,Zc,squeeze(tend(:,center,:)),[-cint*ncontours:cint:-cint],'b-','Linewidth',1.2);
colorbar
%caxis([-max(max(squeeze(Q_mean(:,center,:)))) max(max(squeeze(Q_mean(:,center,:))))])
caxis([-cmax cmax])
%colormap(redblue(15))
colormap(redblue(nint));
title('\rm \epsilon=0.1');
set(gca,'FontSize',12)
xlabel('x')
ylabel('z')

subplot(1,3,3)

%load('DRV_storms_eps04_w_composite.mat')
load('DRV_storms_eps04_q_composite.mat')

eps = 0.4;

%tend = Q_diab_z_mean;
tend = Q_diab_z_mean + Q_adv_z_mean;
%tend = Q_diab_hoskins_mean;
%tend = theta_prime_ensemble;

%ncontours = 12;
%ncontours = 12;
cint = max(max(tend(:,center,:)))/ncontours;
%cint = 120/ncontours;
disp(cint)

%nint = 55;
%nint = 25;
%cmax = max(max(abs(Q_mean(:,center,:))));
%cmax = 15;
cmax = 4;
%cmax = 19;
cinta = cmax/nint;

contourf(Xc,Zc,squeeze(Q_mean(:,center,:)),-cmax:cinta:cmax,'Edgecolor','none'); hold on;
contour(Xc,Zc,squeeze(tend(:,center,:)),[cint:cint:cint*ncontours],'r-','Linewidth',1.2); hold on;
contour(Xc,Zc,squeeze(tend(:,center,:)),[-cint*ncontours:cint:-cint],'b-','Linewidth',1.2);
colorbar
%caxis([-max(max(squeeze(Q_mean(:,center,:)))) max(max(squeeze(Q_mean(:,center,:))))])
caxis([-cmax cmax])
%colormap(redblue(15))
colormap(redblue(nint));
title('\rm \epsilon=0.4');
set(gca,'FontSize',12)
xlabel('x')
ylabel('z')



annotation('textbox',[0.07, 0.94, 0.01,0.03],'String',"(a)",'EdgeColor','none','FontSize',12);
annotation('textbox',[0.36, 0.94, 0.01,0.03],'String',"(b)",'EdgeColor','none','FontSize',12);
annotation('textbox',[0.64, 0.94, 0.01,0.03],'String',"(c)",'EdgeColor','none','FontSize',12);



% path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
% 'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
%    'primitive_equations_constant_strat_r001_eps04_geostrophic/figures_drv_world/drv_composite_eps04_eps001_only_diab.eps'];
%  
% saveas(gca,path_save,'epsc');