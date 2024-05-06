% Plot dedalus output for the tilted/untilted model
close all; clear;

% % r=0.01
 
path = [pwd,'/flux_spec.h5'];

% % r =1
% 
% path = '/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/full_relax_correct/equations_correct_nd/qg_final/qg_r1';
% 
% path = [path,'/flux_spec.h5'];

t = h5read(path,'/spectra/time');
modn = h5read(path,'/spectra/modn');
modn = double(modn);

spec_ke_phi = h5read(path,'/spectra/KE_phi');
spec_ke_tau = h5read(path,'/spectra/KE_tau');
spec_ke_w = h5read(path,'/spectra/KE_w');

% average the spectra 

%ind = length(t);
ind = 100;
%ind = 2;

spec_ke_phi = mean(spec_ke_phi(:,ind:end),2);
spec_ke_tau = mean(spec_ke_tau(:,ind:end),2);
spec_ke_w = mean(spec_ke_w(:,ind:end),2);

% figure
% loglog(modn,spec_ke_phi,'b'); hold on;
% loglog(modn,spec_ke_tau,'r'); hold on;
% loglog(modn,spec_ke_w,'k'); hold on;
% loglog(modn,1e2*(modn.^-3),'k--');
% xlabel('Wavenumber')
% ylabel('Spectrum');
% legend('KE Phi','KE Tau','KE w','k^{-3}'); legend boxoff;
% title(['\rm Energy Spectra r=0.01',' t=',num2str(round(t(ind),2))])
% set(gca,'FontSize',12);

figure
loglog(modn,spec_ke_phi,'b'); hold on;
loglog(modn,spec_ke_tau,'r'); hold on;
loglog(modn,spec_ke_w,'k'); hold on;
loglog(modn,1e2*(modn.^-3),'k--');
xlabel('Wavenumber')
ylabel('Spectrum');
legend('KE Phi','KE Tau','KE w','k^{-3}'); legend boxoff;
title(['\rm Energy Spectra r=0.01 (linear damping)'])
set(gca,'FontSize',12);



