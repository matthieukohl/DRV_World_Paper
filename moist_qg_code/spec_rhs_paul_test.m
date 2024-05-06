% Asymmetry Analysis
close all; clear;


path = [pwd,'/snapshots/snapshots_s1.h5'];

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta_bc');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); ny = h5_dims(3); nt = h5_dims(1);

% read in variables
t = h5read(path,'/scales/sim_time');

sample = nt-7:1:nt;

L = 12*pi; x = linspace(0,L,nx);
[X,Y] = meshgrid(x,x);

list = {'u1','w'};

power_rhs_kx_vs_var = zeros(512,2);
power_rhs_ky_vs_var = zeros(512,2);

for kk = 1:2

if kk == 1
rhs = h5read(path,'/tasks/u1');
else
rhs = h5read(path,'/tasks/w');
end


% plot spectrum of rhs
%x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x); RHS_approx = sin(X) + sin(5*X);
kx = 2*pi/(12*pi) * [0:nx/2,-nx/2+1:-1];
ky = 2*pi/(12*pi)*[0:ny/2,-ny/2+1:-1];
n = [0:nx/2,-nx/2+1:-1];
power_rhs_kx = 0;
power_rhs_ky = 0;
counter = 0;

%RHS_approx = RHS_approx - mean(RHS_approx,1);

% spec - along kx

for ii = 1:size(rhs,1)
   for tt = 1:length(sample)
rhsf = fft(rhs(ii,:,tt));
power = abs(rhsf).^2;
power_rhs_kx = power_rhs_kx + power;
counter = counter + 1;
   end
end
power_rhs_kx = power_rhs_kx/(size(rhs,1)*length(sample));

% spec - along ky
counter = 0;
for ii = 1:size(rhs,2)
   for tt = 1:length(sample)
rhsf = fft(rhs(:,ii,tt));
power = abs(rhsf).^2;
power_rhs_ky = power_rhs_ky + power;
counter = counter + 1;
   end
end
power_rhs_ky = power_rhs_ky/(size(rhs,2)*length(sample));

disp('kx=')
kx_centroid = sum(kx(1:nx/2).*power_rhs_kx(1:nx/2)./sum(power_rhs_kx(1:nx/2)))

disp('ky=')
ky_centroid = sum(ky(1:ny/2)'.*power_rhs_ky(1:ny/2)./sum(power_rhs_ky(1:ny/2)))

power_rhs_kx_vs_var(:,kk) = power_rhs_kx;
power_rhs_ky_vs_var(:,kk) = power_rhs_ky;

end
%figure('renderer','painters','Position',[10 10 1000 400]);

% compare the spectrum of kx^2 * u^2 to w^2

figure

loglog(kx(1:nx/2),kx(1:nx/2).^2.*power_rhs_kx_vs_var(1:nx/2,1)','b'); hold on; %/max(power_rhs_kx(1:nx/2)))
loglog(kx(1:nx/2),power_rhs_kx_vs_var(1:nx/2,2),'r'); hold on; %/max(power_rhs_kx(1:nx/2)))
legend('u1^2*k^2','w'); legend boxoff
ylabel('Magnitude')
xlabel('kx')
xlim([kx(1) kx(nx/2)])
title('\rm 1D-Spectrum x (r=0.01)')

figure

loglog(ky(1:ny/2),ky(1:ny/2).^2.*power_rhs_ky_vs_var(1:ny/2,1)','b'); hold on; %/max(power_rhs_kx(1:nx/2)))
loglog(ky(1:ny/2),power_rhs_kx_vs_var(1:ny/2,2),'r'); hold on; %/max(power_rhs_kx(1:nx/2)))
legend('u1^2*k^2','w^2'); legend boxoff
ylabel('Magnitude')
xlabel('kx')
xlim([ky(1) ky(ny/2)])
title('\rm 1D-Spectrum y')




% subplot(1,2,2)

% semilogx(kx(1:nx/2),(kx(1:nx/2).^2+1).*power_rhs_kx(1:nx/2)); hold on; %/max(power_rhs_kx(1:nx/2)))
% ylabel('Magnitude')
% xlabel('kx')
% xlim([kx(1) kx(nx/2)])
% title('\rm (k^2+1)*1D-Spectrum x')



