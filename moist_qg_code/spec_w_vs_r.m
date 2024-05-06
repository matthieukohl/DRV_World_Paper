% Asymmetry Analysis
close all; clear;

% plot pauls metric

close all; clear;

% loop through all QG runs

list = {'r001','r003','r005','r01_more_hnu','r02','r04_less_hnu','r06','r08','r1'};

r_factor = [0.01,0.03,0.05,0.1,0.2,0.4,0.6,0.8,1.0];

kx_centroid_vs_r = zeros(length(r_factor),1);
ky_centroid_vs_r = zeros(length(r_factor),1);

for kk = 1:length(list)
    
part = '/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/full_relax_correct/equations_correct_nd/qg_final/qg_';

path = [part,list{kk},'/snapshots/snapshots_s1.h5'];



fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta_bc');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); ny = h5_dims(3); nt = h5_dims(1);

% read in variables
t = h5read(path,'/scales/sim_time');

sample = nt-10:1:nt;

L = 12*pi; x = linspace(0,L,nx);
[X,Y] = meshgrid(x,x);


rhs = h5read(path,'/tasks/w');

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
kx_centroid = sum(kx(1:nx/2).*power_rhs_kx(1:nx/2)./sum(power_rhs_kx(1:nx/2)));

disp('ky=')
ky_centroid = sum(ky(1:ny/2)'.*power_rhs_ky(1:ny/2)./sum(power_rhs_ky(1:ny/2)));



kx_centroid_vs_r(kk) = kx_centroid;
ky_centroid_vs_r(kk) = ky_centroid;

end


figure
semilogx(r_factor,kx_centroid_vs_r);
xlabel('Reduction factor r')
ylabel('Wavenumber k')



