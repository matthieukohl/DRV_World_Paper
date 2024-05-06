% Asymmetry Analysis
close all; clear;


path = [pwd,'/snapshots/snapshots_s1.h5'];

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); ny = h5_dims(3); nt = h5_dims(1);

% read in variables
t = h5read(path,'/scales/sim_time');

sample = nt-10:1:nt;

L = 12*pi; x = linspace(0,L,nx);
[X,Y] = meshgrid(x,x);


rhs = h5read(path,'/tasks/rhs');
rhs = squeeze(rhs(5,:,:,:));

% plot spectrum of rhs
%x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x); RHS_approx = sin(X) + sin(5*X);
kx = 2*pi/(6*pi) * [0:nx/2,-nx/2+1:-1]; kx = kx/(2*sqrt(2));
ky = 2*pi/(6*pi)*[0:ny/2,-ny/2+1:-1]; ky = ky/(2*sqrt(2));
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

figure('Position', [10 10 1200 500]);

subplot(1,2,1)
plot(kx(1:nx/2),power_rhs_kx(1:nx/2)/max(power_rhs_kx(1:nx/2)))
ylabel('abs(fft(rhs))^2')
xlabel('kx')
xlim([kx(1) kx(nx/2)])
title('1D-Spectrum lat average')

subplot(1,2,2)
plot(ky(1:ny/2),power_rhs_ky(1:ny/2)/max(power_rhs_ky(1:ny/2)))
ylabel('abs(fft(rhs))^2')
xlabel('ky')
xlim([ky(1) ky(ny/2)])
title('1D-Spectrum lat average')

% Make plots of the pdf over the RHS

rhs = rhs(:,:,end);
%rhs = rhs/max(max(abs(rhs)));

figure(2)
histogram(rhs,1000,'Normalization','pdf')
skew = skewness(rhs(:)); 
kurt = kurtosis(rhs(:));
stdd = std(rhs(:));
title(['std=',num2str(round(stdd,1)),'  skew=',num2str(round(skew,1)),'  kurt=',num2str(round(kurt,1))])
%xlimit([min(min(rhs)) -min(min(rhs))])

% test toy_model_predictions
NN = 1e3;
x = linspace(0,2*pi,NN);
forc = pearsrnd(0,stdd,skew,kurt,NN,1);
forc = pearsrnd(0,1,0,3,NN,1);
%forc = rhs(:); forc = forc(1000:1e3+999); forc = forc-mean(forc);

R = 0.01;

[lambda_w,skew_w] = toy_model_output_arb_rhs(R,forc,x);
% 
% for tt = 1:100
%   forc = pearsrnd(0,1,0,10,NN,1);
% %forc = rhs(:); forc = forc(1000:1e3+999); forc = forc-mean(forc);
% 
% R = 0.01;
% 
% [lambda_w(tt),skew_w(tt)] = toy_model_output_arb_rhs(R,forc,x);  
% 
% end

