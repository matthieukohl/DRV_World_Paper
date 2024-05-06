close all; clear;

list = {'eps0005','eps001','eps005','eps01','eps04','eps06'}; %,'eps05','eps06'};

list_end = {'snapshots.h5','snapshots.h5','snapshots_s1.h5','snapshots_s1.h5','snapshots_s1.h5','snapshots_s1.h5'};

eps = [0.005,0.01,0.05,0.1,0.4,0.6];

ratio_vs_eps = zeros(length(list),1);

path_pe = '/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/primitive_equations_constant_strat_r001_';

for kk = 1:length(list)
kk

path = [path_pe,list{kk},'_geostrophic_with_wdry/snapshots/',list_end{kk}];

%path = [path,'/snapshots/snapshots_s1.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nz = h5_dims(4); nt = h5_dims(1);

x = linspace(0,6*pi,nx); dx = x(2)-x(1);
[X,Y] = meshgrid(x,x);

% only select one of the double planets
X = X(80:end,:); Y = Y(80:end,:);

t = h5read(path,'/scales/sim_time');



ratio_vs_time = zeros(length(t),1);

for tt = 1:length(t)
tt

start = [5 1 1 tt];
count = [1 nx nx 1];
 
w = h5read(path,'/tasks/w',start,count); w = squeeze(w);
w = w(80:end,:);

z = w;
%z = smoothdata(smoothdata(z,1),2);
%z = (X-6*pi).^2.*exp(-((X-6*pi).^2+(Y-6*pi).^2));
%z = randn(size(z));

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

% filter out points that are too close together (multiple maxima within
% same contour)

% diff_peak_x = diff(peak_x); diff_peak_y = diff(peak_y);
% diff_peak = 0.5*(diff_peak_x+diff_peak_y);
% cut_off = 10;
% 
% peak_x = peak_x(diff_peak>cut_off);
% peak_y = peak_y(diff_peak>cut_off);
% peak_xwidth = peak_xwidth(diff_peak>cut_off);
% peak_ywidth = peak_ywidth(diff_peak>cut_off);

% Calculate the ratio of peak distance to peak width for the 5 strongest
% ascents

peak_x = peak_x(1:5); peak_y = peak_y(1:5);

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

%plot w and mark the 5 highest peaks
% figure
% z = z/max(max(z));
% contourf(X,Y,z); hold on;
% caxis([-max(max(z)) max(max(z))])
% colorbar
% colormap(redblue)
% scatter(peak_x,peak_y,'+g');

% hold on

peak_width = 0.5*(peak_xwidth+peak_ywidth)*dx;
%peak_width = 0.5*(peak_xwidth(1:5)+peak_ywidth(1:5))*dx; % convert from grid points to length
ratio = 0;
counter = 0;

for ii = 1:length(peak_x)
    for jj = ii:length(peak_y)
     if ii==jj
         continue % exclude distances from point to itself
     else
    distance = sqrt((peak_x(ii)-peak_x(jj))^2+(peak_y(ii)-peak_y(jj)).^2);
    if distance<1
        continue
    end
    ratio = ratio + distance/peak_width(ii);
    counter = counter + 1;
     end
    end
end

if counter ==0
    ratio = NaN;
else
ratio = ratio/counter;
end


ratio_vs_time(tt) = ratio;
    
end

% plot ratio vs. time

plot(t,ratio_vs_time,'Linewidth',1.6); hold on;
xlabel('time');
ylabel('distance/width');
set(gca,'FontSize',12);
%legend('\eps=0.0005','\eps=0.005','\epsilon =');
legend boxoff

%ratio_vs_r(kk) = mean(ratio_vs_time(end-10:end));
ratio_vs_eps(kk) = nanmean(ratio_vs_time(20:end));


end

figure(2)
semilogx(eps,ratio_vs_eps,'Linewidth',1.6); 
xlabel('r'); ylabel('ratio')
set(gca,'FontSize',12);

for kk = 1:length(list)
kk

path = [path_pe,list{kk},'_geostrophic_with_wdry/snapshots/',list_end{kk}];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nz = h5_dims(4); nt = h5_dims(1);

% start = [5 1 1 1];
% count = [1 nx nx nt];
%  
% w = h5read(path,'/tasks/w',start,count); w = squeeze(w);

w = h5read(path,'/tasks/w');

lambda = 0;

for ii = 1:nz
lambda = lambda + mean(Lambda(squeeze(w(ii,:,:,end-10:end))));
end

lambda_vs_eps(kk) = lambda/nz;

end

figure(3)
semilogx(eps(2:end),lambda_vs_eps,'b-o','Linewidth',1.6); 
xlabel('\epsilon'); ylabel('Asymmetry parameter \lambda')
set(gca,'FontSize',12);