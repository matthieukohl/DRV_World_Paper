% compare the no./prominence of w-maxima across qg r=0.01 simulations with
% different pv gradients

close all; clear;

list = {'eps0005','eps001','eps005','eps04','eps06'}; %,'eps05','eps06'};

list_end = {'snapshots.h5','snapshots.h5','snapshots_s1.h5','snapshots_s1.h5','snapshots_s1.h5'};

eps = [0.005,0.01,0.05,0.4,0.6];

ratio_vs_eps = zeros(length(list),1);

path_pe = '/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/primitive_equations_constant_strat_r001_';



p_peaks_vs_eps = zeros(length(list),1);

for ii = 1:length(list)
ii

path = [path_pe,list{ii},'_geostrophic_with_wdry/snapshots/',list_end{ii}];

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

no_peaks_count = zeros(length(t),1);
p_peaks_count = zeros(length(t),1);

for tt = 1:length(t)
no_peaks = 0;
p_peaks = 0;

start = [5 1 1 tt];
count = [1 nx nx 1];
 
w = h5read(path,'/tasks/w',start,count); w = squeeze(w);
w = w(80:end,:);

for kk = 1:size(w,1)
    
wslice = squeeze(w(kk,:));
wslice = wslice/norm(wslice);

[pks,locs,width,p] = findpeaks(wslice,'MinPeakHeight',0.1);

no_peaks = no_peaks + length(locs);
p_peaks = p_peaks + sum(p);

 end

% if tt==14
%     break
% end
 
no_peaks_count(tt) = no_peaks;
p_peaks_count(tt) = p_peaks/no_peaks;
 
end

% plot the no. and prominence for each simulations

% figure(1)
% plot(t,no_peaks_count,'Linewidth',1.2); hold on;
% xlabel('t'); ylabel('number of peaks')
% %legend('h0 with \beta','h0','h05','h1','h2','h5','Location','SouthEast'); legend boxoff
% set(gca,'FontSize',12);
% 
figure(2)
plot(t,p_peaks_count,'Linewidth',1.2); hold on;
xlabel('t'); ylabel('prominence');
%legend('h0 with \beta','h0','h05','h1','h2','h5'); legend boxoff;
set(gca,'FontSize',12);
ylim([0 0.5])

p_peaks_vs_eps(ii) = mean(p_peaks_count(10:end));



end

figure(5)
semilogx(eps,p_peaks_vs_eps,'Linewidth',1.2);
xlabel('\epsilon'); ylabel('prominence')
set(gca,'FontSize',12)

