% compare the no./prominence of w-maxima across qg r=0.01 simulations with
% different pv gradients

close all; clear;

list = {'eps0005','eps001','eps005','eps01','eps04','eps06'}; %,'eps05','eps06'};

list_end = {'snapshots.h5','snapshots.h5','snapshots_s1.h5','snapshots_s1.h5','snapshots_s1.h5','snapshots_s1.h5'};

eps = [0.005,0.01,0.05,0.1,0.4,0.6];

ratio_vs_eps = zeros(length(list),1);

path_pe = '/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/primitive_equations_constant_strat_r001_';

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

ratio_count = zeros(length(t),1);

for tt = 1:length(t)

start = [5 1 1 tt];
count = [1 nx nx 1];
 
w = h5read(path,'/tasks/w',start,count); w = squeeze(w);
w = w(80:end,:);

ratio = 0;
counter = 0;

for kk = 1:size(w,1)
 

    
wslice = squeeze(w(kk,:));
wslice = wslice/max(wslice);

% if kk==35
%     phi = 5;
% end

[pks,locs,width,p] = findpeaks(wslice,'MinPeakHeight',0.5);

[pks,I] = sort(pks,'descend');
[locs] = locs(I);
[width] = width(I);

if length(locs)==0
    continue
elseif length(locs)==1
ratio = ratio+nx/width;
counter = counter + 1;
else

for i = 1:length(locs)
    for j = i:length(locs)
        if i==j
            continue
        end
       distance = min(abs(locs(i)-locs(j)),nx-abs(locs(i)-locs(j)));
       ratio = ratio + distance/width(i);
       counter = counter + 1;
    end
end
end


end


 

ratio_count(tt) = ratio/counter;
 
end

% plot 1d for each simulation

figure(1)
plot(t,ratio_count,'Linewidth',1.2); hold on;
xlabel('t'); ylabel('1d Ratio');
%legend('h0 with \beta','h0','h05','h1','h2','h5'); legend boxoff;
set(gca,'FontSize',12);
%ylim([0 0.5])

ratio_vs_eps(ii) = mean(ratio_count(10:end));



end

figure(5)
semilogx(eps,ratio_vs_eps,'b-o','Linewidth',1.2);
xlabel('\epsilon'); ylabel('1d Ratio')
set(gca,'FontSize',12)