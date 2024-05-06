close all; clear;

list = {'eps0005','eps001','eps005','eps01','eps04','eps06'}; %,'eps05','eps06'};

list_end = {'snapshots.h5','snapshots.h5','snapshots_s1.h5','snapshots_s1.h5','snapshots_s1.h5','snapshots_s1.h5'};

eps = [0.005,0.01,0.05,0.1,0.4,0.6];

lambda_vs_eps = zeros(length(list),1);

path_pe = '/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/primitive_equations_constant_strat_r001_';


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
semilogx(eps(2:end),lambda_vs_eps(2:end),'b-o','Linewidth',1.6); 
xlabel('Rossby number $\epsilon$'); ylabel('Asymmetry parameter $\lambda$')
set(gca,'FontSize',12);
set(gca,'box','off');

path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
   'primitive_equations_constant_strat_r001_eps04_geostrophic/figures_paper/lambda_vs_eps'];
 
saveas(gca,path_save,'epsc');