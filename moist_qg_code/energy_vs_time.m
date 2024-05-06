% plot pauls metric

close all; clear;

% loop through all QG runs

% with linear damping

list = {'r03','r035','r04','r05','r06','r08','r1'};
r_factor = [0.3,0.35,0.4,0.5,0.6,0.8,1];

for ii = 1:length(list)
    
part = '/nobackup1c/users/mkohl/qg_';

t_total = [];
en_total = [];

for i = 1:4
    
path = [part,list{ii},'_alpha0p15_linear_tau_damping/snapshots/snapshots_s',num2str(i),'.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/w');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nt = h5_dims(1);

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

w = h5read(path,'/tasks/w');

u1 = h5read(path,'/tasks/u1');
u2 = h5read(path,'/tasks/u2');
v1 = h5read(path,'/tasks/v1');
v2 = h5read(path,'/tasks/v2');
tau = h5read(path,'/tasks/tau');

%en = squeeze(mean(mean(w.^2,1),2));

t_total = [t_total;t];

en = u1.^2 + u2.^2 + v1.^2 + v2.^2 + tau.^2;
en = squeeze(mean(mean(en,1),2));

en_total = [en_total;en];

end


semilogy(t_total,en_total,'-o'); hold on;

% semilogy(squeeze(mean(mean(w.^2,1),2)),'-o'); hold on;

%plot(Lambda(w),'-o'); hold on;

end

legend('r03','r035','r04','r05','r06','r08','r1'); legend boxoff;