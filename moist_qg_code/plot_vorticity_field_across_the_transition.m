% plot vorticity field across the transition point

close all; clear;


list = {'r1','r05','r04','r03'};
r_factor = [1,0.5,0.4,0.3];

% list = {'r03','r035','r04','r05'};
% r_factor = [0.3,0.35,0.4,0.5];

% % 
%list = {'r03','r035','r04','r05','r06','r08','r1'};
%r_factor = [0.3,0.35,0.4,0.5,0.6,0.8,1];

% list = flip(list);
% r_factor = flip(r_factor);

colorint = 100;
cmlev = 21;

   
hfig = figure('Position', [10 10 1000 700]);

for ii = 1:length(list)
    
part = '/nobackup1c/users/mkohl/qg_';

path = [part,list{ii},'_alpha0p15_linear_tau_damping/snapshots/snapshots_s1.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/w');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nt = h5_dims(1);

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

% read in the PV and tendency fields

start = [1 1 nt-7];
count = [nx nx 1];
stride = [1 1 1];

zeta_bc = h5read(path,'/tasks/zeta_bc',start,count,stride);
zeta_bt = h5read(path,'/tasks/zeta_bt',start,count,stride);
w = h5read(path,'/tasks/w',start,count,stride);
phi = h5read(path,'/tasks/tau',start,count,stride);
tau = h5read(path,'/tasks/tau',start,count,stride);


psi1 = phi+tau;
zeta1 = zeta_bt+zeta_bc; 
zeta2 = zeta_bt-zeta_bc;
q1 = zeta1 - tau;
q2 = zeta2 + tau;
beta = 0.78;

tt = nt;

ax1 = subplot(2,2,ii);

var = zeta1;

max_variable(ii) = max(var(:));
skew_variable(ii) = Skew(var);

%cmax = max(max(squeeze(var)));
if ismember(ii,[1,2]) == 1
cmax = 10;
else
cmax = 20;
end
%cmax = 5;
nint = 10;
cint = cmax/nint;

imshow(var,'XData',linspace(0,12*pi,512),'Ydata',linspace(0,12*pi,512));
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)


%colormap(ax1,redblue(nint)); 
colorbar; 
caxis([-cmax cmax]); 

title(['\rm r=',num2str(r_factor(ii))]); 

xlabel('x'); ylabel('y')
%t1 = title('\rm Vorticity (Upper Layer)');
%tpos = get(t1,'position');
%set(t1,'position',tpos + [0 dispy 0 ]);
axis on
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
set(gca,'Ydir','normal')


end

set(hfig,'renderer','Painters');

% saveas(gcf,'qg_r1_vs_qg_r001_snapshots','epsc');
% 
% 
saveas(gcf,['/nfs/pool002/users/mkohl/run_correct/'...
'test_runs/adv_correct/full_relax_correct/equations_correct_nd/'...
'primitive_equations_final/primitive_equations_constant_strat_r001_eps04_geostrophic/'...
'figures_paper/qg_simulations_drv_world_transition_alpha_0p15'],'epsc')
