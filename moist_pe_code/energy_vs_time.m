% plot the energy in the Rossby number simulations

list = {'eps04','eps01','eps001'};

en_total = [];
t_total = [];

for ii = 1:3
    
    if ii ==1
    
path_root = ['/nobackup1c/users/mkohl/linear_gradient_formulation/primitive_equations_constant_strat_r001_',list{ii},'_final_high_temp_output_with_advection_more_relaxation'];

    elseif ii == 2
     
path_root = ['/nobackup1c/users/mkohl/linear_gradient_formulation/primitive_equations_constant_strat_r001_',list{ii},'_final_high_temp_output_with_advection_more_relaxation'];

    else
        
path_root = ['/nobackup1c/users/mkohl/linear_gradient_formulation/primitive_equations_constant_strat_r001_',list{ii},'_final_high_temp_output_with_advection'];

    end

path = [path_root,'/snapshots/snapshots_s1.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nz = h5_dims(4); nt = h5_dims(1);

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

w = h5read(path,'/tasks/w');

en = squeeze(mean(mean(mean(w.^2,1),2),3));

plot(t,en); hold on;
set(gca,'FontSize',12)
xlabel('time t')
ylabel('Energe w^2')

end

legend('eps=0.4','eps=0.1','eps=0.01'); legend box off;