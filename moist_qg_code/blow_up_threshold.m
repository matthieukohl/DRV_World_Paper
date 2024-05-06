% plot pauls metric

close all; clear;

% loop through all QG runs

% with linear damping

% list = {'r03','r035','r037','r038','r039','r04','r045'};
% r_factor = [0.3,0.35,0.37,0.38,0.39,0.4,0.45];

list = {'r03','r035','r04','r045'};
r_factor = [0.3,0.35,0.4,0.45];

%t_total = zeros(44,length(r_factor));
%energy = zeros(44,length(r_factor));

t_total = zeros(25,length(r_factor));
energy = zeros(25,length(r_factor));

for ii = 1:length(list)
    
t_ind = [];
energy_ind = [];
 
part = '/nobackup1c/users/mkohl/qg_';

path = [part,list{ii},'_no_linear_tau_damping/snapshots/snapshots_s1.h5'];

if ii==2 || ii==3
    
path = [part,list{ii},'_no_linear_tau_damping_new/snapshots/snapshots_s1.h5'];
 
end

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/w');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nt = h5_dims(1);

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

% read in the w fields

w = h5read(path,'/tasks/w');

u1 = h5read(path,'/tasks/u1');
u2 = h5read(path,'/tasks/u2');
v1 = h5read(path,'/tasks/v1');
v2 = h5read(path,'/tasks/v2');
tau = h5read(path,'/tasks/tau');

%en = squeeze(mean(mean(w.^2,1),2));

en = u1.^2 + u2.^2 + v1.^2 + v2.^2 + tau.^2;
en = squeeze(mean(mean(en,1),2));


t_total(:,ii) = t(1:25);
energy(:,ii) = en(1:25);
[a,index(ii)] = min(abs(t_total(:,ii)-70));

cmap = cool(length(r_factor));
plot(t,en,'color',cmap(ii,:),'Linewidth',1.2); hold on;

end



% figure
% for ii = 1:length(r_factor)
% plot(t_total(:,ii),energy(:,ii),'color',cmap(ii,:),'Linewidth',1.2); hold on;
% end
xlabel('time')
ylabel('Energy')
% legend('r=0.3','r=0.35','r=0.37','r=0.38','r=0.39','r=0.4','r=0.45'); 
%legend('r=0.3','r=0.35','r=0.37','r=0.38','r=0.39','r=0.4','r=0.45','Location','NorthWest','Linewidth',1.5);
legend('r=0.3','r=0.35','r=0.4','r=0.45','Location','NorthWest','Linewidth',1.5);
set(gca,'FontSize',12)
set(gca,'box','off')
legend boxoff;
xlim([0 50])


saveas(gcf,['/nfs/pool002/users/mkohl/run_correct/'...
'test_runs/adv_correct/full_relax_correct/equations_correct_nd/'...
'primitive_equations_final/primitive_equations_constant_strat_r001_eps04_geostrophic/'...
'figures_paper/blow_up_threshold'],'epsc')

figure
plot(t_total(:,2),energy(:,2),'b','Linewidth',1.2); hold on;
plot(t_total(:,6),energy(:,6),'r','Linewidth',1.2); 
set(gca,'FontSize',12)
xlim([0 70]);
legend('r=0.35','r=0.4'); legend boxoff; 
xlabel('time (t)')
ylabel('Energy')




% figure
% plot(t_total,energy); hold on;
% %semilogy(t_total,energy);
% xlabel('time')
% ylabel('energy')
% legend('r=0.3','r=0.35','r=0.37','r=0.38','r=0.39','r=0.4','r=0.45'); 
% %legend('r=0.35','r=0.37','r=0.38','r=0.39','r=0.4');
% legend boxoff;
% xlim([0 70])

figure
semilogx(r_factor,mean(energy,1),'-o','linewidth',1.2);
xlabel('Reduction factor');
ylabel('Mean Energy');

figure

for ii = 1:length(r_factor)

semilogx(r_factor,mean(energy(1:index(ii),:),1),'-o','linewidth',1.2);
xlabel('Reduction factor');
ylabel('Mean Energy');   
   
end



