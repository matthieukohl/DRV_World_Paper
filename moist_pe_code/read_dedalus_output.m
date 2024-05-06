% Plot dedalus output for the tilted/untilted model
close all; clear;

counter = 0;

% hfig = figure('Position', [10 10 800 800]);
% axes
% axis('equal')

hfig = figure('Position', [10 10 1200 500]);
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 2 1])
% 

for ii = 1:1

path = [pwd,'/snapshots/snapshots_s',num2str(ii),'.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nz = h5_dims(4); nt = h5_dims(1);

x = linspace(0,6*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

%hfig = figure('Position', [10 10 1200 500]);
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 2 1])
% 

v = VideoWriter([pwd,'/video/pe_run_r001_eps04.avi']);
%v.FrameRate = 2;
open(v);

for tt = 1:1:nt
 
start_low = [3 1 1 tt];
count_low = [1 nx nx 1];
 
start_mid = [5 1 1 tt];
count_mid = [1 nx nx 1];

start_up = [7 1 1 tt];
count_up = [1 nx nx 1];

zeta1 = h5read(path,'/tasks/zeta',start_up,count_up); zeta1 = squeeze(zeta1);
zeta2 = h5read(path,'/tasks/zeta',start_low,count_low); zeta2 = squeeze(zeta2);
w = h5read(path,'/tasks/w',start_mid,count_mid); w = squeeze(w);

% subplot(1,3,1)
% contourf(X,Y,zeta1,'edgecolor','none'); colorbar; caxis([-max(max(zeta1)) max(max(zeta1))])
% colormap(redblue)
% xlabel('x'); ylabel('y')
% title('Vorticity (Upper Layer)')
% % 
% subplot(1,3,2)
% contourf(X,Y,zeta2,'edgecolor','none'); colorbar; caxis([-max(max(zeta2)) max(max(zeta2))])
% xlabel('x'); ylabel('y')
% title('Vorticity (Lower Layer)')
% colormap(redblue)
% 
% subplot(1,3,3)
% contourf(X,Y,w,'edgecolor','none'); colorbar; caxis([-max(max(w)) max(max(w))])
% colormap(redblue)
% xlabel('x'); ylabel('y')
% title('w')

% axesHandles = findobj(get(hfig,'Children'), 'flat','Type','axes');
% %Set the axis property to square
% axis(axesHandles,'square')

L = 6*pi;

s1 = subplot(1,3,1);
imshow(zeta1,'XData',linspace(0,L,nx),'Ydata',linspace(0,L,nx)); colorbar; 
%contourf(X,Y,zeta,'Edgecolor','none');
%caxis([-max(max(zeta)) max(max(zeta))])
colormap(redblue)
xlabel('x'); ylabel('y')
title('\rm Vorticity (Upper)')
axis on
set(gca,'Ydir','normal')
caxis([-10 10])
set(gca,'FontSize',12)
xticks([0 3*pi 6*pi])
xticklabels({'0','3\pi','6\pi'})
yticks([0 3*pi 6*pi])
yticklabels({'0','3\pi','6\pi'})

s2 = subplot(1,3,2);
imshow(zeta2,'XData',linspace(0,L,nx),'Ydata',linspace(0,L,nx)); colorbar; 
%contourf(X,Y,zeta,'Edgecolor','none');
%caxis([-max(max(zeta)) max(max(zeta))])
colormap(redblue)
xlabel('x'); ylabel('y')
title('\rm Vorticity (Lower)')
axis on
set(gca,'Ydir','normal')
caxis([-10 10])
set(gca,'FontSize',12)
xticks([0 3*pi 6*pi])
xticklabels({'0','3\pi','6\pi'})
yticks([0 3*pi 6*pi])
yticklabels({'0','3\pi','6\pi'})


s3 = subplot(1,3,3);
imshow(w,'XData',linspace(0,L,nx),'Ydata',linspace(0,L,nx)); 
%contourf(X,Y,w,'Edgecolor','none');
colorbar;
%caxis([-max(max(w)) max(max(w))])
colormap(redblue)
xlabel('x'); ylabel('y')
title('\rm Vertical Velocity (Mid Troposphere)')
axis on
set(gca,'Ydir','normal')
caxis([-3 3])
set(gca,'FontSize',12)
set(gca,'Ydir','normal')
xticks([0 3*pi 6*pi])
xticklabels({'0','3\pi','6\pi'})
yticks([0 3*pi 6*pi])
yticklabels({'0','3\pi','6\pi'})
%title(['\rm',num2str(tt)])

% s1 = subplot(1,3,1);
% 
% cmax = 20;
% nint = 100;
% cint = cmax/nint;
% 
% contourf(X,Y,zeta1,-cmax:cint:cmax,'Edgecolor','none');
% caxis([-cmax cmax])
% colormap(redblue(nint))
% colorbar
% xlabel('x'); ylabel('y')
% title('\rm Vorticity (Upper Level)')
% axis on
% set(gca,'Ydir','normal')
% caxis([-20 20])
% set(gca,'FontSize',12)
% xticks([0 3*pi 6*pi])
% xticklabels({'0','3\pi','6\pi'})
% yticks([0 3*pi 6*pi])
% yticklabels({'0','3\pi','6\pi'})
% %axis equal;
% pbaspect([1 1 1])
% 
% s2 = subplot(1,3,2);
% 
% contourf(X,Y,zeta2,-cmax:cint:cmax,'Edgecolor','none');
% colorbar
% colormap(redblue(nint))
% xlabel('x'); ylabel('y')
% title('\rm Vorticity (Lower Level)')
% axis on
% set(gca,'Ydir','normal')
% caxis([-20 20])
% set(gca,'FontSize',12)
% xticks([0 3*pi 6*pi])
% xticklabels({'0','3\pi','6\pi'})
% yticks([0 3*pi 6*pi])
% yticklabels({'0','3\pi','6\pi'})
% %axis equal;
% pbaspect([1 1 1])
% 
% 
% s3 = subplot(1,3,3);
% 
% cmax = 30;
% nint = 100;
% cint = cmax/nint;
% 
% contourf(X,Y,w,-cmax:cint:cmax,'Edgecolor','none');
% colorbar;
% %caxis([-max(max(w)) max(max(w))])
% colormap(redblue(nint))
% xlabel('x'); ylabel('y')
% title('\rm Vertical Velocity (Mid Level)')
% axis on
% set(gca,'Ydir','normal')
% caxis([-30 30])
% set(gca,'FontSize',12)
% set(gca,'Ydir','normal')
% xticks([0 3*pi 6*pi])
% xticklabels({'0','3\pi','6\pi'})
% yticks([0 3*pi 6*pi])
% yticklabels({'0','3\pi','6\pi'})
% %title(['\rm',num2str(tt)])
% %axis equal;
% pbaspect([1 1 1])





%annotation('textbox',[0.05 0.01 0.9 0.12],'String','\rm Moist primitive equation simulations forced by relaxation to a cosinusoidal temperature profile giving rise to two baroclinic zones. The domain is doubly periodic in the horizontal and bounded by plates in the vertical. In this high Rossby number regime, the flow shows the typical zonal wave-like flow pattern along the baroclinic zones despite strong latent heating.','EdgeColor','none','FontWeight','light');
%annotation('textbox',[0.01 0.4 0.1 0.3],'String','Low Rossby Number Regime','FontSize',12,'EdgeColor','none','FontWeight','light');



%set(s1,'Position',get(s1,'Position')+[0 0.05 0 0]);
%set(s2,'Position',get(s2,'Position')+[0 0.05 0 0]);

counter = counter + 1;

saveas(gcf,[pwd,'/images_movie/frame_',num2str(counter)],'png');


frame = getframe(gcf);
writeVideo(v,frame);


pause(1)


end

end

close(v)
close all

% make a plot of the horizontal stratification

theta_z = h5read(path,'/tasks/theta_z');

z = linspace(0,1,10);

theta_z = squeeze(mean(mean(mean(theta_z(:,:,:,100:end),2),3),4));

plot(1+0.01*theta_z,z,'b--','Linewidth',1.4)
xlabel('1+eps*theta_z')
ylabel('z')
set(gca,'FontSize',12);
title('\rm Stratification (eps=0.01)');


% test

% plot the energy vs time

en_total = [];
t_total = [];

for ii = 1:1

path = [pwd,'/snapshots/snapshots_s',num2str(ii),'.h5'];

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

en_total = [en_total;en];
t_total = [t_total;t];

end

%

figure
semilogy(t_total,en_total);
xlabel('time total')
ylabel('energy total');







for ii = 1:nz

start = [ii 1 1 1];
count = [1 nx nx nt];
 
w = h5read(path,'/tasks/w',start,count); w = squeeze(w);

lambda(ii,:) = Lambda(w);


low_layers = [1,2,3];
up_layers = [8,9,19];

if ii <=3
h(ii) = plot(t,lambda(ii,:),'b'); hold on
elseif ii >=8
h(ii) = plot(t,lambda(ii,:),'r'); hold on;
else
h(ii) = plot(t,lambda(ii,:),'g'); hold on;
end
xlabel('time')
ylabel('\lambda')
end

plot(t,mean(lambda,1),'k-.','linewidth',1.4); hold on;

legend([h(1),h(4),h(10)],'bottom','middle','top'); legend boxoff

figure(20)
plot(t,mean(lambda,1),'b','linewidth',1.4); hold on;
plot(t,mean(lambda_qg,1),'r','linewidth',1.4); hold on;
plot(t,mean(lambda_qg_dry,1),'k','linewidth',1.4);
xlabel('t');
ylabel('\lambda')
set(gca,'FontSize',14)
legend('PE','QG','QG dry'); 
legend boxoff
title('\epsilon = 0.4')
ylim([0.5 1])
set(gca,'box','off')
% path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
% 'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
%     'primitive_equations_constant_strat_r001_eps04_geostrophic/figures_paper/lambda_pe_eps04'];
% 
% saveas(gca,path_save,'epsc')

figure(2)
plot(mean(lambda(:,end-10:end),2),linspace(0,1,10));
xlabel('\lambda')
ylabel('z')

w = h5read(path,'/tasks/w');
wz = zeros(size(10,1));
counter = 0;
for ii = 1:nx
    for jj = 1:nx
        for tt = 1:nt
            if w(5,ii,jj,tt)>0
            wz = wz + squeeze(w(:,ii,jj,tt));
            counter = counter + 1;
            end
        end
    end
end

wz = wz/counter;
figure(3)
plot(wz,linspace(0,1,10));

figure(4)
rhs = h5read(path,'/tasks/rhs');
for ii = 1:10
 skew(ii,:) = Skew(squeeze(-rhs(ii,:,:,:)));
if ii <=3
h(ii) = plot(t,skew(ii,:),'b'); hold on
elseif ii >=8
h(ii) = plot(t,skew(ii,:),'r'); hold on;
else
h(ii) = plot(t,skew(ii,:),'g'); hold on;
end
xlabel('time')
ylabel('Skew')
end
plot(t,mean(skew,1),'k-.','linewidth',1.4);
title('rhs=u_z\cdot\nabla\zeta')
legend([h(1),h(4),h(10)],'bottom','middle','top'); legend boxoff

hfig = figure('Position', [10 10 1200 500]);
w = h5read(path,'/tasks/w');
w500 = squeeze(w(5,:,:,end)); hold on;
wqg = h5read(path,'/tasks/wqg');
wqg500 = squeeze(wqg(5,:,:,end)); hold on;

subplot(1,3,1)
imshow(w500); colorbar; caxis([-max(max(w500)) max(max(w500))])
colormap(redblue)
xlabel('x'); ylabel('y')
title('w')

subplot(1,3,2)
imshow(wqg500); colorbar; caxis([-max(max(w500)) max(max(w500))])
colormap(redblue)
xlabel('x'); ylabel('y')
title('wqg')

subplot(1,3,3)
imshow(w500-wqg500); colorbar; caxis([-max(max(w500)) max(max(w500))])
colormap(redblue)
xlabel('x'); ylabel('y')
title('w-wqg')


%
y = linspace(0,6*pi,128);
z = linspace(0,1,10);
[Y,Z] = meshgrid(y,z);
theta_bar = Z/0.4;

theta = h5read(path,'/tasks/theta');
theta_mean = squeeze(mean(mean(theta(:,:,:,end-10:end),3),4));

figure
contourf(Y,Z,theta_mean,'EdgeColor','none'); 
colorbar


theta_test = theta_bar + theta_mean;

theta_test_z = zeros(size(theta_test));
theta_test_z(1:end-1,:) = (theta_test(2:end,:)-theta_test(1:end-1,:))./(Z(2:end,:)-Z(1:end-1,:));

figure
contourf(Y,Z,theta_bar+theta_mean,'EdgeColor','none'); 
colorbar

figure
contourf(Y,Z,theta_test_z,'EdgeColor','none'); 
colorbar

% theta_z = h5read(path,'/tasks/theta_z');
% for tt = 1:size(theta_z,4)
%     test = 1 + 0.01*theta_z(:,:,:,tt);
%     tested = any(test<0);
%     index(tt) = sum(tested(:));
% end
% 
% plot(index)

theta_z = h5read(path,'/tasks/theta_z');
for tt = 1:size(theta_z,4)
    test = 1 + 0.01*theta_z(:,:,:,tt);
    tested = any(test<0);
    index(tt) = sum(tested(:));
end
grid_number = 128*128*10;
plot(index/grid_number)