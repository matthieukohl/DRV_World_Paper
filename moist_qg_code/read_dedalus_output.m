% Plot dedalus output for the tilted/untilted model
close all; clear;

index = 0;
counter = 0;

for ii = 1:1

path = [pwd,'/snapshots/snapshots_s',num2str(ii),'.h5'];


% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta_bc');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nt = h5_dims(1);

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

hfig = figure('Position', [10 10 1200 500]);
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 2 1])
% 

v = VideoWriter([pwd,'/video/run_r001_res512_beta0.78_cfl2.avi']);
v.FrameRate = 1;
open(v);

for tt = 1:1:nt
    
index = index + 1;
    
start = [1 1 tt];
count = [nx nx 1];
 
zeta_bc = h5read(path,'/tasks/zeta_bc',start,count);
zeta_bt = h5read(path,'/tasks/zeta_bt',start,count);
zeta1 = zeta_bt+zeta_bc; 
zeta2 = zeta_bt-zeta_bc;
% zeta1 = h5read(path,'/tasks/u1',start,count);
% zeta2 = h5read(path,'/tasks/v1',start,count);
w = h5read(path,'/tasks/w',start,count);

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

subplot(1,3,1)
imshow(zeta1); colorbar; caxis([-max(max(zeta1)) max(max(zeta1))])
colormap(redblue)
xlabel('x'); ylabel('y')
caxis([-15 15])
title('\rm Vorticity (Upper Layer)')
axis on
set(gca,'Ydir','normal')
% 
subplot(1,3,2)
imshow(zeta2); colorbar; caxis([-max(max(zeta2)) max(max(zeta2))])
xlabel('x'); ylabel('y')
title('\rm Vorticity (Lower Layer)')
caxis([-15 15])
colormap(redblue)
axis on
set(gca,'Ydir','normal')

subplot(1,3,3)
imshow(w); colorbar; caxis([-max(max(w)) max(max(w))])
colormap(redblue)
xlabel('x'); ylabel('y')
caxis([-50 50])
%title(['time=',num2str(t(tt))])
title('\rm Vertical velocity w')
axis on
set(gca,'Ydir','normal')

%print([pwd,'/images_movie/frame',num2str(index)],'-dpng','-r300');
%print([pwd,'/images_movie/frame',num2str(index)]);

counter = counter +1;
saveas(gcf,[pwd,'/images_movie/frame_labels_',num2str(counter)],'png');



frame = getframe(gcf);
writeVideo(v,frame);

t(tt)


pause(0.000001)


end

end

close(v)
close all;

zeta_bc = h5read(path,'/tasks/zeta_bc'); zeta_bt = h5read(path,'/tasks/zeta_bt');
zeta1 = zeta_bt + zeta_bc; zeta2 = zeta_bt-zeta_bc;
plot(t,Skew(zeta1)); hold on; 
plot(t,Skew(zeta2)); 
xlabel('t'); ylabel('Skew'); 
legend('\zeta_1','\zeta_2'); 
legend boxoff; 
title('\rm h=0');



w = h5read(path,'/tasks/w');

plot(t,Lambda(w))

list = {'r001','r005','r01','r02','r04','r06','r08','r1'};
r = [0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0];

% load QG

for ii = 6:8

path_1 = '/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/full_relax_correct/equations_correct_nd/qg_final/';
path_2 = ['qg_',list{ii}];

path = [path_1,path_2,'/time_series/time_series_s1.h5'];

% find the dimension of variables


t = h5read(path,'/scales/sim_time');

en = h5read(path,'/tasks/en');

figure(10)
semilogy(t,squeeze(en)); hold on;
legend



end


for tt = 1:nt
    
start = [1 1 tt];
count = [nx nx 1];

u1 = h5read(path,'/tasks/u1',start,count);
u2 = h5read(path,'/tasks/u2',start,count);
v1 = h5read(path,'/tasks/v1',start,count);
v2 = h5read(path,'/tasks/v2',start,count);

w = h5read(path,'/tasks/w',start,count);

maxu1(tt) = max(w(:))./max(u1(:));
maxu2(tt) = max(w(:))./max(u2(:));
maxv1(tt) = max(w(:))./max(v1(:));
maxv2(tt) = max(w(:))./max(v2(:));
end

plot(t,maxu1); hold on; 
plot(t,maxu2); hold on;
plot(t,maxv1); hold on;
plot(t,maxv2);
xlabel('t')
legend('max(w)/max(u1)','max(w)/max(u2)','max(w)/max(v1)','max(w)/max(v2)');
legend box off

