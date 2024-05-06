% Plot the qg_r1 and qg_r001 fields 
close all; clear;

path = ['/nobackup1c/users/mkohl/qg_r1/snapshots/snapshots_s1.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta_bc');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nt = h5_dims(1);

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

%fig = figure('Position', [10 10 1000 1200]);
fig = figure('Position', [10 10 1000 1200]);

% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 2 1])
% 

v = VideoWriter([pwd,'/video/test.mp4','MPEG-4']);
%v.FrameRate = 2;
open(v);

colorint = 100;
cmlev = 21;
dispy = 0;

mag = 1;

for tt = nt:5:nt
tt    
 
start = [1 1 tt];
count = [nx/mag nx/mag 1];
 
zeta_bc = h5read(path,'/tasks/zeta_bc',start,count);
zeta_bt = h5read(path,'/tasks/zeta_bt',start,count);
phi = h5read(path,'/tasks/phi',start,count);
tau = h5read(path,'/tasks/tau',start,count);
psi1 = phi+tau;
psi2 = phi-tau;
zeta1 = zeta_bt+zeta_bc; 
zeta2 = zeta_bt-zeta_bc;
q1 = zeta1 - tau;
q2 = zeta2 + tau;
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

subplot(3,2,1)

cmax = 20;
nint = 100;
cint = cmax/nint;

%contourf(X,Y,zeta1,-cmax:cint:cmax,'Edgecolor','none');
imshow(zeta1,'XData',linspace(0,12*pi/mag,512),'Ydata',linspace(0,12*pi/mag,512));
caxis([-cmax cmax])
colormap(redblue(nint));
colorbar

xlabel('x'); ylabel('y')
axis on
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
set(gca,'Ydir','normal')


subplot(3,2,3)

cmax = 6;
nint = 100;
cint = cmax/nint;

%contourf(X,Y,zeta2,-cmax:cint:cmax,'Edgecolor','none');
imshow(zeta2,'XData',linspace(0,12*pi/mag,512),'Ydata',linspace(0,12*pi/mag,512));
caxis([-cmax cmax])
colormap(redblue(nint));
colorbar
xlabel('x'); ylabel('y')
axis on
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
set(gca,'Ydir','normal')
caxis([-6 6])

subplot(3,2,5)

cmax = max(abs(w(:)));
nint = 100;
cint = cmax/nint;

% contourf(X,Y,w,-cmax:cint:cmax,'Edgecolor','none');
imshow(w,'XData',linspace(0,12*pi/mag,512),'Ydata',linspace(0,12*pi/mag,512))
caxis([-cmax cmax])
colormap(redblue(nint));

cbr = colorbar;
set(cbr,'YTick',-8:4:8);
xlabel('x'); ylabel('y')

axis on
set(gca,'Ydir','normal')
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})

%print([pwd,'/images_movie/frame',num2str(tt)],'-dpng','-r300');

frame = getframe(gcf);
writeVideo(v,frame);


pause(0.000001)


end

path = [pwd,'/snapshots/snapshots_s1.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta_bc');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nt = h5_dims(1);

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

v = VideoWriter([pwd,'/video/test.mp4','MPEG-4']);
%v.FrameRate = 2;
open(v);

mag = 1;
magw = 1;

for tt = nt:5:nt
tt    
 
start = [1 1 tt];
count = [nx/mag nx/mag 1];

start_w = [1 1 tt];
count_w = [nx/magw nx/magw 1];
 
zeta_bc = h5read(path,'/tasks/zeta_bc',start,count);
zeta_bt = h5read(path,'/tasks/zeta_bt',start,count);
phi = h5read(path,'/tasks/phi',start,count);
tau = h5read(path,'/tasks/tau',start,count);
psi1 = phi+tau;
psi2 = phi-tau;
zeta1 = zeta_bt+zeta_bc; 
zeta2 = zeta_bt-zeta_bc;
q1 = zeta1 - tau;
q2 = zeta2 + tau;
w = h5read(path,'/tasks/w',start_w,count_w);

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

subplot(3,2,2)

cmax = 15;
nint = 100;
cint = cmax/nint;

%contourf(X,Y,zeta1,-cmax:cint:cmax,'Edgecolor','none');
imshow(zeta1,'XData',linspace(0,12*pi/mag,512),'Ydata',linspace(0,12*pi/mag,512))
caxis([-cmax cmax])
colormap(redblue(nint));
colorbar

xlabel('x'); ylabel('y')
%title('\rm Vorticity (Upper Layer)')
axis on
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
set(gca,'Ydir','normal')

% 
subplot(3,2,4)

cmax = 15;
nint = 100;
cint = cmax/nint;

%contourf(X,Y,zeta2,-cmax:cint:cmax,'Edgecolor','none');
imshow(zeta2,'XData',linspace(0,12*pi/mag,512),'Ydata',linspace(0,12*pi/mag,512))
caxis([-cmax cmax])
colormap(redblue(nint));
colorbar

xlabel('x'); ylabel('y')

axis on
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
%title('\rm Vorticity (Lower Layer)')
set(gca,'Ydir','normal')

subplot(3,2,6)

%cmax = max(abs(w(:)));
cmax  = 100;
nint = 100;
cint = cmax/nint;

%contourf(X,Y,w,-cmax:cint:cmax,'Edgecolor','none');
imshow(w,'XData',linspace(0,12*pi/mag,512),'Ydata',linspace(0,12*pi/mag,512));
caxis([-cmax cmax])
colormap(redblue(nint));
colorbar

set(cbr,'YTick',-140:70:140);

xlabel('x'); ylabel('y')
%title('\rm Vertical Velocity')
axis on
set(gca,'Ydir','normal')

if magw ==1
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})

else

xticks([0 3*pi 6*pi])
xticklabels({'0','3\pi','6\pi'})
yticks([0 3*pi 6*pi])
yticklabels({'0','3\pi','6\pi'})

end

%caxis([-140 140])
%caxis([-max(abs(w(:))) max(abs(w(:)))])

%print([pwd,'/images_movie/frame',num2str(tt)],'-dpng','-r300');


frame = getframe(gcf);
writeVideo(v,frame);


pause(0.000001)


end

% a,b,c,d

annotation('textbox',[0.07, 0.91, 0.01,0.03],'String',"(a)",'EdgeColor','none');
annotation('textbox',[0.51, 0.91, 0.01,0.03],'String',"(b)",'EdgeColor','none');

annotation('textbox',[0.07, 0.61, 0.01,0.03],'String',"(c)",'EdgeColor','none');
annotation('textbox',[0.51, 0.61, 0.01,0.03],'String',"(d)",'EdgeColor','none');

annotation('textbox',[0.07, 0.31, 0.01,0.03],'String',"(e)",'EdgeColor','none');
annotation('textbox',[0.51, 0.31, 0.01,0.03],'String',"(f)",'EdgeColor','none');


% textboxes
annotation('textbox',[0.89,0.83,0.15,0.01],'String',['Vorticity' newline '(Upper Layer)'],'EdgeColor','none','FontSize',11);
annotation('textbox',[0.89,0.53,0.15,0.01],'String',['Vorticity' newline '(Lower Layer)'],'EdgeColor','none','FontSize',11);
annotation('textbox',[0.89,0.23,0.05,0.01],'String',"Vertical Velocity",'EdgeColor','none','FontSize',11);

% titles

annotation('textbox',[0.17,0.94,0.3,0.01],'String',"Reduction Factor r=1.0",'EdgeColor','none','FontSize',12);
annotation('textbox',[0.60,0.94,0.3,0.01],'String',"Reduction Factor r=0.01",'EdgeColor','none','FontSize',12);



set(fig,'renderer','Painters');
%set(fig,'renderer','zbuffer');

% saveas(gcf,'qg_r1_vs_qg_r001_snapshots','epsc');
% 
% 
saveas(gcf,['/nfs/pool002/users/mkohl/run_correct/'...
'test_runs/adv_correct/full_relax_correct/equations_correct_nd/'...
'primitive_equations_final/primitive_equations_constant_strat_r001_eps04_geostrophic/'...
'figures_paper/qg_r1_vs_qg_r001_snapshots_linear_tau_damping_alpha_1p7_magnified_different_colormap_test_test_test'],'epsc2')

% saveas(gcf,['/nfs/pool002/users/mkohl/run_correct/'...
% 'test_runs/adv_correct/full_relax_correct/equations_correct_nd/'...
% 'primitive_equations_final/primitive_equations_constant_strat_r001_eps04_geostrophic/'...
% 'figures_paper/qg_r1_vs_qg_r001_snapshots_linear_tau_damping_alpha_1p7_magnified_different_colormap.pdf'])



close all;
