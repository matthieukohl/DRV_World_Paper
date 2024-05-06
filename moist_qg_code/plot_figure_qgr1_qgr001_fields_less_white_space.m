% Plot the qg_r1 and qg_r001 fields 
close all; clear;

path = ['/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/'...
    'full_relax_correct/equations_correct_nd/qg_final/qg_r1/snapshots/snapshots.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta_bc');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nt = h5_dims(1);

x = linspace(0,12*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

%fig = figure('Position', [10 10 1000 1200]);
%fig = figure('Position', [10 10 1000 1200]);
fig = figure('Position', [10 10 5/4*1000 5/4*1200]);
%[ha,pos] = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01]);
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 2 1])
% 

v = VideoWriter([pwd,'/video/test.mp4','MPEG-4']);
%v.FrameRate = 2;
open(v);

colorint = 100;
cmlev = 21;
dispy = 0;

for tt = nt:5:nt
tt    
 
start = [1 1 tt];
count = [nx nx 1];
 
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

%subplot(3,2,1)
%imshow(zeta1)
imshow(zeta1,'XData',linspace(0,12*pi,512),'Ydata',linspace(0,12*pi,512)); colorbar; caxis([-max(max(zeta1)) max(max(zeta1))])
%contourf(X,Y,zeta1,'EdgeColor','none');
%colormap(redblue(colorint))
%colormap(redbluecmap(colorint))
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)
colorbar
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
%caxis([-10 10])
%caxis([-max(abs(zeta1(:))) max(abs(zeta1(:)))])
caxis([-20 20])
axes(ha(1));
% 
%subplot(3,2,3)
%imshow(zeta2)
imshow(zeta2,'XData',linspace(0,12*pi,512),'Ydata',linspace(0,12*pi,512)); colorbar; caxis([-max(max(zeta2)) max(max(zeta2))])
%contourf(X,Y,zeta2,'EdgeColor','none');
xlabel('x'); ylabel('y')
%colormap(redblue(colorint))
%colormap(redbluecmap(colorint))
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)
colorbar
axis on
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
%t2 = title('\rm Vorticity (Lower Layer)');
%tpos = get(t2,'position');
%set(t2,'position',tpos + [0 dispy 0 ]);
set(gca,'Ydir','normal')
%caxis([-10 10])
%caxis([-max(abs(zeta2(:))) max(abs(zeta2(:)))])
caxis([-6 6])
axes(ha(3));

%subplot(3,2,5)
%imshow(w)
imshow(w,'XData',linspace(0,12*pi,512),'Ydata',linspace(0,12*pi,512)); colorbar; 
%caxis([-max(max(w)) max(max(w))])
%contourf(X,Y,w,'EdgeColor','none');
%colormap(redblue(colorint))
%colormap(redbluecmap(colorint))
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)
cbr = colorbar;
set(cbr,'YTick',-8:4:8);
xlabel('x'); ylabel('y')
%t3 = title('\rm Vertical Velocity');
%tpos = get(t3,'position');
%set(t3,'position',tpos + [0 dispy 0 ]);
axis on
set(gca,'Ydir','normal')
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
caxis([-8 8])
caxis([-max(abs(w(:))) max(abs(w(:)))])
axes(ha(5));

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

for tt = nt:5:nt
tt    
 
start = [1 1 tt];
count = [nx nx 1];
 
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

%subplot(3,2,2)
%imshow(zeta1)
imshow(zeta1,'XData',linspace(0,12*pi,512),'Ydata',linspace(0,12*pi,512)); colorbar; caxis([-max(max(zeta1)) max(max(zeta1))])
%axes(ha(4))
%contourf(X,Y,zeta1,'Edgecolor','none');
%colormap(redblue(colorint))
%colormap(redbluecmap(colorint))
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)
colorbar
xlabel('x'); ylabel('y')
%title('\rm Vorticity (Upper Layer)')
axis on
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
set(gca,'Ydir','normal')
%caxis([-10 10])
caxis([-max(abs(zeta1(:))) max(abs(zeta1(:)))])
caxis([-20 20])
axes(ha(2));
% 
%subplot(3,2,4)
%imshow(zeta2)
imshow(zeta2,'XData',linspace(0,12*pi,512),'Ydata',linspace(0,12*pi,512)); colorbar; caxis([-max(max(zeta2)) max(max(zeta2))])
%contourf(X,Y,zeta2,'Edgecolor','none');
xlabel('x'); ylabel('y')
%colormap(redblue(colorint))
%colormap(redbluecmap(colorint))
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)
colorbar
axis on
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
%title('\rm Vorticity (Lower Layer)')
set(gca,'Ydir','normal')
%caxis([-10 10])
%caxis([-max(zeta2(:)) max(zeta2(:))])
%caxis([-max(abs(zeta2(:))) max(abs(zeta2(:)))])
caxis([-20 20])
axes(ha(4))

%subplot(3,2,6)
%imshow(w)
imshow(w,'XData',linspace(0,12*pi,512),'Ydata',linspace(0,12*pi,512)); colorbar; 
%set(h,'Position',pos{6});
%caxis([-max(max(w)) max(max(w))])
%contourf(X,Y,w,'Edgecolor','none');
%colormap(redblue(colorint))
%colormap(redbluecmap(colorint))
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)
cbr = colorbar;
set(cbr,'YTick',-140:70:140);
%caxis([-140 140])
%caxis([-max(w(:)) max(w(:))])
xlabel('x'); ylabel('y')
%title('\rm Vertical Velocity')
axis on
set(gca,'Ydir','normal')
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
%caxis([-100 100])
caxis([-max(abs(w(:))) max(abs(w(:)))])
axes(ha(6));

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

% saveas(gcf,['/nfs/pool002/users/mkohl/run_correct/'...
% 'test_runs/adv_correct/full_relax_correct/equations_correct_nd/'...
% 'primitive_equations_final/primitive_equations_constant_strat_r001_eps04_geostrophic/'...
% 'figures_paper/qg_r1_vs_qg_r001_snapshots_less_white_space'],'epsc')
