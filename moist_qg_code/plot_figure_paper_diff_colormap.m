% plot vorticity and vertical velocity for high and low rossby number

close all; clear;

% r = 1

path = ['/nobackup1c/users/mkohl/qg_r1/snapshots/snapshots_s1.h5'];

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


fig = figure('Position', [10 10 1000 1200]);
%hfig = figure('Position', [10 10 1200 600]);

[ha,pos] = tight_subplot(3,2,[0.03 0.03],[.035 .03],[.035 .1]);

L = 12*pi;
%cmlev = 31;
cmlev = 31;
dispy = 0.01;


for tt = nt:5:nt
    
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


%subplot(3,2,1)
axes(ha(1))

cmax = 20;
nint = 101;
cint = cmax/nint;

contourf(X,Y,zeta1,-cmax:cint:cmax,'Edgecolor','none');
%imshow(zeta1,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128));
caxis([-cmax cmax])
%colormap(redblue(nint));
colorbar
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)

%xlabel('x'); 
ylabel('y');
%axis on
%xticks([0 3*pi 6*pi])
%xticklabels({'0','3\pi','6\pi'})
set(gca,'XTicklabel',[]);
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
set(gca,'Ydir','normal')


%subplot(3,2,3)
axes(ha(3))

cmax = 6;
nint = 101;
cint = cmax/nint;

contourf(X,Y,zeta2,-cmax:cint:cmax,'Edgecolor','none');
%imshow(zeta2,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128));
caxis([-cmax cmax])
%colormap(redblue(nint));
colorbar
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)

%xlabel('x'); 
ylabel('y');
%axis on
%xticks([0 3*pi 6*pi])
%xticklabels({'0','3\pi','6\pi'})
set(gca,'XTicklabel',[]);
yticks([0 6*pi 12*pi])
yticklabels({'0','6\pi','12\pi'})
set(gca,'Ydir','normal')


%subplot(3,2,5)
axes(ha(5))

cmax = 13; %max(abs(w(:)));
%cmax = 3;
nint = 101;
cint = cmax/nint;

contourf(X,Y,w,-cmax:cint:cmax,'Edgecolor','none');
%imshow(w,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128))
caxis([-cmax cmax])
%colormap(redblue(nint));
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)

cbr = colorbar;
set(cbr,'YTick',-cmax:5:cmax);

xlabel('x'); 
ylabel('y');
%axis on
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
yticks([0 3*pi 6*pi])
yticklabels({'0','12\pi','12\pi'})
set(gca,'Ydir','normal')


end

% r = 0.01

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

for tt = nt-15:5:nt-15

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

%subplot(3,2,2)
axes(ha(2))

cmax = 15;
nint = 101;
cint = cmax/nint;

contourf(X,Y,zeta1,-cmax:cint:cmax,'Edgecolor','none');
%imshow(zeta1,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128));
caxis([-cmax cmax])
%colormap(redblue(nint));
colorbar
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)

%xlabel('x'); 
%ylabel('y');
%axis on
%xticks([0 3*pi 6*pi])
%xticklabels({'0','3\pi','6\pi'})
%yticks([0 3*pi 6*pi])
%yticklabels({'0','3\pi','6\pi'})
set(gca,'XTicklabel',[]);
set(gca,'YTicklabel',[]);
%set(gca,'Ydir','normal')

%subplot(3,2,4)
axes(ha(4))

cmax = 15;
nint = 101;
cint = cmax/nint;

contourf(X,Y,zeta2,-cmax:cint:cmax,'Edgecolor','none');
%imshow(zeta2,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128));
caxis([-cmax cmax])
%colormap(redblue(nint));
colorbar
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)
%xlabel('x'); ylabel('y')
%axis on
%xticks([0 3*pi 6*pi])
%xticklabels({'0','3\pi','6\pi'})
%yticks([0 3*pi 6*pi])
%yticklabels({'0','3\pi','6\pi'})
set(gca,'XTicklabel',[]);
set(gca,'YTicklabel',[]);
%set(gca,'Ydir','normal')
caxis([-cmax cmax])

%subplot(3,2,6)
axes(ha(6))

%cmax = max(abs(w(:)));
cmax = 100;
nint = 101;
cint = cmax/nint;

contourf(X,Y,w,-cmax:cint:cmax,'Edgecolor','none');
%imshow(w,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128))
caxis([-cmax cmax])
%colormap(redblue(nint));
cm = redbluecmap;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)

cbr = colorbar;
set(cbr,'YTick',-cmax:20:cmax);
xlabel('x'); 
%ylabel('y')

axis on
set(gca,'Ydir','normal')
xticks([0 6*pi 12*pi])
xticklabels({'0','6\pi','12\pi'})
% yticks([0 3*pi 6*pi])
% yticklabels({'0','3\pi','6\pi'})

set(gca,'YTicklabel',[]);

end

% a,b,c,d

annotation('textbox',[0.001, 0.965, 0.01,0.03],'String',"(a)",'EdgeColor','none');
annotation('textbox',[0.45, 0.965, 0.01,0.03],'String',"(b)",'EdgeColor','none');

annotation('textbox',[0.001, 0.645, 0.01,0.03],'String',"(c)",'EdgeColor','none');
annotation('textbox',[0.45, 0.645, 0.01,0.03],'String',"(d)",'EdgeColor','none');

annotation('textbox',[0.001, 0.32, 0.01,0.03],'String',"(e)",'EdgeColor','none');
annotation('textbox',[0.45, 0.32, 0.01,0.03],'String',"(f)",'EdgeColor','none');


% textboxes
annotation('textbox',[0.885,0.83,0.15,0.01],'String',['Vorticity' newline '(Upper Layer)'],'EdgeColor','none','FontSize',11);
annotation('textbox',[0.885,0.51,0.15,0.01],'String',['Vorticity' newline '(Lower Layer)'],'EdgeColor','none','FontSize',11);
annotation('textbox',[0.885,0.19,0.05,0.01],'String',"Vertical Velocity",'EdgeColor','none','FontSize',11);

% titles

annotation('textbox',[0.12,0.99,0.3,0.01],'String',"High Rossby Regime",'EdgeColor','none','FontSize',12);
annotation('textbox',[0.56,0.99,0.3,0.01],'String',"Low Rossby Regime",'EdgeColor','none','FontSize',12);

%set(gcf,'Units','Normalized','OuterPosition',[0,0,1,1]);
%set(gcf,'Units','pixels');

% save files

%set(fig,'renderer','Painters');

path_save = ['/nfs/pool002/users/mkohl/run_correct/test_runs/'...
'adv_correct/full_relax_correct/equations_correct_nd/primitive_equations_final/'...
   'primitive_equations_constant_strat_r001_eps04_geostrophic/figures_drv_world/test_qg_diff'];
 
saveas(gca,path_save,'epsc');