% plot vorticity and vertical velocity for high and low rossby number

close all; clear;

% eps=0.4

path = [pwd,'/snapshots/snapshots_s1.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nz = h5_dims(4); nt = h5_dims(1);

x = linspace(0,6*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

fig = figure('Position', [10 10 1000 1200]);
%hfig = figure('Position', [10 10 1200 600]);

[ha,pos] = tight_subplot(3,2,[0.03 0.03],[.035 .03],[.035 .1]);

L = 6*pi;
%cmlev = 31;
cmlev = 31;
dispy = 0.01;


for tt = nt:5:nt
    
start_low = [3 1 1 tt];
count_low = [1 nx nx 1];

start_mid = [5 1 1 tt];
count_mid = [1 nx nx 1];

start_up = [7 1 1 tt];
count_up = [1 nx nx 1];
 
zeta1= h5read(path,'/tasks/zeta',start_up,count_up); zeta1 = squeeze(zeta1);
zeta2= h5read(path,'/tasks/zeta',start_low,count_low); zeta2 = squeeze(zeta2);
w = h5read(path,'/tasks/w',start_mid,count_mid); w = squeeze(w);

%subplot(3,2,1)
axes(ha(1))

%cmax = 20;
cmax = 8;
nint = 101;
cint = cmax/nint;

contourf(X,Y,zeta1,-cmax:cint:cmax,'Edgecolor','none');
%imshow(zeta1,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128));
caxis([-cmax cmax])
%colormap(redbluecmap(nint));
colorbar
cm = redbluecmap;
cmlev = nint;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)

%xlabel('x'); 
ylabel('y');
%axis on
%xticks([0 3*pi 6*pi])
%xticklabels({'0','3\pi','6\pi'})
set(gca,'XTicklabel',[]);
yticks([0 3*pi 6*pi])
yticklabels({'0','3\pi','6\pi'})
set(gca,'Ydir','normal')


%subplot(3,2,3)
axes(ha(3))

cmax = 8;
nint = 101;
cint = cmax/nint;

contourf(X,Y,zeta2,-cmax:cint:cmax,'Edgecolor','none');
%imshow(zeta2,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128));
caxis([-cmax cmax])
%colormap(redbluecmap(nint));
colorbar
cm = redbluecmap;
cmlev = nint;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)

%xlabel('x'); 
ylabel('y');
%axis on
%xticks([0 3*pi 6*pi])
%xticklabels({'0','3\pi','6\pi'})
set(gca,'XTicklabel',[]);
yticks([0 3*pi 6*pi])
yticklabels({'0','3\pi','6\pi'})
set(gca,'Ydir','normal')


%subplot(3,2,5)
axes(ha(5))

%cmax = max(abs(w(:)));
cmax = 3;
nint = 101;
cint = cmax/nint;

contourf(X,Y,w,-cmax:cint:cmax,'Edgecolor','none');
%imshow(w,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128))
caxis([-cmax cmax])
%colormap(redbluecmap(nint));
cm = redbluecmap;
cmlev = nint;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)

cbr = colorbar;
set(cbr,'YTick',-cmax:1:cmax);

xlabel('x'); 
ylabel('y');
%axis on
xticks([0 3*pi 6*pi])
xticklabels({'0','3\pi','6\pi'})
yticks([0 3*pi 6*pi])
yticklabels({'0','3\pi','6\pi'})
set(gca,'Ydir','normal')


end

% eps=0.01

path_root = '/nobackup1c/users/mkohl/linear_gradient_formulation/primitive_equations_constant_strat_r001_eps001_final_high_temp_output';
path = [path_root,'/snapshots/snapshots_s1.h5'];

% find the dimension of variables

fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/zeta');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); nz = h5_dims(4); nt = h5_dims(1);

x = linspace(0,6*pi,nx); [X,Y] = meshgrid(x,x);
t = h5read(path,'/scales/sim_time');

L = 6*pi;
%cmlev = 31;
cmlev = 31;
dispy = 0.01;


for tt = nt-15:5:nt-15
    
start_low = [3 1 1 tt];
count_low = [1 nx nx 1];

start_mid = [5 1 1 tt];
count_mid = [1 nx nx 1];

start_up = [7 1 1 tt];
count_up = [1 nx nx 1];
 
zeta1= h5read(path,'/tasks/zeta',start_up,count_up); zeta1 = squeeze(zeta1);
zeta2= h5read(path,'/tasks/zeta',start_low,count_low); zeta2 = squeeze(zeta2);
w = h5read(path,'/tasks/w',start_mid,count_mid); w = squeeze(w);

%subplot(3,2,2)
axes(ha(2))

cmax = 8;
nint = 101;
cint = cmax/nint;

contourf(X,Y,zeta1,-cmax:cint:cmax,'Edgecolor','none');
%imshow(zeta1,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128));
caxis([-cmax cmax])
%colormap(redbluecmap(nint));
colorbar
cm = redbluecmap;
cmlev = nint;
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

cmax = 8;
nint = 101;
cint = cmax/nint;

contourf(X,Y,zeta2,-cmax:cint:cmax,'Edgecolor','none');
%imshow(zeta2,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128));
caxis([-cmax cmax])
%colormap(redbluecmap(nint));
colorbar
cm = redbluecmap;
cmlev = nint;
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
cmax = 10;
nint = 101;
cint = cmax/nint;

contourf(X,Y,w,-cmax:cint:cmax,'Edgecolor','none');
%imshow(w,'XData',linspace(0,6*pi,128),'Ydata',linspace(0,6*pi,128))
caxis([-cmax cmax])
%colormap(redbluecmap(nint));
cm = redbluecmap;
cmlevel = nint;
cm = imresize(cm,[cmlev,3]); cm = min(max(cm,0),1);
cm(0.5*(cmlev+1),:) = [1 1 1];
colormap(cm)

cbr = colorbar;
set(cbr,'YTick',-cmax:5:cmax);
xlabel('x'); 
%ylabel('y')

axis on
set(gca,'Ydir','normal')
xticks([0 3*pi 6*pi])
xticklabels({'0','3\pi','6\pi'})
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
   'primitive_equations_constant_strat_r001_eps04_geostrophic/figures_drv_world/test2'];
 
saveas(gca,path_save,'epsc');