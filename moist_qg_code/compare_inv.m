% Compare w_dedalus to the w_inversion using the 2D-finite difference
% solver

close all; clear;

% path = ['/nfs/pool002/users/mkohl/run_correct/test_runs/adv_correct/'...
% 'full_relax_correct/equations_correct_nd/qg_code_three_layer'];


path = [pwd,'/snapshots/snapshots_s1.h5'];

% read in w and RHS
fid = H5F.open(path);
dset_id = H5D.open(fid,'/tasks/w');
space_id = H5D.get_space(dset_id);
[ndims,h5_dims] = H5S.get_simple_extent_dims(space_id);

nx = h5_dims(2); ny = h5_dims(3); nt = h5_dims(1);

% read in variables
t = h5read(path,'/scales/sim_time');

time = 1:5:nt;

% parameters
Lx = 12*pi;
Ly = 12*pi;
Nx = nx;
Ny = ny;
dx = Lx/Nx;
dy = Ly/Ny;
r = 0.01;

% read in w and RHS

for tt = 1:length(time)
tt

start = [1 1 time(tt)];
count = [ny nx 1];

w1_dedalus = h5read(path,'/tasks/w1',start,count);
w2_dedalus = h5read(path,'/tasks/w2',start,count);
RHS1 = h5read(path,'/tasks/rhs1',start,count);
RHS2 = h5read(path,'/tasks/rhs2',start,count);
RHS1 = 2*RHS1';
RHS2 = 2*RHS2';



[w1_inv,w2_inv] = Omega_Solver(RHS1,RHS2,r,Nx,Ny,dx,dy);
w1_inv = w1_inv'; w2_inv = w2_inv';
[w1_inv_dry,w2_inv_dry] = Omega_Solver_dry(RHS1,RHS2,Nx,Ny,dx,dy);
w1_inv_dry = w1_inv_dry';
w2_inv_dry = w2_inv_dry';


skew_ded1(tt) = Skew(w1_dedalus);
skew_ded2(tt) = Skew(w2_dedalus);
skew_inv1(tt) = Skew(w1_inv);
skew_inv2(tt) = Skew(w2_inv);
skew_inv1_dry(tt) = Skew(w1_inv_dry);
skew_inv2_dry(tt) = Skew(w2_inv_dry);
% skew_inv_dry(tt) = Skew(w_inv_dry);
% skew_inv_1dx(tt) = Skew(w_inv_1dx);
% skew_inv_1dy(tt) = Skew(w_inv_1dy);
% skew_inv_1dxdy(tt) = Skew(w_inv_1dxdy);

lambda_ded1(tt) = Lambda(w1_dedalus);
lambda_ded2(tt) = Lambda(w2_dedalus);
lambda_inv1(tt) = Lambda(w1_inv);
lambda_inv2(tt) = Lambda(w2_inv);
lambda_inv1_dry(tt) = Lambda(w1_inv_dry);
lambda_inv2_dry(tt) = Lambda(w2_inv_dry);




end

% % test residual of omega solver
% rr = ones(size(w_inv)); rr(w_inv>0) = r;
% res = Lap(N,dx)*(rr(:).*w_inv(:))-w_inv(:)-RHS(:);
% res = reshape(res,N,N);

% 2-D comparison w_dedalus, and w_inversion
figure('Renderer','painters','Position',[10 10 900 600])
subplot(1,2,1)
imshow(w2_dedalus); colorbar
caxis([-max(max(w2_inv)) max(max(w2_inv))])
colormap(redblue)
set(gca,'Ydir','normal')
title('w-Dedalus')

subplot(1,2,2)
imshow(w2_inv); colorbar
caxis([-max(max(w2_inv)) max(max(w2_inv))])
colormap(redblue)
set(gca,'Ydir','normal')
title('w-Inversion')

% 1D slice plot w_dedalus vs. w_inversion
figure(2)
plot(w2_dedalus(end/2,:)); hold on;
plot(w2_inv(end/2,:)); hold on;
%plot(w_inv_dry(end/2,:))
legend('dedalus','inversion'); legend boxoff
ylabel('w')

% compare the residuals 

% figure(3)
% plot(w_res(end/2,:)); hold on;
% plot(res(end/2,:));
% legend('dedalus','inversion'); legend boxoff
% ylabel('res')

% Compare the skewness
w2 = h5read(path,'/tasks/w2');

figure(4)
plot(t(time),skew_ded1); hold on;
plot(t(time),skew_ded2); hold on;
plot(t(time),skew_inv1);
plot(t(time),skew_inv2); hold on;
plot(t(time),skew_inv1_dry); hold on;
plot(t(time),skew_inv2_dry);
legend('w1','w2','w1inv','w2inv','w1dry','w2dry')
legend boxoff

% plot(t(time),skew_inv_1dx); hold on;
% plot(t(time),skew_inv_1dy); hold on;
% plot(t(time),skew_inv_1dxdy);
%legend('dedalus','inversion','1dx','1dy','1dxdy'); legend boxoff
xlabel('t')
ylabel('Skew w')
%title([num2str(nx),'x',num2str(nx)])
xlim([t(time(1)) t(time(end))])

figure(40)

plot(t(time),lambda_ded1); hold on;
plot(t(time),lambda_ded2); hold on;
plot(t(time),lambda_inv1); hold on;
plot(t(time),lambda_inv2); hold on;
plot(t(time),lambda_inv1_dry); hold on;
plot(t(time),lambda_inv2_dry);
legend('w1','w2','w1in','w2inv','w1dry','w2dry');
legend boxoff
% plot(t(time),lambda_inv_1dx); hold on;
% plot(t(time),lambda_inv_1dy); hold on;
% plot(t(time),lambda_inv_1dxdx);
%legend('dedalus','inversion','1dx','1dy','1dxdy','Location','SouthEast'); legend boxoff
ylabel('\lambda')
%title([num2str(nx),'x',num2str(nx)])
xlim([t(time(1)) t(time(end))])
xlabel('t')


% Compare the rms between w_dedalus and w_inversion
figure(5)
plot(t(time),rms_diff)
title('rms(w_{ded}-w_{inv})')
xlim([t(time(1)) t(time(end))])

% Compare dry vs moist contribution to skewness

figure(6)
plot(t(time),skew_ded); hold on;
plot(t(time),skew_inv); hold on;
plot(t(time),skew_inv_dry)
legend('w','r+RHS','RHS'); legend boxoff
ylabel('Skew w')
title([num2str(nx),'x',num2str(nx)])
xlim([t(time(1)) t(time(end))])

figure(7)
plot(t(time),lambda_ded); hold on;
plot(t(time),lambda_inv); hold on;
plot(t(time),lambda_inv_dry)
legend('dedalus','r+RHS','RHS'); legend boxoff
ylabel('Skew w')
title([num2str(nx),'x',num2str(nx)])
xlim([t(time(1)) t(time(end))])















