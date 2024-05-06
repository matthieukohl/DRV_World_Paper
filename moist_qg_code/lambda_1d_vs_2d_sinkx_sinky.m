% compare 1d vs 2d inversion with random noise

%close all; clear;

kx = 6.6;
ky = 4.0;

L = 12*pi;
Nx = 512;
Ny = 512;
%Nx = 128;
%Ny = 128;

x = linspace(0,L,Nx);
y = linspace(0,L,Ny);

[X,Y] = meshgrid(x,y);

dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));

draws = 1;

r = [0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0];

[lambda_vs_r,lambda_1d_vs_r] = deal(zeros(length(r),1));

for ii = 1:length(r)
    ii
    
for jj = 1:draws
    
%RHS = randn(Nx,Ny);
%RHS = RHS/max(abs(RHS(:)));
RHS = sin(kx.*X).*sin(ky*Y);


[w] = Omega_Solver(RHS,r(ii),Nx,Ny,dx,dy);

lambda_vs_r(ii) = lambda_vs_r(ii) + Lambda(w);

clear('w')

%RHS_1d = randn(Nx,1);
[w] = Omega_Solver1D(RHS(:),r(ii),Nx*Nx,dx);
%RHS_1d = RHS(1,:)'; 
%[w] = Omega_Solver1D(RHS_1d,r(ii),Nx,dx);

lambda_1d_vs_r(ii) = lambda_1d_vs_r(ii) + Lambda(w);

clear('w')
end
    
lambda_vs_r(ii) = lambda_vs_r(ii)/draws;
lambda_1d_vs_r(ii) = lambda_1d_vs_r(ii)/draws;

end

% lambda_1d_theory = 1./(1+r);
%lambda_2d_theory = 1./(1+r.^2);
lambda_pablo = 1./(1+sqrt(r));

%figure
semilogx(r,lambda_1d_vs_r,'r-o'); hold on;
semilogx(r,lambda_vs_r,'b-o'); hold on;
semilogx(r,lambda_pablo,'k-o');
legend('1d','2d','Pablo'); legend boxoff;
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
set(gca,'FontSize',12);
hold on;
