function [w1,w2] = Omega_Solver(RHS1,RHS2,r,N,Ny,dx,dy)

% Solve 2-D moist omega equation in both layers
% Lap(r(w1)w1)-(2*w1-w2) = 2*RHS1
% Lap(r(w2)w2)-(2*w2-w1) = 2*RHS2
% where r(w1) = r1 = 1, and r(w2) = r2 = r

%tic

RHS1 = RHS1(:); NN = length(RHS1);
RHS2 = RHS2(:);
RHS = [RHS1;RHS2];

w_old = randn(length(RHS),1); % start with random guess

% Start the iterative-Inversion

E = 10; % Initialize Error-term
counter = 0;
while E>=1e-12
counter = counter + 1;
% Define the r-factor based on w1
    
rr = ones(size(w_old(1:NN))); %rr(w_old(1:NN)>0) = r; 
 
% Define the A1-matrix 
   
vp1x = sparse(rr(2:end)'); vp1x(N:N:end) = 0;
vm1x = sparse(rr(1:end-1)'); vm1x(N:N:end) = 0;

vl1x = sparse(zeros(1,NN-(N-1))); vl1x(1:N:end) = rr(N:N:end);
vr1x = sparse(zeros(1,NN-(N-1))); vr1x(1:N:end) = rr(1:N:end);

dA = sparse(diag(sparse(-2*rr'/dx^2-2*rr'/dy^2-2).*sparse(ones(1,NN))));
dAp1x = sparse(diag( sparse(1/dx^2*vp1x), 1 )); 
dAm1x = sparse(diag( sparse(1/dx^2.*vm1x), -1 ));
dAl1x = sparse(diag(1/dx^2*vl1x,(N-1)));
dAr1x = sparse(diag(1/dx^2*vr1x,-(N-1)));

dAp1y = sparse(diag(1/dy^2*sparse(rr(N+1:1:end)'), N ));
dAm1y = sparse(diag(1/dy^2*sparse(rr(1:1:NN-N)'),-N));


dAuy = sparse(diag(1/dy^2*sparse(rr(NN-(N-1):1:end)'),N*(Ny-1))); 
dAly = sparse(diag(1/dy^2*sparse(rr(1:1:N)'),-N*(Ny-1)));

A1 = (dA + dAp1x + dAm1x + dAl1x + dAr1x + dAp1y + dAm1y + dAuy + dAly);

% Define the r-factor based on w2
    
rr = ones(size(w_old(NN+1:end))); rr(w_old(NN+1:end)>0) = r; 
 
% Define the A2-matrix 
   
vp1x = sparse(rr(2:end)'); vp1x(N:N:end) = 0;
vm1x = sparse(rr(1:end-1)'); vm1x(N:N:end) = 0;

vl1x = sparse(zeros(1,NN-(N-1))); vl1x(1:N:end) = rr(N:N:end);
vr1x = sparse(zeros(1,NN-(N-1))); vr1x(1:N:end) = rr(1:N:end);

dA = sparse(diag(sparse(-2*rr'/dx^2-2*rr'/dy^2-2).*sparse(ones(1,NN))));
dAp1x = sparse(diag( sparse(1/dx^2*vp1x), 1 )); 
dAm1x = sparse(diag( sparse(1/dx^2.*vm1x), -1 ));
dAl1x = sparse(diag(1/dx^2*vl1x,(N-1)));
dAr1x = sparse(diag(1/dx^2*vr1x,-(N-1)));

dAp1y = sparse(diag(1/dy^2*sparse(rr(N+1:1:end)'), N ));
dAm1y = sparse(diag(1/dy^2*sparse(rr(1:1:NN-N)'),-N));


dAuy = sparse(diag(1/dy^2*sparse(rr(NN-(N-1):1:end)'),N*(Ny-1))); 
dAly = sparse(diag(1/dy^2*sparse(rr(1:1:N)'),-N*(Ny-1)));

A2 = (dA + dAp1x + dAm1x + dAl1x + dAr1x + dAp1y + dAm1y + dAuy + dAly);



% take into account that there are two layers top and bottom
Abig = blkdiag(A1,A2);

Aw1 = sparse(diag(1*sparse(ones(1,NN)),NN)); 
Aw2 = sparse(diag(1*sparse(ones(1,NN)),-NN));

Abig = Abig+Aw1+Aw2;



% Carry out the inversion to find the new w

w_new = Abig\RHS;


%E = rms(w_new-w_old);
%E = max(w_new-w_old);

w_old = w_new;

%rr = ones(size(w_old)); rr(w_old(NN+1:end)>0) = r; 

r1 = ones(size(w_old(1:NN))); 
r2 = ones(size(w_old(NN+1:end))); r2(w_old(NN+1:end)>0) = r;

rr = [r1;r2];

E1 = rms(Lap(N,Ny,dx,dy)*(rr(1:NN).*w_old(1:NN))-(2*w_old(1:NN)-w_old(NN+1:end))-RHS(1:NN));
E2 = rms(Lap(N,Ny,dx,dy)*(rr(NN+1:end).*w_old(NN+1:end))-(2*w_old(NN+1:end)-w_old(1:NN))-RHS(NN+1:end));
E = 0.5*(E1+E2);
% wtest = reshape(w_old,N,N); 
% plot(wtest(end/2,:)); hold on;

end

%toc
counter
w1 = reshape(w_old(1:NN),N,Ny);
w2 = reshape(w_old(NN+1:end),N,Ny);
end


