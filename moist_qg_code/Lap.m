function [A] = Lap(N,Ny,dx,dy)
% 2d finite difference operator along x

NN = N*Ny;

% Define the A-matrix 
    
dA = sparse(diag((-2/dx^2-2/dy^2)*ones(1,N))); % diagonal matrix
dAp1x = sparse(diag( 1/dx^2*ones(1,N-1), 1 )); % super-diagonal matrix
dAm1x = sparse(diag( 1/dx^2*ones(1,N-1), -1 )); % sub-diagonal matrix

Asmall = dA + dAp1x + dAm1x;
Asmall(1,end) = 1/dx^2;
Asmall(end,1) = 1/dx^2;

Abig = kron(sparse(eye(Ny)),Asmall);

dAp1y = sparse(diag(1/dy^2*sparse(ones(1,NN-N)), N ));
dAm1y = sparse(diag(1/dy^2*sparse(ones(1,NN-N)),-N));


dAuy = sparse(diag(1/dy^2*sparse(ones(1,NN-N*(Ny-1))),N*(Ny-1))); 
dAly = sparse(diag(1/dy^2*sparse(ones(1,NN-N*(Ny-1))),-N*(Ny-1)));

A = (Abig + dAp1y + dAm1y + dAuy + dAly);


end

