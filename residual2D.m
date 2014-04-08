function res = residual2D(beta, u, f, bc0, bc1, bc2, bc3)
% calculates and returns the residual
N = length(f);
dx = 1/(N+1);
[0 1 0; 1 -4 1; 0 1 0]/dx^2+[0 0 0; 0 beta 0; 0 0 0];
beta;
Tdxu = conv2(u,([0 1 0; 1 -4 1; 0 1 0]/dx^2+[0 0 0; 0 beta 0; 0 0 0]),'same');

f(1,:)=f(1,:) - bc0/dx^2;               % all x, y = 0 
f(N,:)=f(N,:) - bc1/dx^2;               % all x, y = 1
f(:,1)=f(:,1) - bc2/dx^2;               % all y, x = 0 
f(:,N)=f(:,N) - bc3/dx^2;               % all y, x = 1

res = - Tdxu - f;

end