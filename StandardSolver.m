function u = StandardSolver(beta, f, bc0, bc1, bc2, bc3)
% solves Poisson equation with backslash operator

N = length(f);
dx = 1/(N+1);
B = bis2D(N);
Tdxx = B/dx^2 ;
f(1,:)=f(1,:) - bc0/dx^2;               % all x, y = 0 
f(N,:)=f(N,:) - bc1/dx^2;               % all x, y = 1
f(:,1)=f(:,1) - bc2/dx^2;               % all y, x = 0 
f(:,N)=f(:,N) - bc3/dx^2;               % all y, x = 1

f = f(:);

u = -(Tdxx+eye(N^2)*beta)\f;

u = vec2mat(u, N);

end