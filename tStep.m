function [unew, vnew] = tStep(uold, vold, dt, bc0, bc1, bc2, bc3, N, k)
% takes one time step dt. Uses full multigrid v-cycle.

N = length(uold);
dx = 1/(N+1);

uConv = [0 1 0; 1 -4 1; 0 1 0]/dx^2*4/dt;
vConv = [0 1 0; 1 -4 + 4/dt^2*dx^2 1; 0 1 0]/dx^2;

f = conv2(uold, uConv, 'same') + conv2(vold,vConv,'same');
if 1
    f_long = reshape(f, [N*N 1]);

    gamma = -4/dt^2;
    tic
    vnew = mxMG_2D(f_long, 1, k, 5, gamma);
    toc
    vnew = reshape(vnew, [N N]);
else
    tic
    vnew = FMGV2D(-4/dt^2, vold, f, bc0, bc1, bc2, bc3);
    toc
end

unew = uold + dt*(vold + vnew)/2;
end