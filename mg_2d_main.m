close all;

MGiterations  = 1;

k = 9;
N = 2^k-1;
dx = 1/(N+1);
dt = 0.05;
gamma = -4/dt^2;

x = linspace(0, 1, N+2);
y = linspace(0, 1, N+2);
[X, Y] = meshgrid(x,y);

% setting up u0 and v0
u0 = sin(2*3.14*x(2:end-1))'*sin(2*3.14*y(2:end-1));
v0 = zeros(size(u0));

x0 = 0.3; 
y0 = 0.4; 

A = 2; 
sigx = 0.005; 
sigy = 0.005;

u0 = A*exp(-( (X-x0).^2/(2*sigx) + (Y-y0).^2/(2*sigy)));
%u0 = A*exp(-( (X-x0).^2/(0.1*sigx) + (Y-y0).^2/(0.1*sigy)));
u0 = u0(2:end-1,2:end-1);

u_kernel = [0 1 0; 1 -4 1; 0 1 0]/dx^2*4/dt;
v_kernel = [0 1 0; 1 -4 + 4/dt^2*dx^2 1; 0 1 0]/dx^2;
ufinal = zeros(N+2);
% boundary conditions
bc0 = 0; bc1 = 0; bc2 = 0; bc3 = 0;

%Create a new figure window
figure('Name', 'Waves', 'Position',[0 0 800 600]);

%Here we draw the initial conditions
ufinal(2:end-1,2:end-1) = u0;
h = surf(X,Y,ufinal);

%We need this to control the color-coding.
set(h,'CDataMapping','direct');
uold = u0;
vold = v0;

% the time stepping loop
cool = A*exp(-( (X(2:end-1,2:end-1)-x0).^2/(2*sigx) + (Y(2:end-1,2:end-1)-y0).^2/(2*sigy)));

plotdata = cell(1, 7);
i=1;
for main = 0:80
    % tic
    %TODO: Solve for the next time step using your multigrid
    % solver. Store the next time step in the matrix u.
   
    f = conv2(uold, u_kernel, 'same') + conv2(vold,v_kernel,'same');
    f_long = reshape(f, [N*N 1]);
    
    %tic
    % Invokes the C-program
    vnew = mxMG_2D(f_long, MGiterations, k, 5, gamma);
    %toc
    vnew = reshape(vnew, [N N]);
    unew = uold + dt*(vold + vnew)/2;
    
    %[unew, vnew] = tStep(uold, vold, dt, bc0, bc1, bc2, bc3, N, k);
    % the Helmholtz equation
    
   if main == 100
        unew = unew + cool;
   end
   if main == 200
       unew = unew + cool;
   end
   if main == 300
       unew = unew + cool;
   end
   if main == 400
       unew = unew + cool;
   end
   if main == 500
       unew = unew + cool;
   end
   if main == 600
       unew = unew + cool;
   end
   if main == 700
       unew = unew + cool;
   end
    
    ufinal(2:end-1,2:end-1) = unew;
    
    h = surf(X,Y,ufinal, 'facecolor','interp', 'edgecolor','none');
    
    %Controls the colors of the drawn mesh
    set(h,'CData',32*ufinal+32)
    
    %Change the height data to be that of the new time step
    set(h,'ZData',ufinal);
    if mod(main,10)==0
        plotdata{i}=ufinal;
        i=i+1;
    end
    %The axis normally wants to follow the data. Force it not to.
    axis([0 1 0 1 -1 1]);
    
    %Draw mesh
    drawnow
    uold = unew;
    vold = vnew;
    %pause;
end
