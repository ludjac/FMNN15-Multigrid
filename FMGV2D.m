function [u, ef]= FMGV2D(beta, u, f, bc0, bc1, bc2, bc3)
% full multigrid V-cycle. Solves: -delta*u - beta*u = f. Uses recursion.

omega= 1;
N = length(u);
dx = 1/(N+1);


if N < 32
    u = StandardSolver(beta, f, bc0, bc1, bc2, bc3);   % solves w. backslash
else
    rf = residual2D(beta, u, f, bc0, bc1, bc2, bc3);  % calculates the residual *changes sign on rf
    u = u - omega*(4/dx^2-beta)^(-1)*rf;         % Pre-smoothing Jacobi * 

    rf = residual2D(beta, u, f, bc0, bc1, bc2, bc3); % calculates the residual

    rf = lowpass2D(rf);
    
    rc = FMGrestrict2D(rf);               % creates the coarse grid

    ec = FMGV2D(beta, zeros(size(rc)), rc, 0, 0, 0, 0); % recursive
    
    ef = FMGprolong2D(ec);               % vector back to ordinary size
    u = u - ef;                         % removes the error

    rf = residual2D(beta, u, f, bc0, bc1, bc2, bc3);  % computes residual
    u = u - omega*(4/dx^2-beta)^(-1)*rf;         % post-smoothing Jacobi

end

end