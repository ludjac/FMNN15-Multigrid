function B = bis2D(N)
% returns a NxN matrix with diagonal values for discrete approximation of
% second degree derivative for a 2D problem
N2 = N^2;
sub = ones((N2-1),1);
sub(N:N:end) = 0;
sup = sub;
main = -4*ones((N2),1);
subb = ones((N2-N),1);
supp = subb;

B = diag(subb,-N) + diag(sub,-1) + diag(main,0) + diag(sup,1) + diag(supp,N);

%B = B*0;
end