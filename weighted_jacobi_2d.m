function [x,r] = weighted_jacobi_2d(A,b,x0,tol,max_it, weight)

omega = weight;

r_nrm = zeros(max_it+1,1);
D = diag(diag(A));
x = x0;

r = b-A*x;
r_nrm(1) = norm(r);
rt = omega*(D\r);

k=1;
while (r_nrm(k) > tol) && (k <= max_it)
  x = x + rt;
  r = b-A*x;  
  r_nrm(k+1) = norm(r);
  rt = omega*(D\r);
  k = k+1;
end
n_it = k; % Save number of iterations, can return optionally