function [xh,rh,r_nrm,e_nrm,nit] = mg(Ah,fh,x0h,tol,maxit, v1)

% number of relaxations
nrel = 1;

% initialize
nh = length(x0h);
xh = zeros(nh,1);
rh = zeros(nh,1);
r_nrm = zeros(maxit+1,1);
e_nrm = zeros(maxit+1,1);

xh = x0h;
rh = fh-Ah*xh;

r_nrm(1) = norm(rh,2)*(nh+1)^(-1/2); % discrete L2-norm
e_nrm(1) = norm(xh,2)*(nh+1)^(-1/2); % discrete L2-norm

fprintf('mg: ||e0|| = %9.5e \n', e_nrm(1));

nit = -1;
for k = 1:maxit
  if r_nrm(k)<tol 
	nit = k-1;
	break
  end
  fprintf('mg: V-cycle %4d ||e|| (before presmoothing) = %9.5e \n', k, norm(xh,2)*(nh+1)^(-1/2));
  [xh,rh] = red_black_gauss_seidel(Ah,fh,xh,tol,nrel); % presmoothing
  fprintf('mg: V-cycle %4d ||e|| (after presmoothing) = %9.5e \n', k, norm(xh,2)*(nh+1)^(-1/2));

  [eh] = cg_corr(Ah,rh,tol,nrel, v1);
  xh = xh+eh;

  fprintf('mg: V-cycle %4d ||e|| (before postsmoothing) = %9.5e \n', k, norm(xh,2)*(nh+1)^(-1/2));
  [xh,rh] = red_black_gauss_seidel(Ah,fh,xh,tol,nrel); % postsmoothing

  r_nrm(k+1) = norm(rh,2)*(nh+1)^(-1/2); % discrete L2-norm
  e_nrm(k+1) = norm(xh,2)*(nh+1)^(-1/2); % discrete L2-norm - assuming 0 solution 
  
  fprintf('mg: V-cycle %4d ||e||  (after postsmoothing) = %9.5e \n', k, e_nrm(k+1));
end

if nit == -1
  nit = -maxit;
else
  r_nrm = r_nrm(1:nit+1);
  e_nrm = e_nrm(1:nit+1);
end

% disp('end mg');