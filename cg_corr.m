function [eh] = cg_corr(Ah,rh,tol,nrel, v1)

global b
global c;

nh = length(rh);
n2h = ceil((nh-1)/2);
x2h = zeros(n2h,1);
tol_int = 1.e-15;

% restriction - injection operator
% f2h = rh(2:2:nh-1);

% restriction - full weighting operator
f2h = (rh(1:2:nh-1) + 2 * rh(2:2:nh-1) + rh(3:2:nh)) / 4;

h = 1/(n2h+1); hs = h^2;
d = ones(n2h,1);

% straightforward discretization
A2h = spdiags([-(1/hs+0.5*b/h)*d, (2/hs+c)*d, -(1/hs-0.5*b/h)*d],-1:1,n2h,n2h); % discretize for larger h

% figure(2)
% plot(A2h\f2h,'b');

% 2-GRID ALGORITHM
% x2h = A2h\f2h;   % exact solution for 2-grid method
% r2h = f2h-A2h*x2h;

% V-CYCLE
X = 5;
if n2h > X % MG if more than X grid points
    [x2h,r2h] = red_black_gauss_seidel(A2h,f2h,x2h,tol_int,nrel); % pre-smoothing
    [e2h] = cg_corr(A2h,r2h,tol,nrel, v1);           % compute coarse grid correction
    x2h = x2h+e2h;                               %   correct approximate solution
    [x2h,r2h] = red_black_gauss_seidel(A2h,f2h,x2h,tol_int,nrel); % post-smoothing
    r2h = f2h-A2h*x2h;
else % otherwise direct solve
    x2h = A2h\f2h;
    r2h = f2h - A2h*x2h;
end

% prolongation - linear interpolation
eh = zeros(nh,1);
eh(2:2:nh-1) = x2h;                          
eh(1:2:nh-2) = 0.5*eh(2:2:nh-1);             
eh(3:2:nh) = eh(3:2:nh) + 0.5*eh(2:2:nh-1);  