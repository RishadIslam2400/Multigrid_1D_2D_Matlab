% Solve by multigrid
% equation -u_xx + bu_x + cu = f
%
global b;
global c;

% Sets coefficient such that we solve Poisson
b = 0;
c = 0;

% Number of grid points, and subsequent step size
ne = 127;
he = 1/(ne+1); hes = he^2;
de = ones(ne,1);

% straightforward discretization 
Ae = spdiags([-(1/hes+0.5*b/he)*de, (2/hes+c)*de, -(1/hes-0.5*b/he)*de],-1:1,ne,ne);

% rhs
fe = zeros(ne,1);

% random initial guess
% rng(314159)
% x0 = 100*rand(ne,1);

% fourier mode initial guess
x0 = zeros(ne, 1);
for j = 1:length(x0)
    x0(j) = (sin(3 * j * pi / ne) + sin(10 * j * pi / ne)) / 2;
end

figure(1)
plot([0 x0']','-r','LineWidth',6);
xlim([1, 128]);
ylim([-1, 1]);
% fprintf('hit spacebar to continue \n\n')
% pause

% relaxation parameter
v1 = 2/3;

maxit = 100;
tol = 1.e-10;
[x1,r1,r_nrm1,e_nrm1,nit1] = mg(Ae,fe,x0,tol,maxit, v1);


figure(2)
plot([0 x1']','-r','LineWidth',6);
xlim([1, 128]);
ylim([-1, 1]);