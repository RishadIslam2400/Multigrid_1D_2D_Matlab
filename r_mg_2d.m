% Solve by multigrid
% equation -u_xx - u_yy = 2[(1-6x^2)y^2(1-y^2) + (1-6y^2)x^2(1-x^2)]
%

% Number of grid points, and subsequent step size
nx = 32; % x direction
ny = 32; % y direction
hx = 1/nx; hx_squared = hx^2;
hy = 1/ny; hy_squared = hy^2;

% Create differential operator with straghtforward finitie difference
% discretizatrion
diagonal_block = eye(ny-1) * (2/hx_squared + 2/hy_squared);
diagonal_block = diagonal_block + diag(-ones(ny-2, 1) / hy_squared, 1);
diagonal_block = diagonal_block + diag(-ones(ny-2, 1) / hy_squared, -1);

Ae = kron(eye(nx-1), diagonal_block);
Ae = Ae + diag(-ones((nx - 2) * (ny - 1), 1) / hx_squared, ny-1);
Ae = Ae + diag(-ones((nx - 2) * (ny - 1), 1) / hx_squared, -(ny-1));

Ae = sparse(Ae);

% Create the rhs meshgrid
x = (1:nx-1) * hx;
y = (1:ny-1) * hy;
[X, Y] = meshgrid(x, y);
F = 2 .* ((1 - 6 .* X.^2) .* Y.^2 .* (1-Y.^2) + (1-6 .* Y.^2) .* X.^2 .* (1 - X.^2));
% Surface plot of the RHS
figure(1);
surf(X, Y, F);
xlabel('X');
ylabel('Y');

% Create the rhs vector
fe = reshape(F, (nx-1)*(ny-1), 1);

% random initial guess
rng(314159)
x0 = 100*rand((nx-1)*(ny-1),1);

% fourier mode initial guess
% x0 = zeros(ne, 1);
% for j = 1:length(x0)
%     x0(j) = (sin(3 * j * pi / ne) + sin(10 * j * pi / ne)) / 2;
% end

figure(2)
plot(x0,'-r','LineWidth',6);
% fprintf('hit spacebar to continue \n\n')
% pause

% Create the exact solution vector
U = (X.^2 - X.^4) .* (Y.^4 - Y.^2);
u = reshape(U, (nx-1)*(ny-1), 1);

% Surface plot of the solution
figure(3);
surf(X, Y, U);
xlabel('X');
ylabel('Y');

% Plot the solution vector
figure(4)
plot(u,'-r','LineWidth',6);

% Compute and plot the initial error
e = u - x0;
figure(5)
plot(e,'-r','LineWidth',6);

% relaxation parameter
% v1 = 2/3;
% 
% maxit = 100;
% tol = 1.e-10;
% [x1,r1,r_nrm1,e_nrm1,nit1] = mg(Ae,fe,x0,tol,maxit, v1);
% 
% 
% figure(2)
% plot([0 x1']','-r','LineWidth',6);
% xlim([1, 128]);
% ylim([-1, 1]);