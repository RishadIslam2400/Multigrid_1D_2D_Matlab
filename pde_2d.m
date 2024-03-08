% Number of grid points in x and y directions
nx = 4;
ny = 4;

% Grid spacing
dx = 1 / nx;
dy = 1 / ny;

% Diagonal block
diag_block = eye(ny-1) * (- 2/dx^2 - 2/dy^2);
diag_block = diag_block + diag(ones(ny-2, 1) / dy^2, 1);
diag_block = diag_block + diag(ones(ny-2, 1) / dy^2, - 1);

% Creating matrix with all the diagonal block
Matrix = kron(eye(nx-1), diag_block);
Matrix = Matrix + diag(ones((nx-2) * (ny-1), 1) / dx^2, ny-1);
Matrix = Matrix + diag(ones((nx-2) * (ny-1), 1) / dx^2, -(ny-1));

% Create RHS
x = (1:nx-1) * dx;
y = (1:ny-1) * dy;
[X, Y] = meshgrid(x, y);
F = sin(X) .* cos(Y);
f = reshape(F, (nx-1)*(ny-1), 1);
figure(1);
surf(X, Y, F);
xlabel('X');
ylabel('Y');

% Solve the system
u = Matrix \ f;
U = reshape(u, ny-1, nx-1);
figure(2);
surf(X, Y, U);
xlabel('X');
ylabel('Y');

% Non zero boundary condition on x = 0
% Boundary condition
ux0 = ones(ny-1, 1);

% RHS
f(1:ny-1) = f(1:ny-1) - ux0 / dx^2;

% Solve the system with new boundary condition
u = Matrix \ f;
U = reshape(u, ny-1, nx-1);
figure(3);
surf(X, Y, U);
xlabel('X');
ylabel('Y');
