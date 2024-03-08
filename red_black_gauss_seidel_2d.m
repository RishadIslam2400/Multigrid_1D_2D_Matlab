function [x, r] = red_black_gauss_seidel_2d(A,b,x0,tol,max_it)
% Gauss-Seidel relxation using lexicographical ordering
%
% Input:    A: Coefficiet matrix of the linear system
%           b: RHS vector of the linear system
%           x0: Initial guess
%           tol: Convergence tolerance
%           max_it: Maximum number of iterations
%
% Output:   x: Approximate solution
%           r: Residual

% Preprocessing
r_nrm = zeros(max_it+1, 1);
e_nrm = zeros(max_it+1, 1);
x = x0;
k = 1;
n = length(b);
x_new = x;

% Computing norms
r = b - A*x;
r_nrm(1) = norm(r, 2);
e_nrm(1) = norm(x, 2);

% Iterate until max_it
while (k <= max_it) && (r_nrm(k) > tol)
    % Gauss-Seidel iteration formula
    % red points
    for i = 2:2:n
        x_new(i) = (b(i) - (A(i, 1:i-1) * x_new(1:i-1)) - (A(i, i+1:n) * x(i+1:n))) / A(i,i); 
    end
    % black points
    for i = 1:2:n
        x_new(i) = (b(i) - (A(i, 1:i-1) * x_new(1:i-1)) - (A(i, i+1:n) * x(i+1:n))) / A(i,i); 
    end

    x = x_new;
    r = b - A*x;
    r_nrm(k+1) = norm(r, 2);
    e_nrm(k+1) = norm(x, 2);
%     if e_nrm(k+1) < tol
%         break;
%     end    
    k = k + 1;
end

n_it = k; % Save number of iterations, can return optionally
end

