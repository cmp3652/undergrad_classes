% Clara Pitkins
% cmp3652@rit.edu
% Numerical Linear Algebra
% Inverse Power Iteration (iterative eigenvalue algorithm)
% Allows one to find an approximate eigenvector when an approximation to a
% corresponding eigenvalue is already known.
% 2023

%Initialize Values
A = [2,1,0;1,3,1;0,1,2];
sigma = 0;
xo = [0;-1;2];
n = size(A,1);

y = xo/norm(xo);
max_iter = 20;

%Expected Value
yex = [0.4082;0.8165;0.4082];

% Iteration
for k = 1:max_iter
    x = (A - sigma*eye(n))\y ;
    y = x/norm(x) ;
    u = transpose(y)*A*y ;
    error = norm(yex - y);
    
    % Uncomment if you want to display for every iteration
    % disp(['iteration ', num2str(k),', eigenvalue:', num2str(1/u,8), sprintf(' eigenvector: (%d, %d, %d)', y), ', error:', num2str(error,8)]);
end

disp(['iteration ', num2str(k),', eigenvalue:', num2str(1/u,8), sprintf(' eigenvector: (%d, %d, %d)', y), ', error:', num2str(error,8)]);
