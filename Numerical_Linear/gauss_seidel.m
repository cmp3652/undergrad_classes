% Clara Pitkins
% cmp3652@rit.edu
% Numerical Linear Algebra
% Gauss Seidel Method (Iterative Method for Solving Linear Systems)
% 2023

% Initialize arrays and values
A = [-3,-1,0,0,0,0,0,0,0,0,0,1/2;-1,3,-1,0,0,0,0,0,0,0,1/2,0;0,-1,3,-1,0,0,0,0,0,1/2,0,0;0,0,-1,3,-1,0,0,0,1/2,0,0,0;0,0,0,-1,3,-1,0,1/2,0,0,0,0;
    0,0,0,0,-1,3,-1,0,0,0,0,0;0,0,0,0,0,-1,3,-1,0,0,0,0;0,0,0,0,1/2,0,-1,3,-1,0,0,0;0,0,0,1/2,0,0,0,-1,3,-1,0,0;0,0,1/2,0,0,0,0,0,-1,3,-1,0;
    0,1/2,0,0,0,0,0,0,0,-1,3,-1;1/2,0,0,0,0,0,0,0,0,0,-1,3];
b = [2.5;1.5;1.5;1.5;1.5;1;1;1.5;1.5;1.5;1.5;2.5];
xo = zeros(12,1);
max_iter = 10;

% Call in-line function
xnew = gauss_seidel_method(A,b,xo,max_iter)

% In-line function
function xnew = gauss_seidel_method(A, b, xo, max_iter)
n = size(A,1);
xold = xo;
xnew = zeros(n,1);

%Actual solution for determining error.
xtrue = A\b;

    for k = 1:max_iter
        for i = 1:n
            s = b(i);
            for j = 1:i-1
                s = s-A(i,j)*xnew(j);
            end
            for j = i+1:n
                s = s-A(i,j)*xold(j);
            end
            xnew(i) = s/A(i,i);
        end
        xold = xnew;
        error = norm(xnew-xtrue,"inf");
        disp(['iteration ', num2str(k),', max norm error:', num2str(error,8)]);
    end
    
end

    
