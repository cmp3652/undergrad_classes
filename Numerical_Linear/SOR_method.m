% Clara Pitkins
% cmp3652@rit.edu
% Numerical Linear Algebra
% Succesive Over Relaxation Method (varient of Gauss-Seidel Method)
% 2023

% Initialize arrays and variables
A = [-3,-1,0,0,0,0,0,0,0,0,0,1/2;-1,3,-1,0,0,0,0,0,0,0,1/2,0;0,-1,3,-1,0,0,0,0,0,1/2,0,0;0,0,-1,3,-1,0,0,0,1/2,0,0,0;0,0,0,-1,3,-1,0,1/2,0,0,0,0;
    0,0,0,0,-1,3,-1,0,0,0,0,0;0,0,0,0,0,-1,3,-1,0,0,0,0;0,0,0,0,1/2,0,-1,3,-1,0,0,0;0,0,0,1/2,0,0,0,-1,3,-1,0,0;0,0,1/2,0,0,0,0,0,-1,3,-1,0;
    0,1/2,0,0,0,0,0,0,0,-1,3,-1;1/2,0,0,0,0,0,0,0,0,0,-1,3];
b = [2.5;1.5;1.5;1.5;1.5;1;1;1.5;1.5;1.5;1.5;2.5];
xo = zeros(12,1);
w = 1.1;
max_iter = 10;

% Run Function
xnew = sor_method(A,b,xo,w,max_iter)

% Function
function xnew = sor_method(A, b, xo, w, max_iter)
n = size(A,1);
xold = xo;
xnew = zeros(n,1);
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
            xnew(i) = (s*w)/A(i,i) + (1-w)*xold(i);
        end
        xold = xnew;
        error = norm(xnew-xtrue,"inf");
        disp(['iteration ', num2str(k),', max norm error:', num2str(error,8)]);
    end
    
end

    
