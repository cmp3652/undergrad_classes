% Clara Pitkins
% cmp3652@rit.edu
% Numerical Linear Algebra Project
% Singular Value Decomposition for Image
% 5/2023
% 
% Download 'Rover2.tiff' from Github

% read the image
rover = imread('Rover2.tiff');
rover = double(rover);
rover = rover(:,:,1);

% Display image with scaled colors
imagesc(rover);
colormap(gray);

% do SVD
[U, S, V] = svd(rover);
%display rank
r = rank(rover)

% number of modes
m = 486;
n = 608;
p = 243;
figure(1);
rover1 = U(:,1:p)*S(1:p,1:p)*V(:,1:p)';
imagesc(rover1);
colormap(gray);
title(['p = ', num2str(p)])
ratio = (p+p*m+p*n)/(r+r*n+r*m)