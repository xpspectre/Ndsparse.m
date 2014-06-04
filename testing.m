% Basic N-dimensional array operations testing script
% Reference: http://www.mathworks.com/help/matlab/math/multidimensional-arrays.html
% Call first 3 dims: row x col x page
clear;clc;close all

A =        [1 2; 3 4; 5 6; 0 1];
A(:,:,2) = [7 8; 9 0; 1 2; 1 0];
A(:,:,3) = [3 4; 5 6; 7 8; 9 3];

B =        [5 7 8 0; 0 1 9 1; 4 3 6 2];
B(:,:,2) = [1 0 4 4; 3 5 6 2; 9 8 7 0];

Ax = tensor(A);
Bx = tensor(B);

As = sptensor(Ax);
Bs = sptensor(Bx);

C = ttt(As,Bs,1,2);
Cs = sptensor(C);

a = Ndsparse(A);
b = Ndsparse(B);

c = ttt(a,b,1,2);