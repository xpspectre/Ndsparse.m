% Basic N-dimensional array operations testing script
% Reference: http://www.mathworks.com/help/matlab/math/multidimensional-arrays.html
% Call first 3 dims: row x col x page
clear; clc; close all
rng('default');

% Add tprod to path if not already present for testing
pathCell = regexp(path, pathsep, 'split');
if ~ismember([pwd '/tprod'], pathCell)
    addpath([pwd '/tprod']);
end
clear pathCell

%% Basic testing
% x = [-1,-2,1,-3];
% y = [2,-1,-2];
% c1 = x(x < 0);
% c2 = y(y < 0);
% ismember(c1,c2)
% ismember(c2,c1)
% setdiff(c1,c2)
% setdiff(c2,c1)

% x = [1 2 3 4];
% n = 1e5;
% tic
% for i = 1:n
%     y = num2str(x);
% %     y = subs2ind(x, x); ~6x faster
% end
% toc
% 
% tic
% for i = 1:n
%     z = str2num(y);
% %     z = ind2subs(x, y);
% end
% toc

%% Dense tensor contraction using tprod
% A dims: [4,2,3]
A =        [1 2; 3 4; 5 6; 0 1];
A(:,:,2) = [7 8; 9 0; 1 2; 1 0];
A(:,:,3) = [3 4; 5 6; 7 8; 9 3];

% B dims: [3,4,2]
B =        [5 7 8 0; 0 1 9 1; 4 3 6 2];
B(:,:,2) = [1 0 4 4; 3 5 6 2; 9 8 7 0];

% C dims: [3,3]
C1 = tprod(A, [-1,-2,1], B, [2,-1,-2]);

C2 = tprod(A, [2,-2,-1], B, [-1,1,-2]);
C3 = tprod(B, [2,-2,-1], B, [1,-2,-1]);

%% Sparse tensor contraction using Ndsparse
a = Ndsparse(A);
b = Ndsparse(B);

assert(isequal(A, a.full))
assert(isequal(B, b.full))

c1 = ttt(a, [-1,-2,1], b, [2,-1,-2]);

c2 = ttt(a, [2,-2,-1], b, [-1,1,-2]);
c3 = ttt(b, [2,-2,-1], b, [1,-2,-1]);

assert(isequal(C1, c1.full))
assert(isequal(C2, c2.full))
assert(isequal(C3, c3.full))
