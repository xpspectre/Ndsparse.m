function results = unittests()
% Run all functions that are tests in this file
rng('default');

% Add tprod to path if not already present
pathCell = regexp(path, pathsep, 'split');
if ~ismember([pwd '/tprod'], pathCell)
    addpath([pwd '/tprod']);
end
clear pathCell

tests = functiontests(localfunctions);
if nargout < 1
    results = tests.run;
end
end

%% Test Ndsparse construction
function testConstructBlank(a)
x = Ndsparse();
a.verifyEqual(x.d, 0)
a.verifyTrue(isempty(x.shape))
a.verifyTrue(isempty(x.entries))
end

function testConstructScalar(a)
v = rand;
x = Ndsparse(v);
a.verifyEqual(x.d, 0)
a.verifyTrue(isempty(x.shape))
a.verifyEqual(x.entries('0'), v)
end

function testConstructListEntries(a)
x = makeNdsparse;
a.verifyEqual(x.d, 4)
a.verifyEqual(x.shape, ones(1,4)*30)
end

function testConstructFromExisting(a)
x = makeNdsparse;
y = Ndsparse(x);
a.verifyEqual(x.d, y.d);
a.verifyEqual(x.shape, y.shape);
a.verifyEqual(x,y); % value equality
a.verifyFalse(x == y) % handle equality - y must refer to different object that x
end

function testConstructFromDense(a)
X = makeDense;
x = Ndsparse(X);
a.verifyEqual(X, x.full)
end

%% Test utility functions
function testNdsparseEquality(a)
x = makeNdsparse;
y = x;
z = makeNdsparse;

% Handle equality
a.verifyTrue(x == y)
a.verifyFalse(x == z)

% Value equality - calls isEqual which checks fields recursively for handle objects
a.verifyEqual(x, z)
a.verifyEqual(x, y)
end

function testAddEntry(a)
x = makeNdsparse;
pos = [9,9,9,9];
val = pi;
x.addEntry(pos, val);
a.verifyEqual(x.getEntry(pos), val)
end

function testRemoveEntry(a)
x = makeNdsparse;
pos = [10 12 19 21]; % hardcode entry of makeNdsparse magic matrix
a.verifyEqual(x.getEntry(pos), 3)

x.removeEntry(pos);
a.verifyEqual(x.getEntry(pos), 0)
end

function testNnz(a)
x = makeNdsparse;
a.verifyEqual(x.nnz, 5)
end

%% Test operations
function testTranspose(a)
x = makeNdsparse;
y = Ndsparse(x); % transpose is inplace so copy first
permutation = [1,3,4,2];
y.transpose(permutation);

% Get an entry and make sure its transpose gives the same value
pos = [10 12 19 21]; % hardcode entry of makeNdsparse magic matrix
xval = x.getEntry(pos);
yval = y.getEntry(pos(permutation));

a.assertEqual(xval, yval)
end

function testContraction(a)
[X,Y] = makeDense;

% C1 dims: [3,3]
C1 = tprod(X, [-1,-2,1], Y, [2,-1,-2]);
C2 = tprod(X, [2,-2,-1], Y, [-1,1,-2]);
C3 = tprod(Y, [2,-2,-1], Y, [1,-2,-1]);

x = Ndsparse(X);
y = Ndsparse(Y);

a.verifyEqual(X, x.full)
a.verifyEqual(Y, y.full)

c1 = ttt(x, [-1,-2,1], y, [2,-1,-2]);
c2 = ttt(x, [2,-2,-1], y, [-1,1,-2]);
c3 = ttt(y, [2,-2,-1], y, [1,-2,-1]);

a.verifyEqual(C1, c1.full)
a.verifyEqual(C2, c2.full)
a.verifyEqual(C3, c3.full)
end

%% Helper functions
% Helper function that makes test object
function x = makeNdsparse
x = Ndsparse(magic(5), ones(1,4)*30);
end

function [X,Y] = makeDense
% X dims: [4,2,3]
X =        [1 2; 3 4; 5 6; 0 1];
X(:,:,2) = [7 8; 9 0; 1 2; 1 0];
X(:,:,3) = [3 4; 5 6; 7 8; 9 3];

% Y dims: [3,4,2]
Y =        [5 7 8 0; 0 1 9 1; 4 3 6 2];
Y(:,:,2) = [1 0 4 4; 3 5 6 2; 9 8 7 0];
end
