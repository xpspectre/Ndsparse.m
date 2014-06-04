function subs = ind2subs(shape,ind)
%IND2SUBS Array of subscripts from linear index.
% SUBS = IND2SUBS(SHAPE,IND)
%
% Adapted from ind2sub from The MathWorks, Inc.

k = [1 cumprod(shape(1:end-1))];
subs = uint64(zeros(size(shape)));
for i = length(shape):-1:1,
    vi = rem(ind-1, k(i)) + 1;
    vj = (ind - vi)/k(i) + 1;
    subs(i) = vj;
    ind = vi;
end