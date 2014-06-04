function ind = subs2ind(shape,subs)
%SUBS2IND Linear index from aray of subscripts.
%
% Adapted from sub2ind from The MathWorks, Inc.

k = [1 cumprod(shape(1:end-1))];
ind = uint64(1);
for i = 1:length(subs)
    v = subs(i);
    if (any(v(:) < 1)) || (any(v(:) > shape(i))) % Verify subscripts are within range
        error('Ndsparse:subs2ind:IndexOutOfRange','Index out of range');
    end
    ind = ind + (v-1)*k(i);
end