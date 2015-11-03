function str = subs2str(subs)
% Array of subscripts -> string representation of subscripts
%   String format is integers separate by a single space, with leading and
%   trailing whitespace trimmed
%
% Inputs:
%   subs [ 1 x d vector of nonnegative integers | n x d matrix of nonnegative integers ]
%       Indices in vector format
%
% Outputs:
%   str [ single string | n x 1 cell vector of strings ]
%       Indices in string format

assert(isnumeric(subs), 'Ndsparse:subs2str:NonNumericSubs', 'subs must be a numeric vector or matrix')

subLine = @(subs) strtrim(sprintf('%d ', subs));

nSubs = size(subs,1);
if nSubs == 1
    str = subLine(subs);
else
    str = cellfun(subLine, num2cell(subs,2), 'UniformOutput', false);
end
