function subs = str2subs(str)
% String representation of subscripts -> array of subscripts
%   Canonical string format is integers separate by a single space, with leading and
%   trailing whitespace trimmed. Here, string format is related, allowing any
%   amount of whitespece between indices and leading and trailing whitespace.
%
% Inputs:
%   str [ single string | n x 1 cell vector of strings ]
%       Indices in string format
%
% Outputs:
%   subs [ 1 x d vector of nonnegative integers | n x d matrix of nonnegative integers ]
%       Indices in vector format

if ischar(str)
    subs = str2num(str);
elseif iscellstr(str)
    subs = cell2mat(cellfun(@str2num, str, 'UniformOutput', false));
else
    error('Ndsparse:str2subs:InvalidStrArg', 'str must be a string or cell vector of strings')
end