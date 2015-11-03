function subs = str2subs(shape, str)
% String representation of subscripts -> array of subscripts

subs = str2num(str);

assert(length(shape) == length(subs), 'Ndsparse:str2subs:SubsLengthMismatch', 'Length of subscripts array (%d) not equal to dimensions in shape(%d).', length(subs), length(shape))