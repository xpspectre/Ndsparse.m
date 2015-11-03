function str = subs2str(shape, subs)
% Array of subscripts -> string representation of subscripts

assert(length(shape) == length(subs), 'Ndsparse:subs2str:SubsLengthMismatch', 'Length of subscripts array (%d) not equal to dimensions in shape(%d).', length(subs), length(shape))

str = num2str(subs);