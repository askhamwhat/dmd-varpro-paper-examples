function p = chebpol(x,n)
%
%   simple chebyshev polynomial evaluator
%   for input in [-1,1]. Returns T_n(x)
%
th = acos(x);
p = cos(n*th);

end