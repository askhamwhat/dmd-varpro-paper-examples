function [xr,err] = best_reconstruction(x,e,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function finds the best reconstruction possible
% of the data x given the eigenvalues e and times t
%
% Input:
%
% x - data to approximate
% e - eigenvalues 
% t - sample times
%
% Output:
%
% Let Phi = exp(e*t') where e and t are column vectors.
% Then
%
% xr - each row is given as a best fit linear combo
%   of the rows of Phi
% err - the value |x-xr|_F/|x|_F, where |.|_F denotes 
%   the Frobenius norm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Phi = varpro2expfun(e,t);

b = Phi\transpose(x);

xr = transpose(Phi*b);

err = norm(x-xr,'fro')/norm(x,'fro');