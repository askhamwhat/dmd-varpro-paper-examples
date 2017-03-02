
function [c,ax,lens] = zconf95(zdat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Obtain parameters of the 95 % confidence
% ellipse for the complex valued data
% in zdat
%
% Input
%
% zdat = column of complex valued data
%
% Output
%
% c = center of ellipse (mean of real
%   and imaginary parts)
% ax = axes of ellipse (each column is 
%   an axis)
% lens = the length of each axis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xdat = [real(zdat) imag(zdat)];
[n,p] = size(xdat);
c = mean(xdat);
[u,s,v] = svd((xdat-repmat(c,n,1))/sqrt(n-1));

F = fstat95(p,n-p);
T = p*(n-1)*F/(n-p);

ax = v;
lens = diag(s)*sqrt(T);

