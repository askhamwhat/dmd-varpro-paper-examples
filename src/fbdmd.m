
function [w,e,atilde,varargout] = fbdmd(X,r,imode,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simple implementation of forward-backward DMD 
% 
% Travis Askham
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ifoverride = 0;

if (nargin > 4)
    ifoverride = varargin{2};
end

if (imode < 3)

    X1 = X(:,1:end-1);
    X2 = X(:,2:end);

    [u1,s1,v1] = svd(X1,'econ');

    u1 = u1(:,1:r);
    v1 = v1(:,1:r);
    s1 = s1(1:r,1:r);

    [u2,s2,v2] = svd(X2,'econ');

    u2 = u2(:,1:r);
    v2 = v2(:,1:r);
    s2 = s2(1:r,1:r);

    % exact DMD (reduced) matrix forward and backward

    fatilde = u1'*X2*v1/s1;

    batilde = u2'*X1*v2/s2;
    
else
    
    [u,~,~] = svd(X,'econ');
    
    u = u(:,1:r);
    
    xu = u'*X;
    
    X1 = xu(:,1:end-1);
    X2 = xu(:,2:end);

    [u1,s1,v1] = svd(X1,'econ');

    u1 = u1(:,1:r);
    v1 = v1(:,1:r);
    s1 = s1(1:r,1:r);

    [u2,s2,v2] = svd(X2,'econ');

    u2 = u2(:,1:r);
    v2 = v2(:,1:r);
    s2 = s2(1:r,1:r);

    % exact DMD (reduced) matrix forward and backward

    fatilde = u1'*X2*v1/s1;

    batilde = u2'*X1*v2/s2;
    
end

if (imode == 1)

    % dawson et. al formula for atilde
    
    atilde = sqrtm(fatilde*pinv(batilde));
    [atilde,errtemp] = besterrnbyn(fatilde,atilde,ifoverride); %sqrt is non-unique
    
    [w,e] = eig(atilde);
    e = diag(e);
    
elseif (imode == 2 || imode == 3)
    
    % use full version of forward and backward matrices
    % for imode = 3 the data has been projected
    
    [wf,ef] = eig(fatilde);
    wf = X2*(v1*(s1\wf));

    [wb,eb] = eig(batilde);
    wb = X1*(v2*(s2\wb));

    fatilde = wf*ef*pinv(wf);
    batilde = wb*eb*pinv(wb);

    atilde = sqrtm(fatilde*pinv(batilde));
    [atilde,errtemp] = besterrnbyn(fatilde,atilde,ifoverride);

    [w,e] = eig(atilde);
    e = diag(e);
    
end

end

function [Cbest,besterr] = besterrnbyn(A,B,ifoverride)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find modification of B which is closest to
% A. We assume that B is computed as a square root.
% The nonuniqueness of matrix square roots
% is up to +/- times each eigenvalue (2^n 
% possibilities for an n x n matrix)
%
% Warning: this is a slow routine! It is exponential
% in n and should not be called for large n
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,n] = size(A);

if (n > 25 && ifoverride ~= 1) 
    error('this will take forever. bomb')
end

C = B;
besterr = norm(A-C,'fro');
Cbest = C;

[w,e] = eig(B);
de = diag(e);
pw = pinv(w);

for i = 1:2^n-1
    v = (-1).^str2num(dec2bin(i,n)');
    C = w*diag(de.*v)*pw;
    err = norm(A-C,'fro');
    if (err < besterr)
        besterr = err;
        Cbest = C;
    end
end

end