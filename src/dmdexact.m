
function [w,e,atilde,varargout] = dmdexact(X,r,imode,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simple implementation of exact DMD 
% 
% Travis Askham
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (imode == 1)
    
    % exact (unprojected) DMD

    X1 = X(:,1:end-1);
    X2 = X(:,2:end);

    [u1,s1,v1] = svd(X1,'econ');

    u1 = u1(:,1:r);
    v1 = v1(:,1:r);
    s1 = s1(1:r,1:r);

    % similar --- but not equal --- to true A

    atilde = u1'*X2*(v1/s1);

    [w,e] = eig(atilde);

    % dmd eigenvalues

    e = diag(e);

    % dmd eigenvectors

    w = X2*(v1*(s1\w));

    if (nargout > 3)

        % build approximation of true A

        varargout{1} = w*diag(e)*pinv(w);
    end
    
elseif (imode == 2)
    
    % compute projected version
    
    [m,~] = size(X);
    
    if (nargin > 3)
        u = varargin{1};
        [mu,nu] = size(u);
        if (m ~= mu || nu ~= r) 
            error('input for projection of wrong dimensions .. bomb')
        end
    else
        [u,~,~] = svd(X,'econ');
        u = u(:,1:r);
    end
    
    % exact (but projected) DMD

    X1 = u'*X(:,1:end-1);
    X2 = u'*X(:,2:end);

    [u1,s1,v1] = svd(X1,'econ');

    u1 = u1(:,1:r);
    v1 = v1(:,1:r);
    s1 = s1(1:r,1:r);

    % similar --- but not equal --- to true A

    atilde = u1'*X2*(v1/s1);

    [w,e] = eig(atilde);

    % dmd eigenvalues

    e = diag(e);

    % dmd eigenvectors

    w = X2*(v1*(s1\w));
    w = u*w;

    if (nargout > 3)

        % build approximation of true A

        varargout{1} = w*diag(e)*pinv(w);
    end
    
    
end

end