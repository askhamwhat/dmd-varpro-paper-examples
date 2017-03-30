function [w,e,atilde,varargout] = tlsdmd(X,r,imode,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simple implementation of total least squares DMD 
% 
% Travis Askham
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (imode == 1)
    
    % exact (unprojected) tlsDMD
    
    Z = [X(:,1:end-1); X(:,2:end)];

    [u1,~,~] = svd(Z,'econ');

    u11 = u1(1:end/2,1:r);
    u21 = u1(end/2+1:end,1:r);

    atilde = u21*pinv(u11);

    [w,e] = eig(atilde);

    e = diag(e);
    
    if (nargout > 3)
        varargout{1} = atilde;
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
    
    ux1 = u'*X;
    ux2 = ux1(:,2:end);
    ux1 = ux1(:,1:end-1);
            
    Z = [ux1; ux2];

    [u1,~,~] = svd(Z,'econ');

    u11 = u1(1:end/2,1:r);
    u21 = u1(end/2+1:end,1:r);

    atilde = u21*pinv(u11); %projected A

    [w,e] = eig(atilde);

    e = diag(e);
    
    if (nargout > 3)
        
        % unprojected A
        
        varargout{1} = u*atilde*u';
    end

    
end

end