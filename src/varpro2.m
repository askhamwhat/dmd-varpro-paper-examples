function [b,alpha,niter,err,imode,alphas] = varpro2stab(y,t,phi,dphi,m,n,is,ia, ...
    alpha_init,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variable projection algorithm for multivariate data
%
% Attempts a fit of the columns of y as linear combinations
% of the columns of phi(alpha,t), i.e.
%
% y_k = sum_j=1^n b_jk phi_j(alpha,t)
%
% Note that phi(alpha,t) is a matrix of dimension
% m x n where m is length (t) and n is number of columns.
%
% phi_j(alpha,t) is the jth column
% y_k is the kth column of the data
%
% Input:
%
% y - M x IS matrix of data
% t - M vector of sample times
% phi(alpha,t) - M x N matrix (or sparse matrix) valued 
%              function with input alpha
% dphi(alpha,t,i) - M x N matrix (or sparse matrix) valued
%                 function of alpha: returns the derivative 
%                 of the entries of phi with respect to the 
%                 ith component of alpha
% m - integer, number of rows of data/number of sample times
% n - integer, number of columns of phi
% is - integer, number of columns of data .. number of 
%      functions to fit
% ia - integer, dimension of alpha
% alpha_init - initial guess for alpha
% opts - options structure. See varpro_opts.m for details. Can
%   be created with default values via 
%       opts = varpro_opts();
%
% Output:
%
% b - N x IS matrix of coefficients .. each column gives
%     the coefficients for one of the functions (columns
%     of data) corresponding to the best fit
% alpha - N vector of values of alpha for best fit
% niter - number of iterations of the Marquardt algorithm
% err - the error for each iteration of the algorithm
% imode - failure mode
%            imode = 0, normal execution, tolerance reached
%            imode = 1, maxiter reached before tolerance
%            imode = 4, failed to find new search direction
%                       at step niter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% various error and warning string formats

mode8str = ['stall detected: residual reduced by less than %e' ...
    ' times residual at previous step. iteration %d' ...
    '. current residual %e'];
mode4str = ['failed to find appropriate step length at iteration %d' ...
    '. current residual %e'];
mode1str = ['failed to reach tolerance after maxiter = %d' ...
    ' iterations. current residual %e'];
optstr1 = ['input opts must be a structure, see varpro_opts.m.' ...
    ' Using default options ...'];
    
% set options, try to catch mistakes

opts_default = varpro_opts();

if (~isstruct(opts))
    warning(optstr1);
    opts = opts_default;
end

lambda0 = varpro2_getfield(opts,opts_default,'lambda0');
maxlam = varpro2_getfield(opts,opts_default,'maxlam');
lamup = varpro2_getfield(opts,opts_default,'lamup');
lamdown = varpro2_getfield(opts,opts_default,'lamdown');
ifmarq = varpro2_getfield(opts,opts_default,'ifmarq');
maxiter = varpro2_getfield(opts,opts_default,'maxiter');
tol = varpro2_getfield(opts,opts_default,'tol');
eps_stall = varpro2_getfield(opts,opts_default,'eps_stall');
iffulljac = varpro2_getfield(opts,opts_default,'iffulljac');

% initialize values

alpha = alpha_init;
alphas = zeros(length(alpha),maxiter);
djacmat = zeros(m*is,ia);
err = zeros(maxiter,1);
res_scale = norm(y,'fro');
scales = zeros(ia,1);

phimat = phi(alpha,t);
[U,S,V] = svd(phimat,'econ');
sd = diag(S);
tolrank = m*eps;
irank = sum(sd > tolrank*sd(1));
U = U(:,1:irank);
S = S(1:irank,1:irank);
V = V(:,1:irank);
b = phimat\y;
res = y - phimat*b;
errlast = norm(res,'fro')/res_scale;

imode = 0;

for iter = 1:maxiter
    
    % build jacobian matrix, looping over alpha indeces
    
    for j = 1:ia
        dphitemp = dphi(alpha,t,j);
        djaca = (dphitemp - sparse(U*(sparse(U'*dphitemp))))*b;
        if (iffulljac == 1)
            % use full expression for Jacobian
            djacb = U*(S\(V'*(sparse(dphitemp'*res))));
            djacmat(:,j) = -(djaca(:) + djacb(:));
        else
            % use approximate expression
            djacmat(:,j) = -djaca(:);
        end
        % the scales give the "marquardt" part of the algo.
        scales(j) = 1;
        if (ifmarq == 1)
            scales(j) = min(norm(djacmat(:,j)),1);
            scales(j) = max(scales(j),1e-6);
        end
    end
    
    % loop to determine lambda (lambda gives the "levenberg" part)
        
    % first check if current or shrunk version works
    
    dtempmat = [djacmat; lambda0*diag(scales)];
    rhs = [res(:); zeros(ia,1)];
    
    % get step
    
    delta0 = dtempmat\rhs;
    
    % new alpha guess
    
    alpha0 = alpha - delta0;
    
    % corresponding residual
    
    phimat = phi(alpha0,t);
    b0 = phimat\y;
    res0 = y-phimat*b0;
    err0 = norm(res0,'fro')/res_scale;
    
    % check if this is an improvement
    
    if (err0 < errlast) 

        % see if a smaller lambda is better
        
        lambda1 = lambda0/lamdown;
        dtempmat = [djacmat; lambda1*diag(scales)];
        delta = dtempmat\rhs;    
        alpha1 = alpha - delta;
        phimat = phi(alpha1,t);
        b1 = phimat\y;
        res1 = y-phimat*b1;
        err1 = norm(res1,'fro')/res_scale;
        
        if (err1 < err0)
            lambda0 = lambda1;
            alpha = alpha1;
            errlast = err1;
            b = b1;
            res = res1;
        else
            alpha = alpha0;
            errlast = err0;
            b = b0;
            res = res0;
        end
    else
    % if not, increase lambda until something works
    % this makes the algorithm more and more like gradient descent
    
        for j = 1:maxlam
        
            lambda0 = lambda0*lamup;
            dtempmat = [djacmat; lambda0*diag(scales)];

            rhs = [res(:); zeros(ia,1)];
            
            delta0 = dtempmat\rhs;
            alpha0 = alpha - delta0;

            phimat = phi(alpha0,t);
            b0 = phimat\y;
            res0 = y-phimat*b0;
            err0 = norm(res0,'fro')/res_scale;
            
            if (err0 < errlast) 
                break
            end

        end
        
        if (err0 < errlast) 
            alpha = alpha0;
            errlast = err0;
            b = b0;
            res = res0;
        else
            
            % no appropriate step length found
            
            niter = iter;
            err(niter) = errlast;
            imode = 4;
            warning(mode4str,iter,errlast);
            return
        end
    end
    
    alphas(:,iter) = alpha;
    
    err(iter) = errlast;
    if (errlast < tol)
        
        % tolerance met
        
        niter = iter;
        return;
    end
    
    if (iter > 1)
        if (err(iter-1)-err(iter) < eps_stall*err(iter-1))
            
            % stall detected
            
            niter = iter;
            imode = 8;
            warning(mode8str,eps_stall,iter,errlast);
            return;
        end
    end
    
    phimat = phi(alpha,t);
    [U,S,V] = svd(phimat,'econ');
    sd = diag(S);
    irank = sum(sd > tolrank*sd(1));
    U = U(:,1:irank);
    S = S(1:irank,1:irank);
    V = V(:,1:irank);
    
end

% failed to meet tolerance in maxiter steps

niter = maxiter;
imode = 1;
warning(mode1str,maxiter,errlast);

end

function out = varpro2_getfield(opts,opts_default,in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get value of field from struct if it exists,
% otherwise set to default value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optstr2 = 'opts struct is missing %s field, using default';
optstr3 = 'opts default struct is missing %s field! bomb';

if (isfield(opts,in))
    out = opts.(in);
else
    warning(optstr2,in);
    if (isfield(opts_default,in))
        out = opts_default.(in);
    else
        error(optstr3,in);
    end
end

end
