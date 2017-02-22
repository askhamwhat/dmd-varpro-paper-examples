
function A = varpro2dexpfun(alpha,t,i)

m = length(t);
n = length(alpha);

if (i < 1 || i > n) 
    error('varpro2dexpfun: invalid index')
end

%A = zeros(m,n);
A = sparse(m,n);

ttemp = reshape(t,m,1);

A(:,i) = ttemp.*exp(alpha(i)*ttemp);

end
