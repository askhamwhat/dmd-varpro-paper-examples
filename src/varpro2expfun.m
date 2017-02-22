
function A = varpro2expfun(alpha,t)

m = length(t);
n = length(alpha);

A = zeros(m,n);

ttemp = reshape(t,m,1);
atemp = reshape(alpha,n,1);

temp = ttemp*transpose(atemp);

A = exp(temp);

%for j = 1:n
%    for i = 1:m
%        A(i,j) = exp(alpha(j)*t(i));
%    end
%end

end
