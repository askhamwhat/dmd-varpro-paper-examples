
function indices = match_vectors(v1,v2)

v1 = v1(:);
v2 = v2(:);
costmat = abs(kron(v1.',ones(length(v2),1)) - kron(v2,ones(1,length(v1))));
[indices,cost] = munkres(costmat);

