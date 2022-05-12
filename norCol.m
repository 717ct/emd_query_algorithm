% 22/5/11
% Normalize each column vector in matrix X

function Y=norCol(X)

colNor=sqrt(sum(X.^2));
Y=bsxfun(@rdivide,X,colNor);