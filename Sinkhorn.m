function [EMD,FM]=Sinkhorn(A,B,WA,WB)
% Squared Euclidean distance matrix
WA=WA/sum(WA);
WB=WB/sum(WB);
disMat=distance(A,B);
disMat(disMat<0)=0;
% disMat=sqrt(disMat);
lambda=60/median(disMat(:));
K=exp(-lambda*disMat);
U=K.*disMat;
[EMD,~,u,v]=Transport(WA',WB',K,U,lambda);
FM=bsxfun(@times,v',bsxfun(@times,u,K));