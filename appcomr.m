% 22/5/11
% compute the approximate radius

function [r, p0] = appcomr(P, p0)
% random select a point Po as a root node 
n=size(P,2);
if nargin==1
    %n=size(P,2);
    R1=randi(n);
    p0=P(:,R1);
end
% approximateLy compute ¡÷~
tempP=sum(bsxfun(@minus,P,p0).^2,1);
r=max(tempP,[],2); % furthest point
r=bsxfun(@power,r,0.5);
%[~,r] = MEB(P,n1+n2,10);

