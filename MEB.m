function [o,r] = MEB(Q,num,rmeb)

init=randi(num);
%rmeb=10;

o=Q(:,init);
for i = 1:rmeb
    [~,index]=max(sum(bsxfun(@power,bsxfun(@minus,Q,o),2)));
    p=Q(:,index);
    c=(o*i+p)/(i+1);
    o=c;
end
r=2*bsxfun(@power,max(sum(bsxfun(@power,bsxfun(@minus,Q,o),2))),0.5);

end

