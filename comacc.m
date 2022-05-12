% 22/5/11
% compute the acc % std

function [acc,res,tab]=comacc(res)
% EMD > T case 1
res(:,12)=(res(:,8)>res(:,2));
res(res(:,12)==0,12)=2;
res(:,13)=(res(:,3)==(res(:,12)));
tab=tabulate(res(:,13))
acc=tab(2,3)