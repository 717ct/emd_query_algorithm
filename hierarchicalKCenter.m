% 22/5/11
% Hierarchical Gonzalez's algorithm

function [id,w,sub,cen]=hierarchicalKCenter(P, WP, num, r)
% inputs:
% P: points sets
% WP: A -> 1 B -> -1
% r: the max radius

% outputs:
% cen: center points id
% sub: subsets points of every cluster
% id: pi belongs to A or B
% w: weights of pi

Num=size(P,2); % the number of points in set P
iniNo=randi(Num); % the first center
disMat=sum(bsxfun(@minus,P,P(:,iniNo)).^2,1);

cen=zeros(num,1);
sub=ones(Num,1);
cen(1)=iniNo; % store the first center

for i=2:num
    i;
    [~,curNo]=max(disMat); % choose the fatherst point as the new cluster center
    cen(i)=curNo; % new cluster center
    tmpMat=sum(bsxfun(@minus,P,P(:,curNo)).^2,1);
    [disMat,tmpSeq]=min([disMat;tmpMat],[],1);% compare to the previous cluster inof, choose the min elements as row vector
    sub(tmpSeq==2)=i; % send points to new cluster
    j=i;
    ri=0;
    while j > 0
        %j;
        [tempr, ~] = appcomr(P(:,sub==j), P(:,curNo));
        ri = max(ri, tempr);
        j = j-1;
    end
    if ri <= r
        %i;
        ri;
        break;
    end
    
    if any(disMat)==0
        break;
    end
    
end

cen(cen==0)=[];
l=length(cen);
CSWA=zeros(1,l);
id=zeros(l,1);
w=zeros(l,1);

for i=1:l
    ithNo=(sub==i);
    CSWA(i)=sum(WP(ithNo)); % weight minus
    % cluster belongs to A if CSWA is positive, else B
    if CSWA(i)>=0
        id(i)=1;
        w(i)=CSWA(i);
    else
        id(i)=-1;
        w(i)=-CSWA(i);
    end
end