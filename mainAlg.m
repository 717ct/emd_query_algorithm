% 22/5/11
% EMD Query algorithm

%clc
%clear
function [EMD0,T1,timeT] = mainAlg(tarD,e,T,A,B,alg)
%% initialization
timeT = zeros(size(T,1),7);
timeT(:,2) = T;
TT = 0;

n1 = size(A,2);
n2 = size(B,2);
n = n1+n2;
%d = size(A,1);
WA = ones(1,n1);
WB = ones(1,n2);

P = [A B];

WP = [WA/sum(WA) -WB/sum(WB)]; % normalize sum of weights

%% compute original emd
tic
disp("original emd:");
if alg == 1
    [EMD0,~] = Sinkhorn(A,B,WA,WB)
elseif alg == 2
    D = distance(A,B);
    EMD0 = maintest(D,WA/sum(WA),WB/sum(WB))
elseif alg == 3
    D = distance(P);
    D = triu(D);
    D(1:n1,1:n1) = 0;
    D(n1+1:n,n1+1:n) = 0;
    extra_mass_penalty = 0;
    flowType = 3;
    WAF = [(1/n1)*ones(n1,1);zeros(n2,1)];
    WBF = [zeros(n1,1);(1/n2)*ones(n2,1)];
    [EMD0] = demo_FastEMD_compute(WAF,WBF,D,extra_mass_penalty,flowType)
end
T1=toc;

%% random select a point Po as a root node and compute ¡÷~
tic
[rA, ~] = appcomr(A);
[rB, ~] = appcomr(B);
delta_0 = max(rA, rB);

%[~, P0] = appcomr(P);
%[~,r] = MEB(P, n1+n2, 10)

%% build tree and visualize its sub-trees
minballs = bsxfun(@power,2,2*tarD); % l <= 2^2¦Ñ balls
cluSub = ones(n1+n2,1); % i=1, cluster 1
cluNum = 1; % i=1, 1 cluster

for i = 1:floor(log2(2/e))+5
%for i=1:3
    % i: run Gonzalez's alg for the i-th layer clusters, and generate
    % add clusters for the i+1-th layer
    i
    idy = 1;
    % radius or dimeter of ball
    delta_i = bsxfun(@rdivide, delta_0, bsxfun(@power, 2, i-3));
    
    cluNum = min(cluNum, n1+n2);
    
    recen = zeros(min(minballs,n1+n2),1); % clusters centers id
    resub = ones(n1+n2,1); % clusters id
    reid = zeros(min(minballs,n1+n2),1); % clusters class: 1->A -1->B
    rew = zeros(min(minballs,n1+n2),1); % clusters centers weight
        
    for j = 1:cluNum
        idn = find(cluSub == j);
        if length(idn) == 1
            %disp("1 point!");
            recen(idy) = idn; % note! center id is different in subset!
            resub(idn) = idy; % one layer
            rew(idy) = abs(WP(idn));
            if WP(idn) > 0
                reid(idy) = 1;
            else
                reid(idy) = -1;
            end
            idy = idy+1;
        continue;
        end
        
        [id, w, sub, cen] = hierarchicalKCenter(P(:,idn), WP(idn), minballs, delta_i);
        l = length(cen);
        recen(idy:idy+l-1) = idn(cen); % note! center id is different in subset!
        resub(idn) = sub+idy-1; % one layer
        reid(idy:idy+l-1) = id;
        rew(idy:idy+l-1) = w;
        
        idy = idy+l;
    end
    cluSub = resub;
    cluSub(cluSub == 0) = [];
    cluNum = idy-1; % add idy-1 new clusters
        
    AA = recen(reid == 1);
    BB = recen(reid == -1);
    WAA = rew(reid == 1);
    WBB = rew(reid == -1);    
       
    AA(WAA == 0) = [];
    BB(WBB == 0) = [];    
    WAA(WAA == 0) = [];
    WBB(WBB == 0) = [];
    if ~isempty(WAA) && ~isempty(WBB)
        if alg == 1
            [EMD1,~]=Sinkhorn(P(:,AA),P(:,BB),WAA',WBB');
        elseif alg == 2
            DD=distance(P(:,AA),P(:,BB));
            EMD1=maintest(DD,WAA/sum(WAA),WBB/sum(WBB));
        elseif alg == 3
            nn1=size(AA,1);
            nn2=size(BB,1);
            PP=[P(:,AA) P(:,BB)]; % fastemd
            nn=nn1+nn2;

            WAF=[WAA; zeros(nn2,1)];
            WBF=[zeros(nn1,1); WBB];
            
            DD=distance(PP);
            DD=triu(DD);
            DD(1:nn1,1:nn1)=0;
            DD(nn1+1:nn,nn1+1:nn)=0;
            [EMD1]=demo_FastEMD_compute(WAF,WBF,DD,extra_mass_penalty,flowType);
        end
    else
        disp("1 cluster A or B!");
    end
    T2=toc;
    
    case1=(EMD1-delta_i)/EMD0;
    case2=(EMD1+delta_i)/EMD0;
       
    m=numel(T);
    % check ending condition
    if m > 0
        for t =1:m
            if T(t) < case1
                timeT(timeT(:,2)==T(t),5)=EMD1; % current EMD
                timeT(timeT(:,2)==T(t),1)=T2; % running time
                timeT(timeT(:,2)==T(t),3)=1; % case
                timeT(timeT(:,2)==T(t),6)=i; % output layer
                timeT(timeT(:,2)==T(t),7)=case1; % case1
                T(t)=0;
            elseif T(t) > case2
                timeT(timeT(:,2)==T(t),5)=EMD1;
                timeT(timeT(:,2)==T(t),1)=T2;
                timeT(timeT(:,2)==T(t),3)=2;
                timeT(timeT(:,2)==T(t),6)=i;
                timeT(timeT(:,2)==T(t),7)=case2; % case3
                T(t)=0;
            %elseif abs(EMD1-EMD0)/EMD0 < min(0.01,power(e,3))
            elseif abs(EMD1-EMD0)/EMD0 < 0.001
                timeT(timeT(:,2)==T(t),5)=EMD1;
                timeT(timeT(:,2)==T(t),3)=3;
                timeT(timeT(:,2)==T(t),6)=i;
                timeT(timeT(:,2)==T(t),7)=abs(EMD1-EMD0)/EMD0; % case3
                T(t)=0;
            end
        end
        T(T==0)=[]; % delete T that already outputs case
        T3=toc;
        TT=T3-T2+TT;
    else
        break;
    end   
end

T4=toc;
timeT(timeT(:,1)==0,6)=i;
timeT(timeT(:,1)==0,3)=3;
timeT(timeT(:,1)==0,1)=T4-TT;
timeT(:,4)=timeT(:,1)./T1;