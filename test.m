% 22/5/11
% test the mainAlg.m

% clear;
% clc;

%% initialization
num = 1; % run 10 times
tarD = 2; % ¶—

% generate T=2^¶»
temp1 = (-10:0.5:0)';
temp2 = (0.1:0.1:2)';
temp3 = (2.5:0.5:10)';
temp = [temp1;temp2;temp3];
T = power(2,temp);

e = [0.01;0.03;0.05]; % epsilon

m = size(T,1);
n = size(e,1);
res = zeros(m*n*num,11); % save results
restot = zeros(m*n,8);

%% synthetic data ~ Low-Dimensional Manifold
% genidmdata.m
% example: 10000

%load("./data/A5000.mat","A1");
%load("./data/B5000.mat","B1");

%% usps
% example: 2000

%load("/home/ct/FNS/usps.mat","dd");
%A1 = dd(:,1:1000);
%B1 = dd(:,2001:3000);

% normalization
%A1 = double(A1);
%B1 = double(B1);

%A1 = mapminmax(A1,0,1);
%B1 = mapminmax(B1,0,1);

%% minist
% example: 2000

images = loadMNISTImages('./data/train-images-idx3-ubyte.txt');
labels = loadMNISTLabels('./data/train-labels-idx1-ubyte.txt');
A = images(:,labels==3);
B = images(:,labels==5);

n1 = size(A,2);
n2 = size(B,2);

A1 = A(:,1:500);
B1 = B(:,1:500);

%% cifar-10
%cifar_1 = load('./cifar-10/data_batch_1.mat');
%cifar_2 = load('./cifar-10/data_batch_2.mat');
%A = cifar_1.data;
%B = cifar_2.data;

%A1 = double(A(1:7500, :));
%B1 = double(B(1:7500, :));

% πÈ“ªªØ
%[A,PSA] = mapminmax(A1,0,1);
%[B,PSB] = mapminmax(B1,0,1);
%A1 = A';
%B1 = B';

T0=tic;
for eid = 1:n
    eid
    res((eid-1)*m*num+1:eid*m*num,10) = e(eid); % current epsilon
    restot((eid-1)*m+1:eid*m,8) = e(eid);
    for j = 1:num
        j
        % 1: sinkhorn, 2: fast network simplex, 3: fast emd
        [EMD0,T1,timeT] = mainAlg(tarD,e(eid),T,A1,B1,3); % T time T1
        res((eid-1)*m*num+(j-1)*m+1:(eid-1)*m*num+j*m,8) = EMD0; % original emd
        res((eid-1)*m*num+(j-1)*m+1:(eid-1)*m*num+j*m,9) = T1; % running time
        res((eid-1)*m*num+(j-1)*m+1:(eid-1)*m*num+j*m,1:7) = res((eid-1)*m*num+(j-1)*m+1:(eid-1)*m*num+j*m,1:7)+timeT;
        restot((eid-1)*m+1:eid*m,1:7) = restot((eid-1)*m+1:eid*m,1:7)+timeT;
    end
end

% compute the acc & res
%[acc,res,tab] = comacc(res);

res;
T5=toc;

%save("./result/minist2000SKHrou2res1.mat","res");
%save("./result/minist2000SKHrou2restot1.mat","restot");