% 22/5/11
% generate the synthetic data

%% Low-Dimensional Manifold
d = 500;
n1 = 1000;
n2 = 1000;
tarD = 5; % the dimensionality of the manifold

disp('two different manifolds');
W = randn(d);
A = zeros(d,n1);
B = zeros(d,n2);
pA = randn(tarD,n1);
pB = randn(tarD,n2);
R1 = randi(d,1,tarD);
R2 = randi(d,1,tarD);
A(R1,:) = pA;
B(R2,:) = pB;
A = W*A;
B = W*B;
