# emd_query_algorithm

#### 1. data sources: USPS (http://www.cad.zju.edu.cn/home/dengcai/Data/MLData.html), MINIST (http://www.cad.zju.edu.cn/home/dengcai/Data/MLData.html), CIFAR-10 (http://www.cs.toronto.edu/~kriz/cifar.html).
#### 2. test.m: initialize parameters, run the EMD Query Algorithm on datasets, and save the results.
#### 3. mainAlg.m, hierarchicalKCenter.m: main code of the EMD Query Algorithm.
#### 4. appcomr.m: randomly select a point, and compute the radius approximately.
#### 5. norCol.m: matrix normalization.
#### 6. distance.m: compute the distance matrix.
#### 7. Sinkhorn.m, Transport.m: Sinkhorn distance.
#### 8. loadMNISTImages.m, loadMNISTLabels.m: load MINIST datasets.
#### 9. emd_hat_mex.mexw64, emd_hat_mex.mexa64, emd_hat_mex.mexmaci64, emd_hat_mex.m, emd_hat_mex_nes.m, emd_hat_mex.cxx: FastEMD (using C++ code with MATLAB on windows, linux and macos).
#### 10. maintest.mexw64, maintest.mexa64: Network Simplex (using C++ code with MATLAB).
