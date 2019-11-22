% Main function for MWLP 
% 2019-05-31
clc;
clear all;

load('data\Cora.mat');
M = MotifAdjacency(sparse(uA),'M4'); %motif adjacency matrix
n = size(M,1); % number of nodes
%lambda = 0.5;
Iter = 10;  % iterations 
para_cluster = cell(10,1);
cnt = 1;
for lambda = 0.1:0.1:1.0  % lambda is a parameter
    disp(lambda);
    clusters = zeros(n, Iter);
    for t = 1:Iter
        clu_sequence = MWLP_new(M,uA,lambda);
        clusters(:,t) = clu_sequence(:,1)'; % needs to revise for different datasets
    end  
    para_cluster{cnt} = clusters;
    cnt = cnt+1;
end

save('result_MWLP\cluster.mat','para_cluster');

%*******Debug*******%
% lambda = 0.8;
% clu = MWLP_new(M,uA,lambda);
