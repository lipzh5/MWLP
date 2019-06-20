function label_sequence = MWLP(M,A,lambda)
% Motif-based weighted label propagation 
% Author: Pei-Zhen Li
% Date: 2019-05-31
% Input: M, the motif adjacency matrix
%        uA, the adjacency matrix, which is symmectric
% Output: L, final label list of nodes

%n = size(M,1); % number of nodes
%M = MotifAdjacency(sparse(uA),'M4'); %motif adjacency matrix
W = M + A;
%lambda = 0.3;
%W = uA;
n = size(W,1);
%*******Initialization********
label = 1:n; %初始化label
MAX_iter = 6;  % Needs to revise later
[e1,e2,s] = find(W);  %symmectric-
%***Pre-processing*****
e = unique([e1,e2]);  % non-isolated node lists
e = sort(e);
num = numel(e); 
nbr_list = cell(num,1);
for i = 1:num
    u = e(i);
    idx = find(e1==u);
    vv = e2(idx); 
    vv = unique(vv);  % a set of unique labels;
    nbr_list{i} = vv;
end
% neighbor关系在每次迭代中不会发生变化
%****nbr_list stores the neighbor list****%
label_sequence = zeros(n,MAX_iter); %保存所有迭代的label结果
for t = 1:MAX_iter
   for i = 1:num  %遍历每个非孤立节点
       %label_tmp = label;
       u = e(i);
       vv = nbr_list{i};
       lvv = label(vv); % labels of vv-这步应该先保存好类标
       num_nbr = numel(lvv); % number of neighbors
       Vscores = zeros(num_nbr,1);
       for j=1:num_nbr
           v = vv(j);
           %lv = label(v); %节点v的label
           lv = lvv(j); %节点v的label 
           C_label = numel(find(label==lv)); %和v相同label的节点数量
           %disp(C_label);
           %C_labels(j) = C_label; %保存每个节点的authority：和它相同label的节点数目
           score = lambda*C_label + (1-lambda)*W(u,v); 
           Vscores(j) = score;
       end
       %disp('Vscores:');
       %disp(Vscores);
       idx = find(Vscores==max(Vscores)); %在neighbor list中的index
       if numel(idx)>1
           %r = randi([1, numel(idx)],1,1);
           r = randi(numel(idx));
           index = idx(r); %tie-breaker:random select
       else
           index = idx;
       end
       picked_label = lvv(index);
       label(u) = picked_label; %更新label-异步更新
   end
   label_sequence(:,t) = label;
   %L = label;
   if t>10
       if label_sequence(:,t)==label_sequence(:,t-1)
           disp('Converge before MAX_iter is reached!');
           %disp(t);
           break;
       end
   end
end
L = label_sequence(:,1);  %选取第**次迭代后的结果
save('tmp.mat','label_sequence');
end





