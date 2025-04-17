%输入图对象和邻接矩阵
%输出度向量的熵

function [h] = degree_entropy(G)

d = degree(G);%提取图的度向量
d0 = d/sum(d);%归一化
d0(d0==0)=[];%删除度向量为0的值，避免无法计算熵
h = -sum(d0.*log2(d0));%熵

end