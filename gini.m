%输入G为图的对象，由之前的函数生成
%输出GI为gini系数
function [GI] = gini(G)
d = degree(G);%生成图各个顶点的度序列
n = length(d);%总的顶点数
b = sort(d,'ascend');%顺序排列后的各顶点度数
T = sum(b);%总的顶点度数
%求基尼指数
s = 0;
for i = 1:n
    s = s + b(i)/T*((n-i+1/2)/n);
end
GI = 1 - 2.*s;
% %洛伦兹曲线
% v = zeros(1,n);
% for j = 1:n
%     for k = 1:j
%         v(j) = v(j) + b(k);
%     end
%     v(j) = v(j)/T;
% end
% v = [0,v];
% t = (0:n)/n;%洛伦兹曲线时间序列
% figure(1);
% plot([1,1],[0,1]);
% hold on
% plot([0,1],[0,1]);
% hold on
% value = spcrv([[t(1),t,t(end)];[v(1),v,v(end)]]);
% plot(value(1,:),value(2,:));
% sgtitle('洛伦兹曲线')
% end