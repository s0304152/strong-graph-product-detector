%����GΪͼ�Ķ�����֮ǰ�ĺ�������
%���GIΪginiϵ��
function [GI] = gini(G)
d = degree(G);%����ͼ��������Ķ�����
n = length(d);%�ܵĶ�����
b = sort(d,'ascend');%˳�����к�ĸ��������
T = sum(b);%�ܵĶ������
%�����ָ��
s = 0;
for i = 1:n
    s = s + b(i)/T*((n-i+1/2)/n);
end
GI = 1 - 2.*s;
% %����������
% v = zeros(1,n);
% for j = 1:n
%     for k = 1:j
%         v(j) = v(j) + b(k);
%     end
%     v(j) = v(j)/T;
% end
% v = [0,v];
% t = (0:n)/n;%����������ʱ������
% figure(1);
% plot([1,1],[0,1]);
% hold on
% plot([0,1],[0,1]);
% hold on
% value = spcrv([[t(1),t,t(end)];[v(1),v,v(end)]]);
% plot(value(1,:),value(2,:));
% sgtitle('����������')
% end