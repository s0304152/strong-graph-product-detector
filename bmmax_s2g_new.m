%�������                        �������  
%sΪ�����ź�                     GΪ��ͨͼ����,RdΪ����Ĳ���c
%KΪfft����
%NΪ������ƽ����ͼ�׵Ķ�������
%glΪ������
%flag_1Ϊ��־������ѡ��ȡȫ���Ĺ����׻���һ�룬0����ȡȫ���Ĺ����ף�1����ȡ�Գƹ����׵�һ��
function [G,Adj] = bmmax_s2g_new(s,K,N,gl,flag_1)

%����������
switch flag_1
    case 0
        R = 1/K*((abs(fft(s,K))).^2);
    case 1
        R = 1/K*((abs(fft(s,K))).^2);
        R = R(end-K/2+1:end);
end

%���÷���ķ����������׷��飬ÿgl��Ϊһ�飬ȡÿ������������
%����gl�У��ܵ�����gl����ȡ���еľ��󣨼�ÿһ��Ϊ�����һ�У�
group = zeros(gl,ceil(length(R)/gl));
for i = 1:length(R)
    group(i) = R(i);
end
R_new = zeros(1,ceil(length(R)/gl));
for j = 1:length(group)
    R_new(j) = max(group(:,j));
end
m = 0:length(R_new)-1;

%��һ��
Rmax = max(R_new);
Rmin = min(R_new);
% Rd = Rmax-Rmin; % ���������Сֵ֮��
U = (R_new-Rmin)/(Rmax-Rmin);%��һ��


%����(���þ�������)
%��������
Qm = m;
for c = 1:length(m)
    for d = 0:N-1
        if U(c)>=d/N && U(c)<(d+1)/N
            Qm(c) = d + 1;
            break
        elseif U(c)==1
             Qm(c) = N;
        end
    end
end

%ת����ͼ���ڽӾ���
Adj= zeros(N);
for h = 1:1:length(Qm)-1
    for g = 1:1:length(m)-h
        if Qm(g)~=Qm(g+h)%����ȡ����
            Adj(Qm(g),Qm(g+h))=1;
            Adj(Qm(g+h),Qm(g))=1;
        end
    end
end

%������ͨͼ�Ķ���
G = graph(Adj);


end

