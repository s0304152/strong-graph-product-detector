%�������                        �������  
%sΪ�����ź�                     GΪ��ͨͼ����
%KΪfft����
%NΪ������ƽ����ͼ�׵Ķ�������
%glΪ������
%flag_1Ϊ��־������ѡ��ȡȫ���Ĺ����׻���һ�룬0����ȡȫ���Ĺ����ף�1����ȡ�Գƹ����׵�һ��
function [G,Adj] = range_s2g_new(s,K,N,gl,flag_1)

%����������
switch flag_1
    case 0
        R = 1/K*((abs(fft(s,K))).^2);
    case 1
        R = 1/K*((abs(fft(s,K))).^2);
        R = R(end-K/2+1:end);
end

%���÷���ķ����������׷��飬ÿgl��Ϊһ�飬ȡÿ�������С֮���������
%����gl�У��ܵ�����gl����ȡ���еľ��󣨼�ÿһ��Ϊ�����һ�У�
group1 = zeros(gl,ceil(length(R)/gl));
group2 = zeros(gl,ceil(length(R)/gl));
group2(:,end) = 10000;%ȡһ���㹻�������ֹ���һ����Сֵȡ����
for i = 1:length(R)
    group1(i) = R(i);
    group2(i) = R(i);
end
R_new1 = zeros(1,ceil(length(R)/gl));
R_new2 = zeros(1,ceil(length(R)/gl));
for j1 = 1:length(group1)
    R_new1(j1) = max(group1(:,j1));
end
for j2 = 1:length(group2)
    R_new2(j2) = min(group2(:,j2));
end
R_new = R_new1-R_new2;
m = 0:length(R_new)-1;


%��һ��
Rmax = max(R_new);
Rmin = min(R_new);
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

