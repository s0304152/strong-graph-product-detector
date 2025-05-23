%输入参数                        输出参数  
%s为输入信号                     G为连通图对象,Rd为输出的参数c
%K为fft点数
%N为量化电平数（图谱的顶点数）
%gl为分组数
%flag_1为标志，用于选择取全部的功率谱还是一半，0代表取全部的功率谱，1代表取对称功率谱的一半
function [G,Adj] = bmmax_s2g_new(s,K,N,gl,flag_1)

%产生功率谱
switch flag_1
    case 0
        R = 1/K*((abs(fft(s,K))).^2);
    case 1
        R = 1/K*((abs(fft(s,K))).^2);
        R = R(end-K/2+1:end);
end

%采用分组的方法将功率谱分组，每gl个为一组，取每组最大组成新谱
%生成gl行，总点数除gl向上取整列的矩阵（即每一组为矩阵的一列）
group = zeros(gl,ceil(length(R)/gl));
for i = 1:length(R)
    group(i) = R(i);
end
R_new = zeros(1,ceil(length(R)/gl));
for j = 1:length(group)
    R_new(j) = max(group(:,j));
end
m = 0:length(R_new)-1;

%归一化
Rmax = max(R_new);
Rmin = min(R_new);
% Rd = Rmax-Rmin; % 计算最大最小值之差
U = (R_new-Rmin)/(Rmax-Rmin);%归一化


%量化(采用均匀量化)
%量化过程
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

%转换成图的邻接矩阵
Adj= zeros(N);
for h = 1:1:length(Qm)-1
    for g = 1:1:length(m)-h
        if Qm(g)~=Qm(g+h)%步长取所有
            Adj(Qm(g),Qm(g+h))=1;
            Adj(Qm(g+h),Qm(g))=1;
        end
    end
end

%创建连通图的对象
G = graph(Adj);


end

