%输入参数                        输出参数  
%s为输入信号                     G为连通图对象
%K为fft点数                      dmax为顶点的最大度
%N为量化电平数（图谱的顶点数）    zeroeig_count为laplacian矩阵零特征值个数                    
function [G,dmax,zeroeig_count] = corr_s2g_zeroeig_new(s,k,N)

s_mean = mean(s);%求s均值
snew = s-s_mean;%去掉均值
%自相关
[R,m] = xcorr(snew,k-1);
R = R(end-k+1:end);%取对称的自相关函数的一半
m = m(end-k+1:end);
%归一化
Rmax = max(R);
Rmin = min(R);
U = (R-Rmin)/(Rmax-Rmin);


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
    for g = 1:1:length(Qm)-h
        if Qm(g)~=Qm(g+h)%步长取所有
            Adj(Qm(g),Qm(g+h))=1;
            Adj(Qm(g+h),Qm(g))=1;
        end
    end
end

%画出连通图
G = graph(Adj);%创建连通图对象

%求所有顶点的最大度
d1 = degree(G);
dmax = max(d1);

%度矩阵
D = zeros(N);
for deg = 1:N
    D(deg,deg)=d1(deg);
end

%laplacian矩阵
 Lap = D-Adj;
%  Lap = laplacian(G);

%0特征值
eigen = eig(Lap);
eigen = round(eigen);
zeroeig_position = find(eigen == 0);
zeroeig_count = length(zeroeig_position);

end

