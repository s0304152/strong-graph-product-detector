%输入参数                        输出参数  
%s为输入信号                     G为连通图对象
%K为fft点数                      Lap为图的laplacian矩阵
%N为量化电平数（图谱的顶点数）    second_eig为对应的laplacian矩阵的第二大特征值                    
%flag_1为标志，用于选择取全部的功率谱还是一半，0代表取全部的功率谱，1代表取对称功率谱的一半
function [G,Lap,second_eig] = signal2graph_newnew(s,K,N,flag_1)


%产生功率谱
switch flag_1
    case 0
        R = 1/K*((abs(fft(s,K))).^2);
    case 1
        R = 1/K*((abs(fft(s,K))).^2);
        R = R(end-K/2+1:end);
end
m = 0:length(R)-1;%功率谱的时间序列



%归一化后的功率谱
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


%转换成均匀图的邻接矩阵
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

%度矩阵
d = degree(G);
D = zeros(N);
for deg = 1:N
    D(deg,deg) = d(deg);
end

%laplacian矩阵
Lap = D-Adj;
if rank(Lap-Lap')~=0
    disp('存在不对称矩阵');
end

%laplacian矩阵的特征值
eigen = eig(Lap);
sort_eig = sort(eigen,'descend');
second_eig = sort_eig(2);

end

