%用2020文中的方法判决
function [G,E]=signal2graph_2020(r,N,theta)

%对信号作量化
%量化过程
u=r;
Qm = zeros(1,length(u));
for c = 1:length(u)
    if u(c) < theta*sqrt(2)*erfinv(2/N-1)
       Qm(c) = 1;
    else 
        for delta = 2:N
            if u(c)>=theta*sqrt(2)*erfinv(2*(delta-1)/N-1)&&u(c)<=theta*sqrt(2)*erfinv(2*(delta)/N-1)
                Qm(c) = delta;
                break
            end
        end
    end
end
m = 0:length(Qm)-1;

%转换成图的邻接矩阵
Adj= zeros(N);
for va = 1:N
    for vb = 1:N
        for h = 1:length(m)-1
            for g = 1:length(m)-h
                if Qm(g) == va && Qm(g+h) == vb
                    Adj(va,vb) = 1;
                    %%去除自环
                    if (va == vb)
                        Adj(va,vb)=0;
                    end
                    %%%%
                    break
                end
            end
            break  %加break只计算步长h等于1，不加所有步长都计算
        end
    end
end
%不对称的邻接矩阵取上三角再做对称
if rank(Adj-Adj') ~= 0
    upper = triu(Adj);
    Adj = upper + upper';
end


%创建连通图的对象
G=graph(Adj);

%计算总边数
E = sum(sum(Adj));
end