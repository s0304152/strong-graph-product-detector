
% 图积函数
% coef:Aa谱的比例
% 输入参数：邻接矩阵Adj1，Adj2; 图积方式flag
% flag = 0：笛卡尔图积
% flag = 1：克罗内克图积
% flag = 2：强图积
% flag = 3：字典积
% 输出参数：图积的邻接矩阵Adj; 图积的拉普拉斯矩阵L

function [Adj_P,L,Afa] = Graph_product(A_G, A_H,flag,coef)
    coef = 0.8; %default
    A1 = A_G;
    A2 = A_H;

    Len1 = length(A1);
    Len2 = length(A2);

    I1 = eye(Len1);
    I2 = eye(Len2);

    % 图积
    switch flag
        % 笛卡尔图积
        case 0
            Adj_P = kron(A1,I2) + kron(I1,A2);

            % 克罗内克图积
        case 1
            Adj_P = kron(A1,A2);

            % 强图积
        case 2
            Adj_P = kron(A1,A2) + kron(A1,I2) + kron(I1,A2);
        case 3%字典积
            n_G = size(A_G, 1);
            n_H = size(A_H, 1);
            Adj_P = zeros(n_G * n_H);  % 初始化字典序乘积的邻接矩阵
    
           % 计算字典序乘积的邻接矩阵
            for i_G = 1:n_G
                for i_H = 1:n_H
                      for j_G = 1:n_G
                           for j_H = 1:n_H
                                u = (i_G-1)*n_H + i_H;
                                v = (j_G-1)*n_H + j_H;
                    
                    % 如果两个节点来自不同的图，则根据边的情况进行赋值
                               if i_G ~= j_G
                                  Adj_P(u, v) = A_G(i_G, j_G);
                               elseif i_H ~= j_H
                                  Adj_P(u, v) = A_H(i_H, j_H);
                               end
                           end
                      end
               end
           end
    end
    D = diag(sum(Adj_P));
    L = D - Adj_P;
%       L = D^(-1/2)*L*D^(-1/2);
%    Afa = D + Adj_P;
     Afa = coef*D+Adj_P*(1-coef);

end