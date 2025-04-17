
function [Gout,Gwhu,Gspe,Gac,G_sig,G_sig1,G_sig2]=signal2graph_Product_sig(r,N,gl,flag,flag_sig)
% function [Gout,Gspe,Gac,Gwhu,G_sig,G_sig1,G_sig2]=signal2graph_Product(r,N,gl,flag)
% [G,Lap,G_sig,Gwhu,Adj,GFT]=signal2graph_PQChu(r,N,flag)
%将L阵的特征值与自环值进行线性组合
%r输入时域样本
%r图的顶点数
%flag 图积的方式flag = 0：笛卡尔图积
% flag = 1：克罗内克图积
% flag = 2：强图积
% flag = 3：字典积
k = length(r);
 % 功率谱及图转换
        R = abs(fft(r)).^2/length(r); % 功率谱
%         gl = 2;
%          R_new = block(R,gl,'max')+block(R,gl,'min');
% %         
         group = zeros(gl,ceil(length(R)/gl));
         for i = 1:length(R)
              group(i) = R(i);
         end
         R_new = zeros(1,ceil(length(R)/gl));
         for j = 1:length(group)
             R_new(j) = sum(group(:,j))/gl;
         end
%            
          R = R_new;  
          m = 0:length(R)-1;
         Rmax = max(R);
       Rmin = min(R);
       U = (R-Rmin)/(Rmax-Rmin);
        %%% 均匀量化    
Qm=m;
for c = 1:length(m)
    for d = 0:N-1
        if U(c)>=d/N && U(c)<(d+1)/N
            Qm(c) = d + 1;
            break
        else
            Qm(c) = N;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%定义新信号
for kk1 = 1:N
    temp_index = find (Qm == kk1);
    if length(temp_index)==0;
        temp_index = kk1;
    end
    Gav(kk1) = sum(U(temp_index));
end 
%   Gav = Gav/sum(Gav);
%%%%%%%%%%%%%%%%%%%%%
%顶点概率向量
 add = rand(length(Qm),1);
 Qm1 = Qm+add;

% h = histogram(Qm,N,'Normalization','probability');
[h,bin] = histcounts(Qm,N,'Normalization','probability');
 Vp = h; % 顶点概率向量
% Vp = sort(Vp,'descend');
% %输出GWAO统计量
% QQm = [1:1:N];
% Gwao = sum(Vp.*QQm);
% Gigt = Vp.*QQm;
%转换成有向图
% %%% 均匀量化    
% Qm=m;
% for c = 1:length(m)
%     for d = 0:N-1
%         if U(c)>=d/N && U(c)<(d+1)/N
%             Qm(c) = d + 1;
%             break
%         else
%             Qm(c) = N;
%         end
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%定义新信号
% for kk1 = 1:N
%     temp_index = find (Qm ==kk1);
%     if length(temp_index)==0;
%         temp_index = kk1;
%     end
%     Gav(kk1) = sum(U(temp_index));
% end 
% 
% %%%%%%%%%%%%%%%%%%%%%
% %顶点概率向量
%  add = rand(length(Qm),1);
%  Qm1 = Qm+add;
% 
% % h = histogram(Qm,N,'Normalization','probability');
% [h,bin] = histcounts(Qm,N,'Normalization','probability');
%  Vp = h; % 顶点概率向量
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Adj= zeros(N);
for h = 1:1:length(Qm)-1
    for g = 1:1:length(Qm)-h
        if Qm(g)~=Qm(g+h)%步长取所有
            Adj(Qm(g),Qm(g+h))=1;
            Adj(Qm(g+h),Qm(g))=1;
        end
    end
end
%创建连通图的对象
%%%%%%%带环图
%       Adjs = zeros(N);
%         for h = 1:1:length(Qm)-1
%             for g = 1:1:length(m)-h
%                 if Qm(g)~=Qm(g+h)
% %                     Adj(Qm(g),Qm(g+h))=Adj(Qm(g),Qm(g+h))+1;
%                      Adjs(Qm(g),Qm(g+h))= 0;
% %                     Adj(Qm(g+h),Qm(g))=Adj(Qm(g+h),Qm(g))+1;
% %                      Adj1(Qm(g+h),Qm(g))= 0;
%                 elseif  Qm(g)==Qm(g+h)
%                        Adjs(Qm(g),Qm(g+h))= 1;
% %                       Adj1(Qm(g),Qm(g+h))=  Adj1(Qm(g),Qm(g+h))+1;
%                 end
%             end
% %              break  %单
%         end
%%%%%%
a = 0;
Adj1 = Adj;
 G1 = graph(Adj1);
  Gspe = G1 ;
%  figure (1);
%  plot(G1);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  G_sig1  = sort(Gav,'descend');
%  G_sig1_q  = Gav; %quantization based 
%  G_sig1_v = Vp;  % VPV
 
G_sig1_q  = sort(Gav,'descend');
G_sig1_v = sort(Vp,'descend');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %自相关图变换   
          
       
  %对信号做自相关
        
        r_mean = mean(r);%求s均值
        rnew = r-r_mean;%去掉均值
        [R,m] = xcorr(rnew,k-1);
        R = R(end-k+1:end);%取对称的自相关函数的一半
        m = m(end-k+1:end);

%归一化
       Rmax = max(R);
       Rmin = min(R);
       U = (R-Rmin)/(Rmax-Rmin); %非线性时%掉
%         V = (R-Rmin)/(Rmax-Rmin);
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%非线性变换
%         delta  = 1;
%         U = delta*sinh(log(2)*2*V/delta);
        
% 再归一化
%         Umax = max(U); Umin = min(U);
%         U = (U-Umin)/(Umax-Umin);
  %%%%%%% %%%%%%   
         
%%% 均匀量化    
Qm=m;
for c = 1:length(m)
    for d = 0:N-1
        if U(c)>=d/N && U(c)<(d+1)/N
            Qm(c) = d + 1;
            break
        else
            Qm(c) = N;
        end
    end
end
 %%%%%%%%%%%%%%%%%%%%%定义新信号
for kk1 = 1:N
    temp_index = find (Qm ==kk1);
    if length(temp_index)==0;
        temp_index = kk1;
    end
    Gav(kk1) = sum(U(temp_index));
end 
%    Gav = Gav/sum(Gav);
%%%%%%%%%%%%%%%%%%%%%
%顶点概率向量
 add = rand(length(Qm),1);
 Qm1 = Qm+add;

% h = histogram(Qm,N,'Normalization','probability');
[h,bin] = histcounts(Qm,N,'Normalization','probability');
 Vp = h; % 顶点概率向量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Adj= zeros(N);
for h = 1:1:length(Qm)-1
    for g = 1:1:length(Qm)-h
        if Qm(g)~=Qm(g+h)%步长取所有
            Adj(Qm(g),Qm(g+h))=1;
            Adj(Qm(g+h),Qm(g))=1;
        end
    end
end

%创建连通图的对象
%%%%%%%带环图
%       Adjs = zeros(N);
%         for h = 1:1:length(Qm)-1
%             for g = 1:1:length(m)-h
%                 if Qm(g)~=Qm(g+h)
% %                     Adj(Qm(g),Qm(g+h))=Adj(Qm(g),Qm(g+h))+1;
%                      Adjs(Qm(g),Qm(g+h))= 0;
% %                     Adj(Qm(g+h),Qm(g))=Adj(Qm(g+h),Qm(g))+1;
% %                      Adj1(Qm(g+h),Qm(g))= 0;
%                 elseif  Qm(g)==Qm(g+h)
%                        Adjs(Qm(g),Qm(g+h))= 1;
% %                       Adj1(Qm(g),Qm(g+h))=  Adj1(Qm(g),Qm(g+h))+1;
%                 end
%             end
% %              break  %单
%         end
%%%%%%
% a = 0;

Adj2 = Adj;
%创建连通图的对象
 G2 = graph(Adj2);
 Gac = G2;
%  figure (2);
%  plot(G2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G_sig2  = sort(Gav,'descend');
% G_sig2_v = sort(Vp,'descend');
% G_sig2_q  = Gav;
% G_sig2_v = Vp;

G_sig2_q  = sort(Gav,'descend');
G_sig2_v = sort(Vp,'descend');

%  G_sig2_v = Vp.*(1-Vp)*10;
%  G_sig2_v = log(Vp);
% G_sig2_v = G_sig2_v./sum(G_sig2_v) ;
% flag = 3;% cardinate product
coef = 0.6;
%% 
%%%% 求补图，对
%  Adj2 = ones(N)-Adj2-eye(N); %自相关图。
% G2 = graph(Adj2);

%   figure (2);
%   plot(G2);

%%%%%%%%%%%%%%
[Adj_P,L,Apa] = Graph_product(Adj1, Adj2, flag,coef);
%%%%%%%  计算统计量
%          G_sig = reshape(kron(G_sig1,G_sig2),[N^2,1]);
% flag_sig = 0; %configuration of the graph signal
switch flag_sig
    case 0
        G_sig1 = G_sig1_v;
        G_sig2 = G_sig2_q;
    case 1
        G_sig1 = G_sig1_q;
        G_sig2 = G_sig2_v;
    case 2
        G_sig1 = G_sig1_v;
        G_sig2 = G_sig2_v;
    case 3
        G_sig1 = G_sig1_q;
        G_sig2 = G_sig2_q;
end
% case 0
%         G_sig1 =  boxcox(G_sig1_v);
%         boxcox
%         G_sig2 =  boxcox(G_sig2_q);
%     case 1
%         G_sig1 =  boxcox(G_sig1_q);
%         G_sig2 =  boxcox(G_sig2_v);
%     case 2
%         G_sig1 =  boxcox(G_sig1_v);
%         G_sig2 =  boxcox(G_sig2_v);
%     case 3
%         G_sig1 =  boxcox(G_sig1_q);
%         G_sig2 =  boxcox(G_sig2_q);
% end

    G_sig = reshape(kron(G_sig1,G_sig2),[N^2,1]);
    % the graph signal on product graph
    G_sig = sort(G_sig,'descend');
    G_sig = G_sig;
%     ./sum(G_sig);
     
      %自相关非线性时好，线性时也好 ，对于APa
%           G_sig = reshape(kron(G_sig1,G_sig2_v),[N^2,1]);
%     G_sig = reshape(kron(G_sig1_v,G_sig2_v),[N^2,1]); %均匀量化
%       L = Apa;  %原来没有用，后来2/11加的
%      L = Adj_P; 
% G_sig = sort(G_sig,'descend');
% 
   Gwhu = G_sig'*L*G_sig;% 老的/
%      Gwhu = max(L*G_sig);  %用单跳滤波器，0209

%  Gwhu = (max(abs(G_sig'*L).^1)).^(1/1);
%  Gwhu = (sum(abs(G_sig'*L).^1)).^(1/1);
 Gout = graph(Adj_P);
%  G_sig1 = G_sig1_v;
%  G_sig2  =  G_sig2;
%  figure (3);
%  plot(Gout);
%%%%%%
end