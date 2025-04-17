
function [G,largest_eig,G_sig,Gwao]=signal2graph_PQChu1(r,N,flag)
%量化值与顶点概率向量融合
%r输入时域样本
%r图的顶点数
%flag 输入信号的形式，0-表示自相关函数 1-表示功率谱
k = length(r);
switch flag
    case 0 
        R = abs(fft(r)).^2/length(r); % 功率谱
%         gl = 2;
% %         R_new = block(R,gl,'max')+block(R,gl,'min');
% % %         
%          group = zeros(gl,ceil(length(R)/gl));
%          for i = 1:length(R)
%               group(i) = R(i);
%          end
%          R_new = zeros(1,ceil(length(R)/gl));
%          for j = 1:length(group)
%              R_new(j) = sum(group(:,j));
%          end
% %            
%           R = R_new;  
          m = 0:length(R)-1;
    case 1
        %对信号做自相关
        
        r_mean = mean(r);%求s均值
        rnew =r-r_mean;%去掉均值
        [R,m] = xcorr(rnew,k-1);
        R = R(end-k+1:end);%取对称的自相关函数的一半
        m = m(end-k+1:end);
end
%归一化
       Rmax = max(R);
       Rmin = min(R);
       U = (R-Rmin)/(Rmax-Rmin);
%非线性变换
%        delta  = 1;
%        U = delta*sinh(log(2)*2*V/delta);
% % 再归一化
%        Umax = max(U); Umin = min(U);
%        U = (U-Umin)/(Umax-Umin);
   %%%%%%   
    
      
      
      
      
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
%%%%%%%%%%%%%%%%%%%%%
%顶点概率向量
% h = histogram(Qm,N,'Normalization','probability');
[h,bin] = histcounts(Qm,N,'Normalization','probability');
Vp = h;  % 顶点概率向量
Vp = sort(Vp,'descend');
%输出GWAO统计量
QQm = [1:1:N];
Gwao = sum(Vp.*QQm);
% %转换成均匀图的邻接矩阵
Adj= zeros(N);
for h = 1:1:length(Qm)-1
    for g = 1:1:length(m)-h
        if Qm(g)~=Qm(g+h)%步长取所有
%              Adj(Qm(g),Qm(g+h))=1;
%              Adj(Qm(g+h),Qm(g))=1;
              Adj(Qm(g),Qm(g+h)) =  Adj(Qm(g),Qm(g+h))+1;
              Adj(Qm(g+h),Qm(g)) =  Adj(Qm(g+h),Qm(g))+1;
           
        end
    end
%      break; %单步长
end
% %%%%
% % %创建连通图的对象
  G = graph(Adj);
%  Lap = laplacian(G);
% %度矩阵
% d = degree(G);
% D = zeros(N);
% for deg = 1:N
%     D(deg,deg) = d(deg);
% end
% d1 = d/sum(d);
% %laplacian矩阵
% Lap = D-Adj;
% Laps = D+Adj;
% if rank(Lap-Lap')~=0
%     disp('存在不对称矩阵');
% end
% 
% %laplacian矩阵的特征值
% % eigen = eig(Lap);
% % Lap = Adj;
%  [eVector,eigen] = eigs(Lap);
% % [eVector1,eigen1] = eig(Laps);
% % % IEL = sum(sqrt(eign)); % 拟拉普拉斯能量
% % % EI = sum(sqrt(eigen1)); %关联能量
%  sort_eig = sort(eigen,'descend');
% % %最大拉普拉斯特征值与顶点和边数的关系
  largest_eig = 1;
% compare = N;
% %%%% 自环图
%      %%%%%%%带环图
%       Adj1 = zeros(N);
%         for h = 1:1:length(Qm)-1
%             for g = 1:1:length(m)-h
%                 if Qm(g)~=Qm(g+h)
% %                     Adj(Qm(g),Qm(g+h))=Adj(Qm(g),Qm(g+h))+1;
%                      Adj1(Qm(g),Qm(g+h))= 0;
% %                     Adj(Qm(g+h),Qm(g))=Adj(Qm(g+h),Qm(g))+1;
%                      Adj1(Qm(g+h),Qm(g))= 0;
%                 elseif Qm(g)==Qm(g+h)
%                      Adj1(Qm(g),Qm(g+h))= Adj1(Qm(g),Qm(g+h))+1;
% %                     Adj(Qm(g),Qm(g+h))= 1;
%                 end
%             end
% %             break  %单
%         end
% %%%%%%
% Adj1;
%  G_sig = sum(Adj1)/sum(sum((Adj1)));
%  G_sig0 = sum(Adj1);
% %  Gwhu = G_sig0;
% %%%%%%
% %输出Hu统计量
%  eigv = diag(eigen)';
%  eigv1 = diag(eigen1)';
%  Gwhu = sum(eigv.*G_sig);
%  Gwhu = sum(d1'.*G_sig);
%  +sum(eigv1.*G_sig);
%  Gwhu = sum(Vp.*G_sig); %顶点概率向量与自环 0.8940,0.8738(BM)
% Gwhu = sum(Vp.*eigv1); %顶点概率向量与特征值 0.9121

%  GFT = (inv(eVector)*G_sig')';
% Gwhu = GFT(2);
%  G_sig = Vp.*QQm; %性能不错0.9129
%  G_sig = Vp.*G_sig; %与自环方法差不多
  G_sig = QQm;%不行
%   Gwhu = sum(abs(inv(eVector1)*G_sig').^2); %%% GFT的能量 0.8920
% *Gwhu ;
% + sum(abs(inv(eVector1)*G_sig').^2);
% ;
% + 
%  Gwhu = sum(eigv.*QQm);
%%%%%%%

end