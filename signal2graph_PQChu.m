
function [G,Lap,G_sig,Gwhu,Adj,GFT]=signal2graph_PQChu(r,N,flag)
% [G,Lap,G_sig,Gwhu,Adj,GFT]=signal2graph_PQChu(r,N,flag)
%GFT L  Vp
%r输入时域样本
%r图的顶点数
%flag 输入信号的形式，1-表示自相关函数 0-表示功率谱
k = length(r); 
gl = 2;
switch flag
    case 0 
        R = abs(fft(r)).^2/length(r); % 功率谱
        
         R_new = block(R,gl,'max')+block(R,gl,'min');
%          plot(R);
%          xlabel('The discrete frequency,k'); % 横坐标标签
%          ylabel('Power Spectrum');    % 纵坐标标
% %         
%          group = zeros(gl,ceil(length(R)/gl));
%          for i = 1:length(R)
%               group(i) = R(i);
%          end
%          R_new = zeros(1,ceil(length(R)/gl));
%          for j = 1:length(group)
%              R_new(j) = sum(group(:,j))/gl;
%          end
% %            
          R = R_new;  
          m = 0:length(R)-1;
    case 1
        %对信号做自相关
        
        r_mean = mean(r);%求s均值
        rnew = r-r_mean;%去掉均值
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
%%%%%%%%%%%%%%%%%%%%%定义新信号
% for kk1 = 1:N
%     temp_index = find (Qm ==kk1);
%     if length(temp_index)==0;
%         temp_index = kk1;
%     end
%     Gav(kk1) = max(U(temp_index))+min(U(temp_index));
% end 
%  Gav;



%%%%%%%%%%%%%%%%%%%%%
%顶点概率向量
% h = histogram(Qm,N,'Normalization','probability');
[h,bin] = histcounts(Qm,N,'Normalization','probability');
% Qm;
Vp = h;
% % subplot(2,1,1);
% % stem(Vp, 'LineWidth', 1.5, 'Color', [0 0.447 0.741]); % 设置线宽和颜色
% % xlabel('The Number of Vertices '); % 横坐标标签
% % ylabel('Vertex Probability)');    % 纵坐标标签
% % title('Vertex Probability Vector Distribution'); % 图形标题
% % grid on; % 显示网格线
% 
% % 可选：添加坐标轴范围限制（根据实际数据调整）
% xlim([1 length(Vp)]);
% ylim([0 max(Vp)*1.1]);
% Vp = h.Values; % 顶点概率向量
%  Vp = sort(Vp,'descend');
 
%  add = rand(length(Qm),1);
% Qm1 = Qm+add;
% [h1,bin1] = histcounts(Qm1,N,'Normalization','probability');
% Vpn = h1; % 顶点概率向量
% Vpn = sort(Vpn,'descend');
% Vps = cumsum(Vp);
%输出GWAO统计量
% QQm = [1:1:N];
% Gwao = sum(Vp.*QQm);
% Gigt = Vp.*QQm;
%转换成均匀图的邻接矩阵
Adj= zeros(N);
for h = 1:1:length(Qm)-1
    for g = 1:1:length(m)-h
        if Qm(g)~=Qm(g+h)%步长取所有
            %对于正弦干扰
%                 Adj(Qm(g),Qm(g+h)) =  Adj(Qm(g),Qm(g+h))+1;
%                 Adj(Qm(g+h),Qm(g)) =  Adj(Qm(g+h),Qm(g))+1;
                 Adj(Qm(g),Qm(g+h))=1;
                 Adj(Qm(g+h),Qm(g))=1;
        elseif Qm(g)==Qm(g+h)
                 Adj(Qm(g),Qm(g+h))= 0;
%                   Adj(Qm(g+h),Qm(g))= Adj(Qm(g+h),Qm(g))+1;
        end
    end
%      break; %单步长
end
%%%%
% Adj = Adj -diag(Adj);
%创建连通图的对象

G = graph(Adj);
% subplot(2,1,2);
% plot(G);
% plot(G)
%度矩阵
d = degree(G);
D = zeros(N);
for deg = 1:N
    D(deg,deg) = d(deg);
end
% d1 = d/sum(d);
%laplacian矩阵
Lap = D-Adj;

% Laps = D+Adj;
% if rank(Lap-Lap')~=0
%     disp('存在不对称矩阵');
% end

%  Lap = laplacian(G);
%laplacian矩阵的特征值
%  eigen = eig(Lap);
%   IEL = sum(sum(sqrt(abs(eigen)))); % 拟拉普拉斯能量
% EI = sum(sum(sqrt(eigen1))); %关联能量
% Lap = Adj;
%  Lap = Laps;
%  Lap00 = D^(-1/2)*Lap*D^(1/2)
% [eVector,eigen] = eig(Lap);
% [eVector1,eigen1] = eig(Laps);

% 
% sort_eig = sort((eigen),'descend');
% %最大拉普拉斯特征值与顶点和边数的关系
% largest_eig = sort_eig(1);
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
%                 elseif  Qm(g)==Qm(g+h)
% %                       Adj1(Qm(g),Qm(g+h))= 1;
%                       Adj1(Qm(g),Qm(g+h))=  Adj1(Qm(g),Qm(g+h))+1;
%                 end
%             end
% %              break  %单
%         end
%%%%%%
% Adj1
% Lg = Adj1+Laps; %分组信号 更好
%   Lg = Adj1+Laps;
 Lg = Lap;
%  Lg = Adj;
  [eVector,eigen] = eig(Lg);
%  G_sig = sum(Adj1)/sum(sum((Adj1)));
%  G_sig = sum(Adj1); %gao
% G_sig = sort(G_sig,'descend');
% stest = (sum(Adj1));
%%%%%%
%输出Hu统计量
 eigv = diag(eigen)';
%  eigv = sort((eigv),'descend');
%  
%  eigv1 = diag(eigen1)';
%  IEL = (sum(sqrt(abs(eigv))));
%   Gwhu = sum(eigv1.*G_sig);  % 不重排，分组和
%  Gwhu = sum(d1'.*G_sig);
%  +sum(eigv1.*G_sig);
%  Gwhu = sum(Vp.*G_sig);
%  G_sig = Gav;
% Gwhu = sum(eigv1.*G_sig);%顶点概率向量与自环 0.8940,0.8738(BM) 比较好 不排 单重图
% eigv1 = eigv1./max(eigv1);
% eigv = eigv./max(eigv);
%    Gwhu = sum(Vp.*eigv); %顶点概率向量与特征值 0.9121 分组求和时了高 多重权重科
%     G_sig = Vp;
    G_sig = Vp;
%     G_sig = (Gav-mean(Gav))./std(Vp); %GFT要用拉普拉斯阵来进行GFT变换。
D = inv(eVector);
% D = D^0.8;
  GFT = (D*G_sig')';
% Gwhu = GFT(2);
%
% %  .*QQm;
%  .*QQm;
%  Gwhu = sum(abs(inv(eVector)*G_sig').^2); %%% GFT的能量 0.8920 BS
% * ;
%       Gwhu = sum(abs(inv(eVector)*G_sig').^2);%可以BM 94
       Gwhu = abs(sum(diag(eigen).*GFT.^2')); %简单图，多重图不行
%    Lg = Lap;
% Lg = inv(D)*Lap;
%           Gwhu = G_sig*Lg*G_sig';  %简单图不能用多重图。
% %         Gwhu = G_sig*Adj*G_sig';
%   G_sig = Gav;
%     Gwhu = G_sig*Lg*G_sig';
% % %   Gwhu = IEL;
% %%%%%%%
% [eVectorA,eigenA] = eig(Adj);
% sort_eigA = sort(diag(eigenA),'descend');
% AdjN = Adj/sort_eigA(1);
% TV = sum(abs(G_sig'-AdjN*G_sig')); %加了绝对值反不好
% TV = 0;
end