%用于生成经过MIMO信道的接收信号
%输入参数
%s为MIMO系统的发射信号
%Mt,Mr分别为发射天线数和接收天线数
%Q=Mr是接收信号的过采样因子
%RB为发射信号的码元速率
%f= RB*Q为一个码元周期的采样点数
%k为基带信号的码元个数
%A为载波振幅
%L是信道长度为冲击响应延续时间Th与码元周期TB的商的向上取整，L=Th/TB
function [r0,r1,sigma]=MIMO_new2(s,Mt,Mr,L,snr,flag)
%第i个叠加序列右移i位
% 20221017赵修
S=s;%为生成序列矩阵做准备
H=[];%为生成增益矩阵做准备
H=[H,1];%未移位序列增益
for i=1:L
    eval(['s',num2str(i),'= circshift(s,i);']);%生成变量si并将当前移位后序列存入
    s_Rshift = ['s',num2str(i)];%调用si
    s_Rshift = eval(s_Rshift);%把s_Rshift里存储的字符串名si变为变量si里存储的值
    s_Rshift(1,1:i) = 0;
    S=[S;s_Rshift];%生成序列矩阵
    h=(exp(-4*i/L)).^0.5;%生成对应每个序列的增益
    H=[H,h];%生成增益序列
end
H = repmat(H,Mt,1);%生成增益矩阵
Rmt1 = H*S;
Rmt = awgn(Rmt1,snr,'measured');%添加信道噪声
noise = Rmt - Rmt1;%单独求信道噪声，新增
rsingle = sum(Rmt);%单个接收天线接收值
noise_single = sum(noise);%单个天线噪声，新增
R = repmat(rsingle,Mr,1);
N = repmat(noise_single,Mr,1);%噪声，新增
%flag对应天线不同的接收端叠加方式
%0代表天线相加，1代表天线移位相乘，其他值代表天线移位延续
if flag == 0
    r1 = sum(R);
    r0 = sum(N);%噪声，新增
elseif flag == 1
    r1 = R(1,:).*conj(circshift(R(2,:),1));
    r0 = N(1,:).*conj(circshift(N(2,:),1));%噪声，新增
    % sp = fft(r1);
    % figure(1)
    % plot(abs(sp));
else 
    [w1,w2]=size(R);
    r1 = reshape(R,1,w1*w2);%将扩充后的矩阵整形成1行作为对每一个码元做Q个采样
    r0 = reshape(N,1,w1*w2);
end
sigma = std(r0);
end

