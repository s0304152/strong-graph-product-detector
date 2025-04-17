function [noisy_signal, noise] = add_laplacian_noise(s, snr_dB, signal_power)
% ADD_LAPLACIAN_NOISE 添加拉普拉斯噪声至信号，并根据SNR控制噪声功率
%   [noisy_signal, noise] = add_laplacian_noise(s, snr_dB) 
%   输入:
%       s: 原始信号（实数或复数）
%       snr_dB: 信噪比（dB）
%       signal_power (可选): 信号功率（若未提供，自动计算）
%   输出:
%       noisy_signal: 加噪后的信号
%       noise: 生成的拉普拉斯噪声

    % 计算信号功率
    if nargin < 3 || isempty(signal_power)
        signal_power = mean(abs(s).^2); % 默认计算信号平均功率
    end
    
    % 将SNR从dB转为线性值，并计算噪声功率
    snr_linear = 10^(snr_dB/10);
    noise_power = signal_power / snr_linear;
    
    % 生成拉普拉斯噪声（参数b由噪声功率决定）
    b = sqrt(noise_power / 2); % 拉普拉斯噪声方差为2b²
    
    % 生成标准均匀分布随机数
    u = rand(size(s)) - 0.5; 
    
    % 通过逆变换法生成拉普拉斯噪声
    noise = -b * sign(u) .* log(1 - 2*abs(u));
    
    % 可选：强制调整噪声功率（确保数值稳定性）
    current_noise_power = mean(noise.^2);
    if current_noise_power ~= 0
        noise = noise * sqrt(noise_power / current_noise_power);
    end
    
    % 添加噪声至信号
    noisy_signal = s + noise;
end