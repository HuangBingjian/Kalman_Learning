% SIR粒子滤波的应用，详见http://blog.csdn.net/heyijia0327/article/details/40899819

clear all 
close all 
clc

% 初始化变量
x = 0.1;    % 初始真实状态
x_N = 1;    % 系统过程噪声的协方差(一维对应的是方差)
x_R = 1;    % 测量的协方差
T = 75;     % 共进行75次
N = 100;    % 粒子数，粒子越多效果越好，计算量也越大

V = 2;      % 初始分布的方差
x_P = [];   % 粒子

% 高斯分布随机产生初始粒子
for i = 1:N
    x_P(i) = x+ sqrt(V) * randn;
end

z_out = [x^2/20 + sqrt(x_R) * randn];   % 实际测量值
x_out = [x];                            % 实际值
x_est = [x];                            % 预测值
x_est_out = [x_est];                    % 

for t = 1:T
    x = 0.5*x + 25*x/(1+x^2) + 8*cos(1.2*(t-1)) + sqrt(x_N)*randn;
    z = x^2/20 + sqrt(x_R)*randn;
    for i = 1:N
        %从先验p(x(k)|x(k-1))中采样  
        x_P_update(i) = 0.5*x_P(i) + 25*x_P(i)/(1 + x_P(i)^2) + 8*cos(1.2*(t-1)) + sqrt(x_N)*randn;  
        %计算采样粒子的值，为后面根据似然去计算权重做铺垫  
        z_update(i) = x_P_update(i)^2/20;  
        %对每个粒子计算其权重，这里假设量测噪声是高斯分布。所以 w = p(y|x)对应下面的计算公式  
        P_w(i) = (1/sqrt(2*pi*x_R)) * exp(-(z - z_update(i))^2/(2*x_R));  
    end  
    % 归一化
    P_w = P_w./sum(P_w);  
        
    % Resampling没有用histc函数，不过目的和效果是一样的  
    for i = 1 : N  
        x_P(i) = x_P_update(find(rand <= cumsum(P_w),1));   % 粒子权重大的将多得到后代  
    end                                                     % find( ,1) 返回第一个符合前面条件的数的下标  
      
    % 状态估计，重采样以后，每个粒子的权重都变成了1/N  
    x_est = mean(x_P);  
      
    % Save data in arrays for later plotting  
    x_out = [x_out x];  
    z_out = [z_out z];  
    x_est_out = [x_est_out x_est];  
      
end  
  
t = 0:T;  
figure(1);  
clf  
plot(t, x_out, '.-b', t, z_out, '-.g', t, x_est_out, '-.r','linewidth',3);  
set(gca,'FontSize',12); set(gcf,'Color','White');  
xlabel('time step'); ylabel('flight position');  
legend('True flight position', 'Measurement', 'Particle filter estimate');  
