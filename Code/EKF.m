% 扩展卡尔曼滤波(EKF)算法

kx = .01;                   % x方向阻尼系数
ky = .05;                   % y方向阻尼系数
g = 9.8;                    % 重力
t = 10;                     % 仿真时间
Ts = 0.1;                   % 采样周期
len = fix(t/Ts);            % 仿真步数，向0靠拢取整

% 真实轨迹模拟（类似平抛运动，坐标(0,500),水平初始速度50，垂直初始速度0）
dax = 1.5;                  % 系统噪声
day = 1.5;
X = zeros(len,4);
X(1,:) = [0, 50, 500, 0];   %初始化，[x坐标，x速度，y坐标，y速度]
for k = 2:len
    x =  X(k-1,1);
    vx = X(k-1,2);
    y =  X(k-1,3);
    vy = X(k-1,4);
    
    x = x + vx * Ts;                                    % 水平位置
    vx = vx + (-kx * vx^2 + dax * randn(1,1)) * Ts;     % 减速+噪声
    y = y + vy * Ts;                                    % 垂直位置
    vy = vy + (ky * vy^2 - g + day* randn(1,1)) * Ts;   % 加速+噪声
    
    X(k,:) = [x, vx, y, vy];    
end

figure(1),hold off,
plot(X(:,1),X(:,3),'-b'),grid on
figure(2),
plot(X(:,2:2:4))            % 平抛曲线
    
% 测量
mrad = 0.001;
dr = 10;
dafa = 10*mrad;             % 测量噪声
for k = 1:len
    r = sqrt(X(k,1)^2 + X(k,3)^2) + dr*randn(1,1);      % 距离
    a = atan(X(k,1)/X(k,3)) + dafa*randn(1,1);          % 角度
    Z(k,:) = [r, a];
end

figure(1),hold on,          % 测量所得散点图
plot(Z(:,1).*sin(Z(:,2)), Z(:,1).*cos(Z(:,2)),'*')

% EKF滤波
Qk = diag([0; dax; 0; day])^2;
Rk = diag([dr; dafa])^2;
Xk = zeros(4,1);
Pk = 100*eye(4);                        % 初始协方差矩阵
X_est = X;
for k = 1:len
    Ft = JacobianF(X(k,:), kx, ky);     % 一阶偏导数
    Hk = JacobianH(X(k,:));
    fX = fff(X(k,:), kx, ky, g, Ts);
    hfX = hhh(fX);
    [Xk, Pk, Kk] = ekf(eye(4) + Ft * Ts, Qk, fX, Pk, Hk, Rk, Z(k,:)'-hfX);
    X_est(k,:) = Xk';
end

figure(1),
plot(X_est(:,1), X_est(:,3), '+r')
xlabel('X');ylabel('Y');title('ekf simulation');
legend('real', 'measurement', 'ekf estimated');


%--------------------------------------------------------------------
function F = JacobianF(X, kx, ky)        % 系统状态雅可比函数(一阶导数)
    vx = X(2);
    vy = X(4);
    F = zeros(4,4);
    F(1, 1) = 1;
    F(2, 2) = -2 * kx * vx;
    F(3, 4) = 1;
    F(4, 4) = 2 * ky * vy;
end
    
function H = JacobianH(X)               % 量测雅可比函数  
    x = X(1); 
    y = X(3);  
    H = zeros(2,4);  
    r = sqrt(x^2+y^2);  
    H(1,1) = 1/r; 
    H(1,3) = 1/r;  
    xy2 = 1+(x/y)^2;  
    H(2,1) = 1/xy2*1/y; 
    H(2,3) = 1/xy2*x*(-1/y^2);  
end  
    
function fX = fff(X, kx, ky, g, Ts)     % 系统状态非线性函数  
    x = X(1); 
    vx = X(2); 
    y = X(3); 
    vy = X(4);   
    x1 = x + vx*Ts;  
    vx1 = vx + (-kx*vx^2)*Ts;  
    y1 = y + vy*Ts;  
    vy1 = vy + (ky*vy^2-g)*Ts;  
    fX = [x1; vx1; y1; vy1];  
end  
    
function hfX = hhh(fX)                  % 量测非线性函数  
    x = fX(1); 
    y = fX(3);  
    r = sqrt(x^2+y^2);  
    a = atan(x/y);  
    hfX = [r; a];  
end
    
function [Xk, Pk, Kk] = ekf(Phikk_1, Qk, fXk_1, Pk_1, Hk, Rk, Zk_hfX) % ekf 滤波函数  
    Pkk_1 = Phikk_1*Pk_1*Phikk_1' + Qk;   
    Pxz = Pkk_1*Hk';    
    Pzz = Hk*Pxz + Rk;     
    Kk = Pxz*Pzz^-1;  
    Xk = fXk_1 + Kk*Zk_hfX;  
    Pk = Pkk_1 - Kk*Pzz*Kk';  
end 
    
    
    
    
    
