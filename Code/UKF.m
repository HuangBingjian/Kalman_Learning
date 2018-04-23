% ----------------------------------------
% 该文件用于编写无迹卡尔曼滤波算法及其测试
% 注解：主要子程序包括：轨迹发生器、系统方程
%       测量方程、UKF滤波器
% 来源：网上资源
% ---------------------------------------

clc;
tic;    % 用于计算程序运行时间
global Qf n;                  %定义全局变量
% 初始化
stater0=[220; 1;55;-0.5];     %标准系统初值
state0=[200;1.3;50;-0.3];     %测量状态初值
p=[0.005 0 0 0;0 0.005 0 0;
   0 0 0.005 0;0 0 0 0.005];  %状态误差协方差初值                             
n=4; 
T=3;
Qf=[T^2/2 0;0 T;T^2/2 0;0 T];

stater=stater0;state=state0; xc=state;
staterout=[]; stateout=[];xcout=[];
errorout=[];tout=[];
t0=1; h=1; tf=1000;          %仿真时间设置

% 滤波算法
for t=t0:h:tf
    [state,stater,yc]=track(state,stater); %轨迹发生器：标准轨迹和输出
    [xc,p]=UKFfiter(@systemfun,@measurefun,xc,yc,p);
    error=xc-stater;              %滤波处理后的误差
    staterout=[staterout,stater];
    stateout=[stateout,state];
    errorout=[errorout,error];
    xcout=[xcout,xc];  
    tout=[tout,t];
end 

% 显示结果
figure;
plot(tout,xcout(1,:),'r',tout,staterout(1,:),'g',...
     tout,stateout(1,:),'black');
legend('滤波后','真实值','无滤波');
grid on;
xlabel('时间 t（s）');
ylabel('系统状态A');
title('无迹卡尔曼滤波');
figure;
plot(tout,xcout(2,:),'r',tout,staterout(2,:),'g',...
     tout,stateout(2,:),'black');
grid on;
legend('滤波后','真实值','无滤波');
xlabel('时间 t（s）');
ylabel('系统状态B');
title('无迹卡尔曼滤波');
figure;
plot(tout,xcout(3,:),'r',tout,staterout(3,:),'g',...
     tout,stateout(3,:),'black');
grid on;
legend('滤波后','真实值','无滤波');
xlabel('时间 t（s）');
ylabel('系统状态C');
title('无迹卡尔曼滤波');
figure;
plot(tout,xcout(4,:),'r',tout,staterout(4,:),'g',...
     tout,stateout(4,:),'black');
grid on;
legend('滤波后','真实值','无滤波');
xlabel('时间 t（s）');
ylabel('系统状态D');
title('无迹卡尔曼滤波');
figure;
plot(tout,errorout(1,:),'r',tout,errorout(2,:),'g',...
     tout,errorout(3,:),'black',tout,errorout(4,:),'b');
grid on;
legend('A','B','C','D');
xlabel('时间 t（s）');
ylabel('滤波后的状态误差');
title('无迹卡尔曼滤波误差');

toc;  %计算仿真程序运行时间

%---------------------------------------------

function [state,stater,yout]=track(state0,stater0)
%轨迹发生函数

T=3;
F=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];
G=[T^2/2 0;0 T;T^2/2 0;0 T];
V=0.005*randn(2,1);
W=0.008*randn(1,1);
state=F*state0+G*V;                     
stater=F*stater0;                       
yout=atan(stater0(3)/stater0(1))+W;     
%用真实值得到测量值，在滤波时结果才会与真实值重合。
end
 

function  state=systemfun(state0)
% 系统方程
    
T=3;
F=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];
state=F*state0;
end


function  yout=measurefun(state0)
% 测量方程
yout=atan(state0(3)/state0(1));
end


function [xc,p]=UKFfiter(systemfun,measurefun,xc0,yc,p0)
%此程序用于描述UKF（无迹kalman滤波）算法

global Qf n;
%----------------参数注解-------------------
% xc0---状态初值(列向量)  yc---系统测量值
% p0----状态误差协方差    n----系统状态量数
% systemfun---系统方程 measurefun--测量方程
%---------------滤波初始化------------------
alp=1;                                    % default, tunable
kap=-1;                                   % default, tunable
beta=2;                                   % default, tunable
lamda=alp^2*(n+kap)-n;                    % scaling factor
nc=n+lamda;                               % scaling factor
Wm=[lamda/nc 0.5/nc+zeros(1,2*n)];        % weights for means
Wc=Wm;
Wc(1)=Wc(1)+(1-alp^2+beta);               % weights for covariance
ns=sqrt(nc);
%-------------------------------------------
  sxk=0;spk=0;syk=0;pyy=0;pxy=0; p=p0;     
    
%--------------构造sigma点-----------------
pk=ns*chol(p); % A=chol(B);meant:A'*A=B;
sigma=xc0;
for k=1:2*n
    if(k<=n)
       sigma=[sigma,xc0+pk(:,k)];
    else
       sigma=[sigma,xc0-pk(:,k-n)];
    end
end
%-------------时间传播方程----------------
for ks=1:2*n+1
    sigma(:,ks)=systemfun(sigma(:,ks));%利用系统方程对状态预测
    sxk=Wm(ks)*sigma(:,ks)+sxk;
end 
%----------完成对Pk的估计
for kp=1:2*n+1
    spk=Wc(kp)*(sigma(:,kp)-sxk)*(sigma(:,kp)-sxk)'+spk;   
end
    spk=spk+Qf*0.005*Qf';
%-----------------------
for kg=1:2*n+1
    gamma(kg)=measurefun(sigma(:,kg));
end
for ky=1:2*n+1
    syk=syk+Wm(ky)*gamma(ky);
end
%--------------测量更新方程--------------
for kpy=1:2*n+1
    pyy=Wc(kpy)*(gamma(kpy)-syk)*(gamma(kpy)-syk)'+pyy;
end
    pyy=pyy+0.008;
for kxy=1:2*n+1
    pxy=Wc(kxy)*(sigma(:,kxy)-sxk)*(gamma(kxy)-syk)'+pxy;
end
   kgs=pxy/pyy;                    %修正系数
   xc=sxk+kgs*(yc-syk);            %测量信息修正状态
   p=spk-kgs*pyy*kgs';             %误差协方差阵更新
end
