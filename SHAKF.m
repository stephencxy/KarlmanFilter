%SHAKF Sage-Husa自适应卡尔曼滤波(对于高维系统变量该算法的稳定性和收敛性不能保证）
clear
clc

dt = 0.01;
t = 0:dt:1;
L = length(t);

x = zeros(1,L);
y = zeros(1,L);

%生成真实信号和观测信号
% x_k = exp(2t);
% y_k = x_k;
for i = 1:L
    x(i) = exp(2*t(i));
    y(i) = x(i)+normrnd(0,0.01);
end

%SHAKF
Fhi = 1+2*dt;%状态转移矩阵
H = 1;%观测转移矩阵
%初始化
Xplus = zeros(1,L);
Pplus = zeros(1,L);
Qplus = zeros(1,L);
Rplus = zeros(1,L);
qplus = zeros(1,L);
rplus = zeros(1,L);
Pplus(1) = 0.01;
Qplus(1) = 0.1;
Rplus(1) = 1;
qplus(1) = 0.01;
rplus(1) = 0.1;
Xplus(1) = 1;
b = 0.4;%遗忘因子

for i = 2:L
    %预测步
    Xminus = Fhi*Xplus(i-1)+qplus(i-1);
    Pminus = Fhi*Pplus(i-1)*Fhi'+Qplus(i-1);
    %更新步
    dk = (1-b)/(1-b^(i+1));
    rplus(i) = (1-dk)*rplus(i-1)+dk*(y(i)-H*Xminus);
    yerr = y(i)-H*Xminus-rplus(i);%观测误差
    Rplus(i) = (1-dk)*Rplus(i-1)+dk*(yerr*yerr'-H*Pminus*H');
    K = Pminus*H'*inv(H*Pminus*H'+Rplus(i));%卡尔曼增益
    Pplus(i) = (eye(1)-K*H)*Pminus;
    Qplus = (1-dk)*Qplus+dk*(K*yerr*yerr'*K'+Pplus(i)-Fhi*Pplus(i-1)*Fhi');
    Xplus(i) = Xminus+K*yerr;
    qplus(i) = (1-dk)*qplus(i-1)+dk*(Xplus(i)-Fhi*Xplus(i-1));
end

plot(t,x,t,Xplus)
