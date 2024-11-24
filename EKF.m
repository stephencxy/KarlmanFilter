%EKF算法

dt = 0.01;
t = 0:dt:1;
L = length(t);

x = zeros(1,L);
y = zeros(1,L);
x(1) = 0.1;
y(1) = 0.1^2;
%生成真实信号和观测信号
% x_k = sin(3*x_k-1);
% y_k = x_k^2;
for i = 2:L
    x(i) = sin(3*x(i-1));
    y(i) = x(i)^3+normrnd(0,0.01);
end

%初始化
Xplus = zeros(1,L);
Pplus = 0.1;
Xplus(1) = 0.1;
Q = 0.001;
R = 10;
for i = 2:L
    %预测步
    A = 3*cos(3*Xplus(i-1));%线性化 x_k=f(x_k-1)+Q; A=f`(x_k-1);
    Xminus = sin(3*Xplus(i-1));
    Pminus = A*Pplus*A.'+Q;
    %更新步
    C = 3*Xminus;%线性化 y_k=h(x_k);C=h`(x_k);
    K = Pminus*C*inv(C*Pminus*C.'+R);
    Xplus(i) = Xminus+K*(y(i)-Xminus^2);
    Pplus = (eye(1)-K*C)*Pminus;
end

plot(t,x,t,Xplus)