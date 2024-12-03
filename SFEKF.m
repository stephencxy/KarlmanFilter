clear
clc
%SFEKF算法 min E(eta_j+1 eta_j)

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
    y(i) = x(i)^2+normrnd(0,0.01);
end
d=zeros(1,L);
%初始化
Xplus = zeros(1,L);
Pplus = 0.1;
Xplus(1) = 0.1;
Q = 0.01;
R = 10;
beta = 1;%弱化因子
lamda = 1;%渐消因子
fhi = 2;%迭代步长
for i = 2:L
    %预测步
    A = 3*cos(3*Xplus(i-1));%线性化 x_k=f(x_k-1)+Q; A=f`(x_k-1);
    Xminus = sin(3*Xplus(i-1));
    C = 3*Xminus^2;%线性化 y_k=h(x_k);C=h`(x_k);
    %残差协方差矩阵
    if i == 2
        V = (y(i)-Xminus^2)*(y(i)-Xminus^2).';
    else
        V = (0.95*V+(y(i)-Xminus^2)*(y(i)-Xminus^2).')/(1+0.95);
    end
    a(i) = V;
    %次优渐消因子
    N = V-C*Q*C.'-beta*R;
    M = C*A*Pplus*A'*C';
    lamda0 = trace(N)/trace(M);
    if lamda0 >= 1
        lamda = lamda0;
    end
    Pminus = lamda*A*Pplus*A.'+Q;
    %更新步
    K = Pminus*C*inv(C*Pminus*C.'+R);
    Xplus(i) = Xminus+K*(y(i)-Xminus^2);
    Pplus = (eye(1)-K*C)*Pminus;
    %预测结束，计算残差
%最优渐消因子
%     if i == 2
%         V = (y(i)-Xminus)*(y(i)-Xminus).';
%     else
%         V = (0.95*V+(y(i)-Xminus)*(y(i)-Xminus).')/(1+0.95);
%         W = Pminus*C.'-K*V;
%         [n,m] = size(W);%n矩阵行数，m矩阵列数
%         for k = 1:n
%             for j = 1:m
%                 B = C*Pminus*C.'+R;
%                 diff_g_lamda = A*Pplus*A.'*C.'*(eye(1)-inv(B)*V)+K*C*A*Pplus*A.'*C.'*inv(B)*V;
%                 g_lamda = 2*W(k,j)*diff_g_lamda;
%                 lamda  = lamda-fhi*g_lamda;
%                 d(i)=lamda;
%             end
%         end
%     end
end

plot(t,x,t,Xplus)