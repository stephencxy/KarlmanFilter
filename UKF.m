%UKF算法
%命名 plus->+; minus->-; hat->^;
clear
clc

dt = 0.01;
t = 0:dt:1;
L = length(t); 

x = zeros(2,L);
z = zeros(2,L);
x(1,1) = 0.2;
x(2,1) = 0.2;

for i = 2:L
    x(1,i) = sin(x(1,i-1))+2*cos(x(2,i-1));%真实数据
    x(2,i) = 3*cos(x(1,i-1))+sin(x(2,i-1));
    z(1,i) = x(1,i)+x(2,i)+normrnd(0,1);%观测数据
    z(2,i) = x(1,i)+x(2,i)+normrnd(0,1);
end
%初始化
X = zeros(2,L);
X(1,1) = 0.1;
X(2,1) = 0.2;
Pplus = eye(2);
Q = 1e-2*eye(2);
R = 1*eye(2);
%设置权重
n = 2;%X的维度
w = zeros(1,2*n+1);
lamda = 1;

for i = 1:2*n+1
    w(i) = 1/(2*(n+lamda));
end
w(1) = lamda/(n+lamda);

%UKF
for i = 2:L
    xsigma = zeros(n,2*n+1);
    UL = chol(Pplus);%矩阵分解
    xsigma(:,1) = X(:,i-1);
    for j = 1:n
        xsigma(:,j+1) = xsigma(:,1)+sqrt(n+lamda)*UL(:,j);
        xsigma(:,j+1+n) = xsigma(:,1)-sqrt(n+lamda)*UL(:,j);
    end
    %预测步
    xsigmaminus = zeros(n,2*n+1);
    for j = 1:2*n+1
        xsigmaminus(1,j) = sin(xsigma(1,j))+2*cos(xsigma(2,j));
        xsigmaminus(2,j) = 3*cos(xsigma(1,j))+sin(xsigma(2,j));
    end
    %求期望方差
    xhatminus = zeros(n,1);
    Pminus = zeros(n,n);
    for j = 1:2*n+1
        xhatminus = xhatminus + w(j)*xsigmaminus(:,j);
    end
    for j = 1:2*n+1
        Pminus = Pminus+w(j)*(xsigmaminus(:,j)-xhatminus)*(xsigmaminus(:,j)-xhatminus).';
    end 
    Pminus = Pminus +Q;
    %更新步
    xsigma = zeros(n,2*n+1);
    xsigma(:,1) = xhatminus;
    UL1 = chol(Pminus);
    for j = 1:n
        xsigma(:,j+1) = xsigma(:,1)+sqrt(n+lamda)*UL1(:,j);
        xsigma(:,j+1+n) = xsigma(:,1)-sqrt(n+lamda)*UL1(:,j);
    end
    yhat = zeros(n,1);
    y = zeros(n,1);
    for j = 1:2*n+1
        y(1,j) = xsigma(1,j)+xsigma(2,j);
        y(2,j) = xsigma(1,j)+xsigma(2,j);
        yhat = yhat+w(j)*y(:,j);
    end
    Py = zeros(n,n);
    Pxy = zeros(n,n);
    for j = 1:2*n+1
        Pxy = Pxy+w(j)*(xsigma(:,j)-xhatminus)*(y(:,j)-yhat).';
        Py = Py+w(j)*(y(:,j)-yhat)*(y(:,j)-yhat).';
    end
    Py = Py+R;
    %卡尔曼增益
    K = Pxy*inv(Py);
    %观测数据
    Y = zeros(n,1);
    Y(1,1) = z(1,i);
    Y(2,1) = z(2,i);
    
    X(:,i) = xhatminus+K*(Y-yhat);
    
    Pplus = Pminus+K*Py*K.';
end

plot(t,x(1,:),t,X(1,:))