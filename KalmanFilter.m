dt = 0.01;
t = 0:dt:5;
L = length(t);

signal_true = zeros(1,L);%真实轨迹
signal_observe1 = zeros(1,L);%观测信号1
signal_observe2 = zeros(2,L);%观测信号2
Xplus1 = zeros(1,L);%后验信号
Xplus2 = zeros(3,L);%后验信号

for i = 1:L
    signal_true(i) = 1/2*t(i)^2;
    signal_observe1(i) = signal_true(i) + normrnd(0,0.1);
    signal_observe2(i) = signal_true(i) + normrnd(0,0.01);
end

%X(k)=X(k-1)+Q
%Y(k)=X(k)+R
%初始化
F1 = 1;
H1 = [1;1];
Q1 = 0.01;
R1 = diag([0.1,0.01]);
%初始化X0
Xplus1(1) = 0.01;
Pplus1 = 0.01^2;
%卡尔曼滤波算法
%X(k)minus = F*X(k-1)plus
%P(k)minus = F*P(k-1)plus*F'+Q
%K=P(k)minus*H'*inv(H*P(k)minus*H'+R)
%X(k)plus = X(k)minus + K*(y(k)-H*X(k)minus)
%P(k)plus = (1-K*H)*P(k)minus

for i=2:L
    Xminus1 = F1*Xplus1(i-1);
    Pminus1 = F1*Pplus1*F1'+Q1;
    K1 = Pminus1*H1'*inv(H1*Pminus1*H1'+R1);
    Xplus1(i) = Xminus1+K1*(signal_observe1(i)-H1*Xminus1);
    Pplus1 = (1-K1*H1)*Pminus1;
end
X=movmean(signal_observe1,10);

F2 = [1,dt,0.5*dt^2;0,1,dt;0,0,1];
H2 = [1,0,0;1,0,0];
Q2 = [1,0,0;0,0.01,0;0,0,0.0001];
R2 = diag([20,50]);
Xplus2(:,1) = [0.01;0;0];
Pplus2 = diag([0.01,0.01,0.0001]);
for i=2:L
    Xminus2 = F2*Xplus2(:,i-1);
    Pminus2 = F2*Pplus2*F2.'+Q2;
    K2 = Pminus2*H2.'*inv(H2*Pminus2*H2.'+R2);
    Xplus2(:,i) = Xminus2+K2*(signal_observe1(i)-H2*Xminus2);
    Pplus2 = (eye(3)-K2*H2)*Pminus2;
end

 plot(t,signal_true,t,Xplus1,t,Xplus2(1,:),'g')