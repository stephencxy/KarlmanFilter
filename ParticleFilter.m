t = 0:0.01:1;
L = length(t);
x = zeros(1,L);
y = zeros(1,L);

x(1) = 0.1;
y(1) = 0.01^2;

for i = 2:L%生成数据
    x(i) = sin(x(i-1))+5*x(i-1)/(x(i-1)^2+0.5);%状态方程
    y(i) = x(i)^2 + normrnd(0,0.01);%观测方程
end

n = 100;%粒子数
%开辟空间
xold = zeros(1,n);
xnew = zeros(1,n);
xplus = zeros(1,n);
w = zeros(1,n);
%初始化
for i = 1:n
    xold(i) = 0.1;
    w(i) = 1/n; 
end

for i = 2:L
    for j = 1:n
        xold(j) = sin(xold(j))+5*xold(j)/(xold(j)^2+0.5)+normrnd(0,0.1);%更新粒子位置
    end
    for j = 1:n
        w(j) = exp(-((y(i)-xold(j)^2)^2/(2*0.001)));%更新权重
    end
    w = w/sum(w);
    %重采样
    c = zeros(1,n);
    c(1) = w(1);
    for j = 2:n
        c(j) = c(j-1)+w(j);
    end

    for j =1 :n
        a=unifrnd(0,1);
        for k = 1:n
            if(a<c(k))
                xnew(j) = xold(k);
                break;
            end
        end
    end
    %将重采样后的粒子赋予粒子
    xold = xnew;
    for j= 1:n
        w(j)=1/n;
    end
    %获取当前估计值（通过期望）
    xplus(i) = sum(xnew)/n;
end

plot(t,x,t,xplus)