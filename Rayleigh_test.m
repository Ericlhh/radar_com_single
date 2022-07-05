%Rayleigh_test
%该matlab程序主要为是为了检测一个瑞利分布的随机数再乘以一个相位均匀分布的余弦值后，其分布是否是正态分布
clear all;
close all;
N = 10^5;                                         %测试样本点数
sigma=0.5;
r=raylrnd(sigma,1,N);   %生成N个瑞利分布随机数
fai=2*pi*rand(1,N);   %生成N个（0，2pi）均匀分布的相位
a=cos(fai);
b=a.*r;    %b即为相乘以后生成的随机数
k=-2:0.01:2;
dk = 0.01;                                        %组距
p = [];
for i = 1:(length(k)-1)
    num = length(find(b >= k(i) & b < k(i+1)));   %找到b里面介于k(i)和k(i+1)的元素的个数
    p = [p num/length(b)];                        %得到b中的值介于k(i)和k(i+1)的概率
end
pdf = p/dk;
figure;
k_new = -1.995:0.01:(2-0.005);                     %取一组的中值表征这个区间
plot(k_new, pdf,'-.');
hold on;
y=normpdf(k,0,sigma);
plot(k, y,'-.');
grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
xlabel('b');
ylabel('PDF');
legend('仿真值sigma=0.5','理论值sigma=0.5的标准正太分布');
