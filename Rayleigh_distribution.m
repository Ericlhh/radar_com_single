%-----------------------Rayleigh信道衰落系数的计算-------------------------%
clc
clear
close all

sigma_ = [0.5 1 2 3 4];
figure;
for i = 1:5
    sigma = sigma_(i);
    V = 0:0.01:6;
    pdf_rayleigh = (2*V/sigma^2).*exp(-V.^2/sigma^2);
    plot(V, pdf_rayleigh,'-');
    hold on;

    N = 10^5;                                         %测试样本点数
    x = sqrt(sigma^2/2)*randn(2,N)+0;  %产生均值为0,方差为sigma^2/2的高斯随机变量
    xx = x(1,:) + sqrt(-1)*x(2,:);
    R = abs(xx);                                      %信号包络(理论上服从瑞利分布)

    k = 0:0.01:6;                                     %统计区间
    dk = 0.01;                                        %组距
    p = [];
    for i = 1:(length(k)-1)
        num = length(find(R >= k(i) & R < k(i+1)));   %找到R里面介于k(i)和k(i+1)的元素的个数
        p = [p num/length(R)];                        %得到R中的值介于k(i)和k(i+1)的概率
    end
    pdf = p/dk;                                       %概率除以组距即可得到概率密度
    k_new = 0.005:0.01:(6-0.005);                     %取一组的中值表征这个区间
    plot(k_new, pdf,'-.');
    grid on;
    set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
    xlabel('R');
    ylabel('PDF');
    title('Rayleigh衰落的概率密度函数');
    legend('理论值sigma=0.5','仿真值sigma=0.5','理论值sigma=1','仿真值sigma=1',...
        '理论值sigma=2','仿真值sigma=2','理论值sigma=3','仿真值sigma=3',...
        '理论值sigma=4','仿真值sigma=4');
end
% ————————————————
% 版权声明：本文为CSDN博主「凝望，划过星空.scut」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
% 原文链接：https://blog.csdn.net/weixin_44586473/article/details/114484088