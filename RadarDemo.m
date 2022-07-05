clear all;
clc;
L=100;
SNR=20;
target_num=3;
t_work=500;

%% 生成信号
t=0.001:0.001:0.001*L;
%x=2*randi([0,1],1,L)-1;
x=chirp(0.001:0.001:0.001*L,0,0.1,100);
y=zeros(1,t_work);
max_delay=t_work-L;
delay=sort(randi([0,max_delay],1,target_num));
amp=unifrnd(1,3,1,target_num);
for m=1:target_num
    y(delay(m)+1:delay(m)+L)=y(delay(m)+1:delay(m)+L)+amp(m)*x;
end
y=awgn(y,SNR,'measured');
figure;
stem(x);
figure;
stem(y);

h_ans=zeros(max_delay+1,1);
for m=1:target_num
    h_ans(delay(m)+1)=h_ans(delay(m)+1)+amp(m);
end
figure;
stem(h_ans);
title('Answer');

%% 建立方程
X=zeros(t_work,max_delay+1);
for m=1:max_delay+1
    X(m:m+L-1,m)=x;
end

%% 匹配滤波
h_MF=X.'*y'/sum(abs(x).^2);
figure;
stem(abs(h_MF));
title('Match Filter');

%% 最小二乘
h_LS=(X.'*X)^(-1)*X.'*y';
figure;
stem(abs(h_LS));
title('Least Square');

%% 压缩感知
E=1/t_work*sum(abs(y).^2);
sigma=sqrt(E/(1+10^(SNR/10)));
cvx_begin
    variable h_CS(max_delay+1) complex
    minimize(norm(h_CS,1))
    subject to
        norm(X*h_CS-y')<=2*sqrt(t_work)*sigma;
        %X*h_CS==y';
cvx_end
figure;
stem(abs(h_CS));
title('Compressed Sensing');

%% OMP
h_OMP=zeros(max_delay+1,1);
k=3;
r=y;
index=[];
for m=1:k
    [~,I]=max(abs(X'*r'));
    index=[index,I];
    r=y-(X(:,index)*(X(:,index)'*X(:,index))^(-1)*X(:,index)'*y')';
end
h_OMP_index=(X(:,index)'*X(:,index))^(-1)*X(:,index)'*y';
h_OMP(index)=h_OMP_index;
figure;
stem(abs(h_OMP));
title('OMP');

%% lasso
h_lasso=lasso(X,y','Lambda',0.01);
figure;
stem(abs(h_lasso));
title('lasso');
