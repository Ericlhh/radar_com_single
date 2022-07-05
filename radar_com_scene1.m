clear all;
close all;
clc;
K=6;     %发射天线数
L=8;     %接收天线数
Es=1;    %初始化白噪声能量
SNR_dB=10;
N0_dB=10*log10(K*Es)-SNR_dB;  %信噪比用对数形式表示时，SNR_dB=S_dB-N_dB,注意这里的信号功率要乘以K，因为有K道发射信号
N0=10.^(N0_dB/10);
N_sym0=128;                   %每根天线处理的QPSK符号数
fc=3e9;      %载波频率
H=sqrt(0.5)*(randn(L,K)+1i*randn(L,K)); %%产生L、K路QPSK信号，H信道响应
% %QPSK
% Dt=round(rand(K,N_sym0)*3); 
% Dt=2*Dt+1;
% modDt=exp(1i*Dt/4*pi);
%BPSK
Dt=round(rand(K,N_sym0)*1);
modDt=2*Dt-1;
HS=H*modDt;                                  %接收端通信信号

%接收端雷达信号
M = 4;
N = 128;%N=N_sym0
epi = 0.02;
FAR_model = zeros(N,M*N);
%Cn = randperm(M)-1
block_sparsity=10;
for n = 0 : N-1
    Cn = floor(rand()*M);
    for q = 0 : N-1
        for p = 0:M-1
            FAR_model(n+1,q*M+p+1) = exp(1i*2*pi*p/M*Cn+1i*2*pi*q/N*n*(1+Cn*epi));
        end
    end
end
sparse_signal = zeros(M,N);
block = randperm(N,block_sparsity);
sparse_signal(:,block) = exp(1i*2*pi*rand(M,block_sparsity));
y_radar= FAR_model * sparse_signal(:);
cvx_begin
variable x(M,N) complex
norm21 = 0;
for i = 1:N
    norm21 = norm21 + norm(x(:,i));
end
minimize(norm21)
subject to
FAR_model * x(:) == y_radar
cvx_end

y_radar_antenna=ones(L,1)*y_radar.';%天线上接收到的雷达信号
y_antenna_r=y_radar_antenna+HS;%天线实际上收到的通信信号加上雷达信号，暂不考虑噪声
y_radar_con=zeros(N,1);%代表即将利用凸优化方法求解接收信号中混杂的雷达信号
for i=1:N
    y_i=y_antenna_r(:,i);
  
    cvx_begin
       variable xx complex 
       middle=pinv(H)*(y_i-ones(L,1)*xx);
       middle=angle(middle);
       %minimize(sum(abs(abs(real(middle))-abs(imag(middle))))) 非凸
       minimize(abs(abs(real(middle2))-abs(imag(middle2))))
%        minimize(sum(abs(imag(middle))))%针对BPSK
    cvx_end
    y_radar_con(i,1)=xx;
end
error_radar=norm(y_radar_con-y_radar);