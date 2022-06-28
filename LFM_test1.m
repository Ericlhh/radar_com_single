%%线性调频信号
clear all;
close all;
T=10e-6;                                  %p脉冲持续时间10us
B=30e6;                                   %线性调频信号的频带宽度30MHz
K=B/T;                                      %调频斜率
Fs=2*B;Ts=1/Fs;                      %采样频率和采样间隔
N=T/Ts;
t=linspace(-T/2,T/2,N);
St=exp(j*pi*K*t.^2);                    %线性调频信号
subplot(211)
plot(t*1e6,real(St));
xlabel('时间/us');
title('线性调频信号的实部');
grid on;axis tight;
subplot(212)
freq=linspace(-Fs/2,Fs/2,N);
plot(freq*1e-6,fftshift(abs(fft(St))));
xlabel('频率/MHz');
title('线性调频信号的幅频特性');
grid on;axis tight;
