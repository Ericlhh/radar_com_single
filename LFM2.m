%%%%%%%  利用频域处理方法进行脉冲压缩  %%%%%%%
clear all
clc
clf
eps = 1e-10;
B=100e-6;       %信号周期
Fm=1e6;           %带宽
k=Fm/B;          %调频斜率
Ts=1/(5*Fm);        %采样周期
Ns=fix(B/Ts);        %采样点数
Nf=1024;               % fft点数
t=0:Ts:B-Ts; 
y=exp(j*pi*k*t.^2);   %脉冲压缩前的线形调频信号
yfft = fft(y,Nf) ;
h=zeros(1,Ns);
%% %%%%%%%%%%%%%%%%%%%%%%%%%Hamming窗%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for i=1:Ns
    h(i)=conj(y(Ns-i+1));
end
hfft= fft(h,Nf);     % 匹配滤波器的频域响应
lfm =abs(ifft(yfft .*hfft)); %脉冲压缩    
maxval = max (lfm);
lfm = eps + lfm ./ maxval;    % 利用最大值归一化
lfm_db=20*log10(lfm);   %取对数
%%%%%%%%%%%%%% 加窗处理 %%%%%%%
win = hamming(Ns)';
h_w=h.*win;       % 加窗
hfft_w=fft(h_w,Nf);     % 加窗的匹配滤波器的频域响应
lfm_w = abs(ifft(yfft .*hfft_w)); %脉冲压缩 
maxval1 = max(lfm_w);
val=lfm_w ;
lfm_w = eps + lfm_w ./ maxval;    % 利用lfm的最大值归一化
lfm_w1 = eps + val./ maxval1;    % 利用lfm_w的最大值归一化
lfm_w_db=20*log10(lfm_w);   %取对数
lfm_w1_db=20*log10(lfm_w1);   %取对数
%%%%%%%%%%%%%%%%
tt =0:Ts:2*B-Ts;
figure(1)
plot (tt,lfm_db(1:2*Ns),'b')
axis([.2*B 1.8*B -60 0] )
xlabel ('t - seconds ');
ylabel(' db')
title('没有加Hamming窗的脉冲压缩输出')
grid on
figure(2)
plot (tt,lfm_w1_db(1:2*Ns),'r')
axis([.2*B 1.8*B -60 0] )
xlabel ('t - seconds ');
ylabel(' db')
title('加Hamming窗的脉冲压缩输出')
grid on
figure(3)
plot (tt,lfm_db(1:2*Ns),'b',tt,lfm_w_db(1:2*Ns),'r')
axis([.7*B  1.3*B -60 0] )
xlabel ('t - seconds ');
ylabel(' db')
legend('未加Hamming窗','加Hamming窗');
title('脉冲压缩输出对比')
grid on
