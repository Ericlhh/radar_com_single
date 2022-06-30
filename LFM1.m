%% 线性调频与脉冲压缩
%eric
clear,clc,close all
set(0,'defaultfigurecolor','w')
%% Chirp信号参数设置
lhh=10000;
Tr = 1e-6;%时宽
Br = 200e6;%带宽
Fs = 4*Br;%采样率
%% Chirp信号参数导出
Kr = Br/Tr;%调频率
N =  round( Tr / (1/Fs) );%采样点数
t = linspace( -Tr/2 , Tr/2 , N);%在[-Tp/2,Tp/2]选取采样点
%% Chirp信号生成
st = ( abs(t) < Tr/2 ) .* exp( 1j * pi * Kr * t.^2 ); 
f_chirp= Kr * t; %信号频率
phase_chirp = pi * Kr * t.^2;%信号相位
%% 频谱
freq = linspace(-Fs/2,Fs/2,N);%频域采样
Sf = fftshift( fft(st) );
%% 时域匹配滤波
ht = conj( fliplr(st) ); %时域匹配滤波为发射信号时间反褶再取共轭
s1 = conv(st,ht); %线性调频信号经过匹配滤波器后的输出(时域卷积)
N1 = N+N-1 ;%线性卷积后信号长度变为 N1+N2-1
t1 = linspace( -Tr/2 , Tr/2 , N1);
% 时域匹配滤波
figure,plot( t1*1e6 , abs(s1) ),xlabel('t /us'),ylabel('幅度谱'),title('时间反褶取共轭，时域卷积');
%% 频域匹配滤波1 (复制发射脉冲进行时间反褶并取共轭，计算补零DFT)
N2 = 2*N; %循环卷积长度 （N2应当>=N+N-1，其中弃置区位于长度大于N+N-1的部分）
t2 = linspace( -Tr/2 , Tr/2*(N2/N-1) , N2);
Hf2 = fft(ht,N2); %频域匹配滤波器
Sf2 = fft(st,N2);%频域信号
S2 = Sf2 .* Hf2;%频域乘积
s2 = ifft(S2);
% 频域匹配滤波1
figure,plot( t2*1e6, abs(s2) ),xlabel('t /us'),ylabel('幅度谱'),title('时间反褶取共轭，补零FFT，频域乘积，IFFT');
%弃置区再后面，可以通过将N2设置的大一些观察出来
% ————————————————
% 版权声明：本文为CSDN博主「格桑蓝莲」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
% 原文链接：https://blog.csdn.net/weixin_44566643/article/details/107508520
% %% 频域匹配滤波2（复制脉冲补零后进行DFT，对结果取复共轭（无时间反褶））
N3 = 2*N; %循环卷积长度
t3 = linspace( -Tr/2 , Tr/2 , N3);
Hf3 = conj( fft(st,N3) );% 复制脉冲补零后进行DFT，对结果取复共轭
Sf3 = fft(st,N3);
S3 = Sf3 .* Hf3;%频域乘积
s3 = fftshift(ifft(S3));
% 频域匹配滤波2
figure,plot( t3*1e6 , abs(s3) ),xlabel('t /us'),ylabel('幅度谱'),title('复制脉冲补零后FFT，取共轭，频域乘积，IFFT');
% ————————————————
% 版权声明：本文为CSDN博主「格桑蓝莲」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
% 原文链接：https://blog.csdn.net/weixin_44566643/article/details/107508520
% 第二种方法存疑，因为当存在补零的情况是，这个时候对时域做补零求FFT再共轭得到的频谱和原本时域反褶再共轭再补零求FFT得到的频谱不一样