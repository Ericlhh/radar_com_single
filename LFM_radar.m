%% LFM 脉冲多普勒雷达
%==================================================================
function LFM_radar(T,B,Rmin,Rmax,R,RCS)
if nargin==0
    T=10e-6;                                          %脉冲宽度10us
    Tr=100e-6;                                          %脉冲重复周期
    B=30e6;                                           %频带宽度30MHz
    Rmin=10000;Rmax=15000;              %测距范围
%     R=[13000,13500];%目标点的位置，每一个目标相对于雷达的斜距
%     RCS=[1,1];                           %雷达截面积，一维数组
    RCS=[1,1,1,1,1,1];
    R=[10500,11000,12000,12008,13000,13002];  %目标点的位置，每一个目标相对于雷达的斜距
    RCS=[1 1 1 1 1 1];   %雷达截面积，一维数组
end
%====================================
%%
C=3e8;                                            %光速
K=B/T;                                             %调频斜率
Rwid=Rmax-Rmin;                           %最大测距长度
Twid=2*Rwid/C;                               %回波窗的长度
Fs=5*B;Ts=1/Fs;                             %采样频率与采样时间
Nwid=ceil(Twid/Ts);                         %采样窗内的采样点数
%==================================================================
%%产生回拨    
t=linspace(2*Rmin/C,2*Rmax/C,Nwid); %回波窗
                                                            %open window when t=2*Rmin/C
                                                            %close window when t=2*Rmax/C                            
M=length(R);                                        %目标的个数                                       
td=ones(M,1)*t-2*R'/C*ones(1,Nwid);
Srt=RCS*(exp(j*pi*K*td.^2).*(abs(td)<T/2));%从点目标来的回波  
%==================================================================
%%数字信号处理  脉冲压缩
Nchirp=ceil(T/Ts);                               %脉冲宽度离散化
Nfft=2^nextpow2(Nwid+Nwid-1);          %方便使用FFT算法，满足2的次方形式 
                                                          
Srw=fft(Srt,Nfft);                                  %回波做FFT
t0=linspace(-T/2,T/2,Nchirp); 
St=exp(j*pi*K*t0.^2);                            %线性调频信号原始信号作为参考信号  
% %%加窗处理
% win=blackman(Nwid)';
% St_w=St.*win';
% %%
Sw=fft(St,Nfft);                                    %参考信号做FFT
Sot=fftshift(ifft(Srw.*conj(Sw)));           %脉冲压缩后的信号
% %考虑另外一种匹配滤波，时域反褶共轭再补零fft
% Sw2 =fft(conj( fliplr(St)),Nfft); %时域匹配滤波为发射信号时间反褶再取共轭
% Sot2=ifft(Srw.*Sw2);
% Z2=abs(Sot2);
% Z2=Z2/max(Z2);
% Z2=20*log10(Z2+1e-6);
% %
%==================================================================
N0=Nfft/2-Nchirp/2;
Z=abs(Sot(N0:N0+Nwid-1));
Z=Z/max(Z);
Z=20*log10(Z+1e-6);
%figure
subplot(211)
plot(t*1e6,real(Srt));axis tight;
xlabel('Time/us');ylabel('幅度')
title('雷达回波没经过脉冲压缩');
subplot(212)
plot(t*C/2,Z)
axis([10000,15000,-60,0]);
xlabel('距离/m');ylabel('幅值/dB')
title('雷达回波经过脉冲压缩');
