%QPSK复数基带误码率仿真
clear all;close all;     
Br=5e6;                 %雷达带宽
fs=10*Br;               %采样频率
Tp=5e-6;                %雷达脉宽
t=0:1/fs:15*Tp+Tp-1/fs;      %采样时间
N_r=length(t);          %采样点数
fc=3e9;                 %载波频率
%% 通信信号 采用QPSK调制
Tc=2*1/fs;                %码元周期
Bc=1/Tc;
N=16*fix(Tp/Tc);
bitstream=randi([0,1],1,2*N);%随机产生的比特数0、1
bitstream=2*bitstream-1;%单极性变为双极性（0到-1；1到1）
I=[];Q=[];
%奇数进I路,偶数进Q路
for i=1:2*N
    if mod(i,2)~=0
        I=[I,bitstream(i)];
    else
        Q=[Q,bitstream(i)];
    end
end
bit_t=0:1/fs:Tc-1/fs;%定义一个码元的时间轴
carrier=[];
for i=1:N
    carrier=[carrier,(I(i)+j*Q(i))*exp(j*2*pi*fc*(bit_t+(i-1)*Tc))];%Q路载波信号
end
%传输信号
QPSK_signal=carrier.*exp(-j*2*pi*fc*t);
snr=0:1:10;
for k=1:length(snr)
    QPSK_receive=awgn(QPSK_signal,snr(k),'measured');%awgn()添加噪声
    I_recover=[];
    Q_recover=[];
    %% 简单信号处理
    for i=1:N
        if real(QPSK_receive(2*i)+QPSK_receive(2*i-1))>0 %积分器求和，大于0为1，否则为-1
            I_recover=[I_recover,1];
        else
            I_recover=[I_recover,-1];
        end
        if imag(QPSK_receive(2*i)+QPSK_receive(2*i-1))>0 %积分器求和，大于0为1，否则为-1
            Q_recover=[Q_recover,1];
        else
            Q_recover=[Q_recover,-1];
        end
    end
    recover_signal=I_recover+j*Q_recover;
    send_signal=I+j*Q;
    err=sum(abs(I_recover-I)/2+abs(Q_recover-Q)/2);
    ber(k)=err/(2*N);
end
figure;
plot(snr,ber);hold on;
plot(ebno0,ber0_thero);