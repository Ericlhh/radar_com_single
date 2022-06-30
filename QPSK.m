clear all;close all;
fc=2.06e9;                 %载波频率   太关键了，要保证低采样时，载波频率与采样频率之间不存在倍数关系，否则会出大问题
Br=5e6;                 %雷达带宽
fs=10*Br;               %采样频率,满足带通采样定理
Tc=1e-7;                %码元周期
Bc=1/Tc;
Tp=5e-6;                %雷达脉宽
t=0:1/fs:15*Tp+Tp-1/fs;      %采样时间
N_r=length(t);          %采样点数
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
I_carrier=[];Q_carrier=[];
bit_t=0:1/fs:Tc-1/fs;%定义一个码元的时间轴
for i=1:N
    I_carrier=[I_carrier,I(i)*cos(2*pi*fc*(bit_t+(N-1)*Tc))];%I路载波信号
    Q_carrier=[Q_carrier,Q(i)*cos(2*pi*fc*(bit_t+(N-1)*Tc)+pi/2)];%Q路载波信号
end
%传输信号
QPSK_signal=I_carrier+Q_carrier;%复数只需要给后面加个j即可
%绘图
figure();%产生一个新图
subplot(3,1,1)
plot(t,I_carrier);legend('I signal')%I路信号
subplot(3,1,2)
plot(t,Q_carrier);legend('Q signal')%Q路信号
subplot(3,1,3)
plot(t,QPSK_signal);legend('QPSK signal')%I路、Q路和的信号

snr=0;%信躁比
%接收信号
QPSK_receive=awgn(QPSK_signal,snr,'measured');%awgn()添加噪声
%解调
% for i=1:N
%     I_output=QPSK_receive(1,(i-1)*length(bit_t)+1:i*length(bit_t)).*cos(2*pi*fc*bit_t);
%     if sum(I_output)>0 %积分器求和，大于0为1，否则为-1
%         I_recover(i)=1;
%     else
%         I_recover(i)=-1;
%     end
%      Q_output=QPSK_receive(1,(i-1)*length(bit_t)+1:i*length(bit_t)).*cos(2*pi*fc*bit_t+ pi/2);
%     if sum(Q_output)>0
%         Q_recover(i)=1;
%     else
%         Q_recover(i)=-1;
%     end
% end
II_output=[];
QQ_output=[];
for i=1:N
    II_output=[II_output,QPSK_receive(1,(i-1)*length(bit_t)+1:i*length(bit_t)).*cos(2*pi*fc*(bit_t+(N-1)*Tc))];
    QQ_output=[QQ_output,QPSK_receive(1,(i-1)*length(bit_t)+1:i*length(bit_t)).*cos(2*pi*fc*(bit_t+(N-1)*Tc)+ pi/2)];
end
% err=sum(abs(I-I_recover))/2+sum(abs(Q-Q_recover))/2;
% err_ber=err/(2*N);
figure;
plot(abs(fft(II_output)));
<<<<<<< HEAD
II_fliter=lowpass(II_output,0.85*Bc,fs);
QQ_fliter=lowpass(QQ_output,0.85*Bc,fs);
=======
II_fliter=lowpass(II_output,Bc,fs);
QQ_fliter=lowpass(QQ_output,Bc,fs);
>>>>>>> cb96433325c15d808a6e0e9b471bc8db2345bdf1
figure;
plot(abs(fft(II_fliter)));
figure;
plot(abs(fft(QQ_fliter)));
for i=1:N
    I_middle=sum(II_fliter(1,(i-1)*5+1:i*5));
    if sum(I_middle)>0 %积分器求和，大于0为1，否则为-1
        I_recover(i)=1;
    else
        I_recover(i)=-1;
    end
    Q_middle=sum(QQ_fliter(1,(i-1)*5+1:i*5));
    if sum(Q_middle)>0 %积分器求和，大于0为1，否则为-1
        Q_recover(i)=1;
    else
        Q_recover(i)=-1;
    end
end
err=sum(abs(I-I_recover))/2+sum(abs(Q-Q_recover))/2;
err_ber=err/(2*N);
figure;
plot(abs(fft(QPSK_receive)));