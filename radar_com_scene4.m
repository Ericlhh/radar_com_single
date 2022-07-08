%radar_com_scene4
%LFM测速测距代码
clear all;close all;  
fc=2.06e9;                 %载波频率    
Br=5e6;                 %带宽
fs=10*Br;               %采样频率
Tp=5e-6;                %脉宽
Kr=Br/Tp;               %频率变化率
c=3e8;                  %光速
lamda=c/fc;             %波长
N_mc=1;
PRF=2000;
Tr=1/PRF;
t=0:1/fs:15*Tp+Tp-1/fs;      %采样时间
N_r=length(t);          %采样点数
N_target=1;             %目标个数
Rmax=c/2*15*Tp;                             %目标最大距离（本来应该是1/2*c*Tr，但是采样时间限制了不可能那么大）
% R_t=Rmax*abs(rand(1,N_target));             %目标的距离（这样以来目标的距离一定是小于最大距离的）
R_t=1/2*Rmax;      %固定目标距离
RCS_t=1*(exp(i*2*pi*rand(1,N_target)));    %目标RCS，幅度为1，相位在（0,2pi）之间随机分布
Vmax=lamda*PRF/2;                           %目标最大速度，最大测速范围满足在第一盲速之内
v=Vmax*((1+rand(1,N_target))/2);            %目标速度（这样以来目标的速度一定是小于第一盲速的），每一个目标都有一个自己的速度，对应一个矩阵
radar_power=1;
%% 生成目标矩阵
sr=zeros(N_mc,N_r); %N_mc 脉冲个数   N_r 采样点数
for i=1:N_mc
    ta=(i-1)*Tr;
    sri=0;%每一次从内层for循环出来之后，我们认为上一个脉冲的回波不会干扰到下一个脉冲的回波
    %%内层for循环，一个目标一个目标来研究，对应每一个回波脉冲是由每一个目标回波之和组成
    for k=1:N_target
        tao=2*(R_t(k)-v(k).*(ta))/c;
        srj=RCS_t(k).*rectpuls(t-tao-Tp/2,Tp).*exp(-1j*2*pi*fc*tao+1j*pi*Kr.*(t-tao-Tp/2).^2);
        sri=sri+srj;
    end
    %%外层for循环，不同的脉冲，对应的ta是不同值，再代入来计算回波
    sr(i,:)=sri;
end

                    
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
com_power=1;%通信功率
for i=1:N
    carrier=[carrier,sqrt(1/2)*com_power*(I(i)+j*Q(i))*exp(j*2*pi*fc*(bit_t+(i-1)*Tc))];%Q路载波信号
end
%传输信号
QPSK_signal=carrier.*exp(-j*2*pi*fc*t);
snr=7;
P_noise=10^(-snr/10)*com_power;
sinr=10*log10(com_power/(1+P_noise));
% QPSK_receive=QPSK_signal;
QPSK_receive=awgn(QPSK_signal,snr,'measured');%awgn()添加噪声
receive_signal=QPSK_receive+sr;
%% 未处理信号脉冲压缩
st=rectpuls(t-Tp/2,Tp).*exp(1i*pi*Kr*(t-Tp/2).^2);%参考信号 时域 也就是匹配滤波器的时域
stf=conj(fft(st));%匹配滤波器的频域特性
for i=1:N_mc
    signal_yasuo(i,:)=ifft(fft(receive_signal(i,:)).*stf);  %分别对每一行脉冲压缩 频域脉冲压缩          
end               
figure;
plot(t*c/2,abs(signal_yasuo(1,:)))    
figure;
plot(abs(fft(QPSK_signal)));
radar_reflect=radar_power*abs(RCS_t);
if(com_power>radar_reflect) %%如果通信信号功率较强
    I_recover=[];
    II_recover=[]
    Q_recover=[];
    QQ_recover=[];
    %% 简单信号处理
    for i=1:N
        if real(receive_signal(2*i)+receive_signal(2*i-1))>0 %积分器求和，大于0为1，否则为-1
            I_recover=[I_recover,1];
            II_recover=[II_recover,1,1];
        else
            I_recover=[I_recover,-1];
            II_recover=[II_recover,-1,-1];
        end
        if imag(receive_signal(2*i)+receive_signal(2*i-1))>0 %积分器求和，大于0为1，否则为-1
            Q_recover=[Q_recover,1];
            QQ_recover=[QQ_recover,1,1];
        else
            Q_recover=[Q_recover,-1];
            QQ_recover=[QQ_recover,-1,-1];
        end
    end
    recover_signal=I_recover+j*Q_recover;
    send_signal=I+j*Q;
    err_all=sum(abs(I_recover-I)/2+abs(Q_recover-Q)/2);
    err_radar=sum(abs(I_recover(1,938:1062)-I(1,938:1062))/2+abs(Q_recover(1,938:1062)-Q(1,938:1062))/2);
    ber=err_radar/(250);
    %雷达信号影响的信号区域在【1876-2125】
    receive_angle=angle(receive_signal);
    a=abs(receive_signal).*(cos(pi/4-receive_angle));
    b=abs(receive_signal).*(sin(pi/4-receive_angle));
    
    signal_com_recover=com_power*sqrt(1/2)*(II_recover+j*QQ_recover);
    radar=receive_signal-signal_com_recover;
    for i=1:N_mc
        radar_yasuo(i,:)=ifft(fft(radar(i,:)).*stf);  %分别对每一行脉冲压缩 频域脉冲压缩
    end
    figure;
    plot(t*c/2,abs(radar_yasuo(1,:)));
    else
    %% 当雷达信号功率比较大的时候，先去剪掉雷达信号，这里已经知道了目标的距离
    radar_decect=RCS_t(1).*rectpuls(t-tao-Tp/2,Tp).*exp(-1j*2*pi*fc*tao+1j*pi*Kr.*(t-tao-Tp/2).^2);
    signal_eli_radar=receive_signal-radar_decect;
    I_eli_recover=[];
    II_eli_recover=[];
    Q_eli_recover=[];
    QQ_eli_recover=[];
    for i=1:N
        if real(signal_eli_radar(2*i)+signal_eli_radar(2*i-1))>0 %积分器求和，大于0为1，否则为-1
            I_eli_recover=[I_eli_recover,1];
            II_eli_recover=[II_eli_recover,1,1];
        else
            I_eli_recover=[I_eli_recover,-1];
            II_eli_recover=[II_eli_recover,-1,-1];
        end
        if imag(signal_eli_radar(2*i)+signal_eli_radar(2*i-1))>0 %积分器求和，大于0为1，否则为-1
            Q_eli_recover=[Q_eli_recover,1];
            QQ_eli_recover=[QQ_eli_recover,1,1];
        else
            Q_eli_recover=[Q_eli_recover,-1];
            QQ_eli_recover=[QQ_eli_recover,-1,-1];
        end
    end
    recover_eli_signal=I_eli_recover+j*Q_eli_recover;
    err_all=sum(abs(I_eli_recover-I)/2+abs(Q_eli_recover-Q)/2);
    err_radar=sum(abs(I_eli_recover(1,938:1062)-I(1,938:1062))/2+abs(Q_eli_recover(1,938:1062)-Q(1,938:1062))/2);
    ber=err_radar/(250);
    signal_com_recover=com_power*sqrt(1/2)*(II_eli_recover+j*QQ_eli_recover);
    radar=receive_signal-signal_com_recover;
    for i=1:N_mc
        radar_yasuo(i,:)=ifft(fft(radar(i,:)).*stf);  %分别对每一行脉冲压缩 频域脉冲压缩
    end
    figure;
    plot(t*c/2,abs(radar_yasuo(1,:)));
end
%% 仅雷达信号的脉冲压缩图
snr_radar=10*log10(radar_reflect/(P_noise));
radar_noise=awgn(sr,snr_radar,'measured');%awgn()添加噪声
for i=1:N_mc
        radar_only_yasuo(i,:)=ifft(fft(radar_noise(i,:)).*stf);  %分别对每一行脉冲压缩 频域脉冲压缩
end
figure;
plot(t*c/2,abs(radar_only_yasuo(1,:)));
title("pulse compression(only radar signal)");
xlabel("distance(m)");