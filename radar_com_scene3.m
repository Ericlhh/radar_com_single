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
N_target=5;             %目标个数
Rmax=c/2*15*Tp;                             %目标最大距离（本来应该是1/2*c*Tr，但是采样时间限制了不可能那么大）
R_t=Rmax*abs(rand(1,N_target));             %目标的距离（这样以来目标的距离一定是小于最大距离的）
RCS_t=ones(1,N_target);    %目标RCS，幅度为10，相位在（0,2pi）之间随机分布
Vmax=lamda*PRF/2;                           %目标最大速度，最大测速范围满足在第一盲速之内
v=Vmax*((1+rand(1,N_target))/2);            %目标速度（这样以来目标的速度一定是小于第一盲速的），每一个目标都有一个自己的速度，对应一个矩阵
%% 生成目标矩阵
sr=zeros(N_mc,N_r); %N_mc 脉冲个数   N_r 采样点数
for i=1:N_mc
    ta=(i-1)*Tr;
    sri=0;%每一次从内层for循环出来之后，我们认为上一个脉冲的回波不会干扰到下一个脉冲的回波
    %%内层for循环，一个目标一个目标来研究，对应每一个回波脉冲是由每一个目标回波之和组成
    for k=1:N_target
        tao=2*(R_t(k)-v(k).*(ta))/c;
        srj=RCS_t(k).*rectpuls(t-tao-Tp/2,Tp).*cos(2*pi*fc*(t-tao)+pi*Kr.*(t-tao-Tp/2).^2);
        sri=sri+srj;
    end
    %%外层for循环，不同的脉冲，对应的ta是不同值，再代入来计算回波
    sr(i,:)=sri;
end

                    
%% 通信信号 采用QPSK调制
Tc=1e-7;                %码元周期
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
    carrier=[carrier,(I(i)+j*Q(i))*exp(j*2*pi*fc*(bit_t+(N-1)*Tc))];%Q路载波信号
end
%传输信号
QPSK_signal=real(carrier);
%% 联合雷达信号通信信号处理
receive_signal=QPSK_signal+sr;
figure;
plot(abs(fft(sr)));
figure;
plot(abs(fft(receive_signal)));
figure;
plot(abs(fft(QPSK_signal)));
rece_signal_down=receive_signal.*(cos(2*pi*fc*t));
figure;
plot(abs(fft(rece_signal_down)));
rece_fliter_signal=lowpass(rece_signal_down,Bc,fs);