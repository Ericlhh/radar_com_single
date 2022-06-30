%LFM测速测距代码
clear all;clc;close all;
fc=3e9;                 %载波频率
PRF=2000;       
Br=5e6;                 %带宽
fs=10*Br;               %采样频率
Tp=5e-6;                %脉宽
Kr=Br/Tp;               %频率变化率
c=3e8;                  %光速
lamda=c/fc;             %波长
Tr=1/PRF;               %脉冲重复周期
N_mc=1.5/60*PRF;        %脉冲个数
t=0:1/fs:15*Tp+Tp;      %采样时间
N_r=length(t);          %采样点数
N_target=5;             %目标个数
Rmax=c/2*15*Tp;                             %目标最大距离（本来应该是1/2*c*Tr，但是采样时间限制了不可能那么大）
R_t=Rmax*abs(rand(1,N_target));             %目标的距离（这样以来目标的距离一定是小于最大距离的）
RCS_t=10*(exp(i*2*pi*rand(1,N_target)));    %目标RCS，幅度为10，相位在（0,2pi）之间随机分布
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
        srj=RCS_t(k).*rectpuls(t-tao-Tp/2,Tp).*exp(-1j*2*pi*fc*tao+1j*pi*Kr.*(t-tao-Tp/2).^2);
        sri=sri+srj;
    end
    %%外层for循环，不同的脉冲，对应的ta是不同值，再代入来计算回波
    sr(i,:)=sri;
end

%% 脉冲压缩前的回波
tm=(1:N_mc)/PRF;
R=c*t/2;
figure(1);
image(R,tm,255*abs(sr)/max(max(abs(sr))));
figure(2);
plot(t*c/2,abs(sr(1,:)));     %画图我们只反映了第一个脉冲的回波情况

%% 脉冲压缩
st=rectpuls(t-Tp/2,Tp).*exp(1i*pi*Kr*(t-Tp/2).^2);%参考信号 时域 也就是匹配滤波器的时域
stf=conj(fft(st));%匹配滤波器的频域特性
for i=1:N_mc
    sr(i,:)=ifft(fft(sr(i,:)).*stf);  %分别对每一行脉冲压缩 频域脉冲压缩          
end
figure(3);
image(R,tm,255*abs(sr)/max(max(abs(sr))))                
figure;
plot(t*c/2,abs(sr(1,:)))                       
             
sr=fft(sr,[],1);
V=linspace(0,PRF,51)*lamda/2;
figure;image(R,V,255*abs(sr)/max(max(abs(sr))))  
