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
N_target=1;             %目标个数
Rmax=c/2*15*Tp;                             %目标最大距离（本来应该是1/2*c*Tr，但是采样时间限制了不可能那么大）
R_t=Rmax*abs(rand(1,N_target));             %目标的距离（这样以来目标的距离一定是小于最大距离的）
RCS_t=1*(exp(i*2*pi*rand(1,N_target)));    %目标RCS，幅度为1，相位在（0,2pi）之间随机分布
Vmax=lamda*PRF/2;                           %目标最大速度，最大测速范围满足在第一盲速之内
v=Vmax*((1+rand(1,N_target))/2);            %目标速度（这样以来目标的速度一定是小于第一盲速的），每一个目标都有一个自己的速度，对应一个矩阵
Na=8;%接收天线数目
spacing = lamda/2;
RxAntLoc = spacing*[0:Na-1];
angles=fix(rand(1,N_target)*180)-90;
angles=50;
angles=[angles;zeros(1,length(angles))];
B_L = zeros(Na,size(angles,2));%观测向量的设置
for i = 1:size(angles,2)
    B_L(:,i) =steervec(RxAntLoc/lamda,angles(:,i));
end
%% 生成目标矩阵
sr=zeros(N_mc*Na,N_r); %N_mc 脉冲个数   N_r 采样点数
for i=1:N_mc
    ta=(i-1)*Tr;
    sri=0;%每一次从内层for循环出来之后，我们认为上一个脉冲的回波不会干扰到下一个脉冲的回波
    %%内层for循环，一个目标一个目标来研究，对应每一个回波脉冲是由每一个目标回波之和组成
    for k=1:N_target
        tao=2*(R_t(k)-v(k).*(ta))/c;
        srj=RCS_t(k).*rectpuls(t-tao-Tp/2,Tp).*exp(-1j*2*pi*fc*tao+1j*pi*Kr.*(t-tao-Tp/2).^2);
        srj=B_L(:,k)*srj;
        sri=sri+srj;
    end
    %%外层for循环，不同的脉冲，对应的ta是不同值，再代入来计算回波
    sr((i-1)*Na+1:i*Na,:)=sri;
end
%生成通信信号
K=6;     %发射用户数目
L=Na;     %接收天线数
AoD1 = [-30;0]; AoD2 = [-10;0]; AoD3 = [8;0]; AoD4 = [25;0];AoD5 = [30;0];AoD6 = [70;0];
Uangles = [AoD1, AoD2, AoD3, AoD4, AoD5, AoD6];
H=zeros(L,K);
for i=1:K
    H(:,i) =steervec(RxAntLoc/lamda,Uangles(:,i));
end
Es=1/K;    %初始化能量
SNR_dB=10;
N0_dB=10*log10(K*Es)-SNR_dB;  %信噪比用对数形式表示时，SNR_dB=S_dB-N_dB,注意这里的信号功率要乘以K，因为有K道发射信号
N0=10.^(N0_dB/10);
N_sym0=length(t);                   %每根天线处理的QPSK符号数
% H=sqrt(0.5)*(randn(L,K)+1i*randn(L,K)); %%产生L、K路QPSK信号，H信道响应,此为瑞利衰落信道的H
Rece_data=zeros(L*N_mc,length(t));
for i=1:N_mc
%QPSK
Dt=round(rand(K,N_sym0)*3); 
Dt=2*Dt+1;
modDt=sqrt(Es)*exp(1i*Dt/4*pi);
% %BPSK
% Dt=round(rand(K,N_sym0)*1);
% modDt=2*Dt-1;
HS=H*modDt;                                  %接收端通信信号
Noise=sqrt(0.5)*(randn(L,N_sym0)+1i*randn(L,N_sym0));  %L路接收端的噪声
RxDt=HS+sqrt(N0)*Noise;   %接收信号（考虑L路噪声影响）

Rece_data((i-1)*Na+1:i*Na,:)=RxDt+sr((i-1)*Na+1:i*Na,:);
end
res1=0;
res2=0;
res3=0;
Rece_data_1=Rece_data(1:Na,:);
for i=1:4000
    tt1=modDt(:,i);
    tt2=RxDt(:,i);
    tt3=Rece_data_1(:,i);
    res1=res1+tt1*tt1';
    res2=res2+tt2*tt2';
    res3=res3+tt3*tt3';
end
Rxcol=res1/4000;
Rxcol2=res2/4000;
Rxcol3=res3/4000;
% Rnn=Rxcol2;%可以利用获得的通信信号加噪声，来获得Rnn
%同时也可以利用发射信号两两正交以及抑制噪声功率的情况下，获得MVDR要利用的矩阵Rnn
Rnn=H*H'/K+N0*eye(L);
angleSpace = linspace(-pi/2, pi/2, 360);
angleSpaceDeg = linspace(-90, 90, 360);
a = zeros(Na, length(angleSpace));
for j = 1:Na
    a(j,:) = exp((-1i * 2 * pi * RxAntLoc(j) / lamda) .* sin(angleSpace));
end
d = zeros(size(angleSpace));
for i = 1:length(angleSpace)
    d(i) = transpose(a(:,i)) * Rxcol3 * conj(a(:,i));
end
plot(angleSpaceDeg, abs(d));
xlabel('Angle Space [-90^\circ,90^\circ]'); 
%% 
%MVDR获得雷达信号，提出通信信号
w=inv(Rnn)*B_L(:,1)/(B_L(:,1)'*inv(Rnn)*B_L(:,1));
for i=1:N_mc
    r_w(i,:)=w'*Rece_data((i-1)*Na+1:i*Na,:);
end
% r_w_1=r_w(1,:);
% % r_w_1=w'*Rece_data_1;
% err=r_w_1-srj(1,:);
% r_w_abs=abs(r_w_1);
% srj_abs=abs(srj(1,:));
%处理后的雷达信号脉冲压缩
st=rectpuls(t-Tp/2,Tp).*exp(1i*pi*Kr*(t-Tp/2).^2);%参考信号 时域 也就是匹配滤波器的时域
stf=conj(fft(st));%匹配滤波器的频域特性
% r_w_stf=ifft(fft(r_w_1).*stf); 
for i=1:N_mc
    r_w_stf(i,:)=ifft(fft(r_w(i,:)).*stf);  %分别对每一行脉冲压缩 频域脉冲压缩  
end
r_dis=Rece_data(1,:);%未进行处理的接收信号
r_dis_stf=ifft(fft(r_dis).*stf);
R=c*t/2;
figure;
plot(t*c/2,abs(r_w_stf(1,:)));
figure;
plot(t*c/2,abs(r_dis_stf(1,:)));
r_w_v=fft(r_w_stf,[],1);
V=linspace(0,PRF,51)*lamda/2;
figure;image(R,V,255*abs(r_w_v)/max(max(abs(r_w_v))))  ;
%% 剔除雷达信号，获得通信信号,我们目前只处理雷达脉冲最后一段的通信信号，前面也可以处理
%首先我们可以根据处理过后的雷达信号得到雷达信号干扰的区域
%MVDR参与的信号能量大约为1/sqrt((B_L(:,1)'*inv(Rnn)*B_L(:,1))),剩下的雷达信号的能量模大概为接收到能量，即处理过后的雷达信号尽可能保留了雷达信号本来的原貌
threshold=3*1/sqrt((B_L(:,1)'*inv(Rnn)*B_L(:,1)));%阈值的选择决定了干扰的区域，这个必须谨慎
r_w_abs=abs(r_w(N_mc,:));
inter_area=find(r_w_abs>threshold);
Rece_inter_area=Rece_data(393:400,inter_area(1):inter_area(end));
Rece_non_area_1=Rece_data(393:400,1:inter_area(1)-1);
Rece_non_area_2=Rece_data(393:400,inter_area(end)+1:end);
%采用MMSE的思想
w_mmse_inter=H'*inv(H*H'+N0*eye(Na)/Es+B_L*B_L'/Es);
w_mmse_nointer=H'*inv(H*H'+N0*eye(Na)/Es);
Rece_com_mmse_inter=w_mmse_inter*Rece_inter_area;
estDt_mmse_inter=(sign(real(Rece_com_mmse_inter))+1i*sign(imag(Rece_com_mmse_inter)));
Rece_com_mmse_non_1=w_mmse_nointer*Rece_non_area_1;
Rece_com_mmse_non_2=w_mmse_nointer*Rece_non_area_2;
Rece_com_mmse=[Rece_com_mmse_non_1,Rece_com_mmse_inter,Rece_com_mmse_non_2];
estDt_mmse=(sign(real(Rece_com_mmse))+1i*sign(imag(Rece_com_mmse)));
modDt_t=(sign(real(modDt))+1i*sign(imag(modDt)));
b_err_mmse=estDt_mmse-modDt_t;
err_mmse=find(b_err_mmse~=0);
%采用ZF的思想
w_zf=pinv(H);
Rece_com_zf=w_zf*Rece_data(393:400,:);
estDt_zf=(sign(real(Rece_com_zf))+1i*sign(imag(Rece_com_zf)));
b_err_zf=estDt_zf-modDt_t;
err_zf=find(b_err_zf~=0);
% Rnn_c=B_L*B_L'+N0*eye(L);
% w_c=inv(Rnn_c)*B_L(:,1)/(B_L(:,1)'*inv(Rnn_c)*B_L(:,1));