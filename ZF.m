%ZF
clear all;
close all;
SNR_dB=0:2:20;
K=6;     %发射天线数
L=8;     %接收天线数
Es=1;    %初始化白噪声能量
len_SNR=length(SNR_dB);
N0_dB=10*log10(K*Es)-SNR_dB;  %信噪比用对数形式表示时，SNR_dB=S_dB-N_dB,注意这里的信号功率要乘以K，因为有K道发射信号
N0=10.^(N0_dB/10);
%
count=zeros(1,len_SNR);       %错误接受码元计数
BER=zeros(1,len_SNR);         %误码率
N_block=5000;                 %分块数目
N_sym0=100;                   %每根天线处理的QPSK符号数
N_err=2000;                   %最小错误数
n_init=1;                     %信噪比指针
while (n_init<=len_SNR)&&(count(len_SNR)<N_block)
	H=sqrt(0.5)*(randn(L,K)+1i*randn(L,K));      %%产生L、K路QPSK信号，H信道响应
    Dt=round(rand(K,N_sym0)+1i*rand(K,N_sym0));  %发送K路QPSK信号
    modDt=sqrt(Es/2)*(Dt*2-(1+1i));
    HS=H*modDt;                                  %接收端信号
	Noise=sqrt(0.5)*(randn(L,N_sym0)+1i*randn(L,N_sym0));  %L路接收端的噪声
    for n=n_init:len_SNR
        count(n)=count(n)+1
        n0=N0(n);
        RxDt=HS+sqrt(n0)*Noise;   %接收信号（考虑L路噪声影响）
        if K==L
            W=inv(H);             %inv函数表示矩阵求逆变换，W为加权系数
        else
            W=pinv(H);
       %pinv函数返回矩阵A的伪逆矩阵。如果矩阵A是可逆（非奇异）的，那么pinv(A)与inv(A)的结果是一样的。
       %但如果矩阵A是奇异矩阵，则inv(A)不存在；
       %但pinv(A)仍然存在，并表现出一些与逆矩阵类似的性质。在pinv函数中，A不一定是方阵。
        end;    
        zt=W*RxDt;
        estDt=(sign(real(zt))+1i*sign(imag(zt))+1+1i)/2;
        err=abs(round(Dt-estDt)).^2;  %错误接收码元数
        BER(n)=BER(n)+sum(sum(err));  %计算误码率
    end;
    if mean(BER(n_init))>=N_err       
        n_init=n_init+1;
    end;   
end;
format short e;                       %format short e格式控制指令，表示5字长浮点数           
BER=BER./(2*K*N_sym0*count);
semilogy(SNR_dB, BER, '-*');
strtitle=['ZF for a ',num2str(K),'x',num2str(L),' QPSK System'];   %strtitle表示字符标题
title(strtitle);
xlabel('Rx SNR per antenna (dB)');
ylabel('BER');
