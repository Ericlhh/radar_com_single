%多天线估计到达角度，利用压缩感知方法,
clear all;
close all;
clc
Nr=8;%接收天线数目
fc = 24e9;
c = 3e8;
lamb = c/fc;

spacing = lamb/2;
RxAntLoc = spacing*[0:Nr-1];
angles=-90:1:90;angles=[angles;zeros(1,length(angles))];
B_L = zeros(Nr,length(angles));%观测向量的设置
for i = 1:length(angles)
    B_L(:,i) =steervec(RxAntLoc/lamb,angles(:,i));
end
sparse_signal=zeros(length(angles),1);
num=randperm(181,4);
tgt=num-91;
sparse_signal(num,1)=exp(1i*2*pi*rand(4,1));
y=B_L*sparse_signal;
result=OMP(B_L,[],y,4);
find(result)